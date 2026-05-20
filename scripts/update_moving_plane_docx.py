#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import io
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET
import re

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SOURCE_DOCX = DOCS / "MironovMV_RK5-81B_moving_plane_update.docx"
OUTPUT_DOCX = DOCS / "MironovMV_RK5-81B_moving_plane_results.docx"
FIGURES = DOCS / "figures_numbered_moving_plane"
MAIN_RESULTS = ROOT / "results" / "main_scale_hyperelastic_reference_triplet_coarse_moving_plane"
COARSE_RESULTS = ROOT / "results" / "main_scale_hyperelastic_reference_triplet_very_coarse_moving_plane"

EMU_PER_CM = 360_000
MAX_WIDTH_EMU = int(15.6 * EMU_PER_CM)
MAX_HEIGHT_EMU = int(10.6 * EMU_PER_CM)

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
REL = "http://schemas.openxmlformats.org/package/2006/relationships"
MC = "http://schemas.openxmlformats.org/markup-compatibility/2006"

for prefix, uri in {
    "w": W,
    "m": M,
    "r": R,
    "wp": WP,
    "a": A,
    "rel": REL,
    "mc": MC,
}.items():
    ET.register_namespace(prefix, uri)


def q(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


def register_namespaces_from_xml(xml_bytes: bytes) -> None:
    for _, namespace in ET.iterparse(io.BytesIO(xml_bytes), events=("start-ns",)):
        prefix, uri = namespace
        if prefix in ("", "xml"):
            continue
        try:
            ET.register_namespace(prefix, uri)
        except ValueError:
            pass


@dataclass(frozen=True)
class CaseStats:
    load_steps: int
    converged_load_steps: int
    nonlinear_iterations: int
    max_penetration: float
    total_time_seconds: float
    contact_force_norm: float
    contact_patch_length: float | None


@dataclass(frozen=True)
class TripletStats:
    penalty: CaseStats
    augmented: CaseStats
    surrogate: CaseStats
    max_pressure: float
    half_angle_rad: float
    half_angle_deg: float
    half_width: float


def text_nodes(element: ET.Element) -> list[ET.Element]:
    return [node for node in element.iter() if node.tag in (q(W, "t"), q(M, "t"))]


def element_text(element: ET.Element) -> str:
    return "".join((node.text or "") for node in text_nodes(element)).strip()


def set_element_text(element: ET.Element, value: str) -> None:
    nodes = text_nodes(element)
    if not nodes:
        raise RuntimeError(f"Cannot set text on element without text nodes: {element.tag}")
    nodes[0].text = value
    for node in nodes[1:]:
        node.text = ""


def run_text(run: ET.Element) -> str:
    return "".join((node.text or "") for node in text_nodes(run))


def make_text_run(text: str, template_run: ET.Element | None = None) -> ET.Element:
    run = ET.Element(q(W, "r"))
    if template_run is not None:
        run_properties = template_run.find(q(W, "rPr"))
        if run_properties is not None:
            run.append(ET.fromstring(ET.tostring(run_properties, encoding="utf-8")))
    text_node = ET.SubElement(run, q(W, "t"))
    if text.startswith(" ") or text.endswith(" "):
        text_node.set(q("http://www.w3.org/XML/1998/namespace", "space"), "preserve")
    text_node.text = text
    return run


def run_field_char_type(run: ET.Element) -> str | None:
    field_char = run.find(q(W, "fldChar"))
    if field_char is None:
        return None
    return field_char.get(q(W, "fldCharType"))


def collapse_ref_fields(paragraph: ET.Element) -> int:
    children = list(paragraph)
    collapsed = 0
    output: list[ET.Element] = []
    index = 0
    while index < len(children):
        child = children[index]
        if child.tag != q(W, "r") or run_field_char_type(child) != "begin":
            output.append(child)
            index += 1
            continue

        field_runs = [child]
        index += 1
        while index < len(children):
            field_runs.append(children[index])
            if children[index].tag == q(W, "r") and run_field_char_type(children[index]) == "end":
                index += 1
                break
            index += 1

        instruction = "".join(
            (node.text or "") for run in field_runs for node in run.findall(q(W, "instrText"))
        ).strip()
        if not instruction.startswith("REF "):
            output.extend(field_runs)
            continue

        result_started = False
        result_runs: list[ET.Element] = []
        for run in field_runs:
            field_type = run_field_char_type(run)
            if field_type == "separate":
                result_started = True
                continue
            if field_type == "end":
                break
            if result_started and run_text(run):
                result_runs.append(run)

        result_text = "".join(run_text(run) for run in result_runs)
        if result_text:
            output.append(make_text_run(result_text, result_runs[0]))
        collapsed += 1

    if collapsed:
        paragraph[:] = output
    return collapsed


def collapse_all_ref_fields(body_children: list[ET.Element]) -> int:
    return sum(
        collapse_ref_fields(paragraph)
        for paragraph in body_children
        if paragraph.tag == q(W, "p")
    )


FIGURE_RENUMBERING = {
    "2.6": "2.5",
    "2.9": "2.6",
    "2.8": "2.7",
    "4.2": "4.1",
    "4.3": "4.2",
    "4.4": "4.3",
    "4.5": "4.4",
    "4.6": "4.5",
    "4.7": "4.6",
    "5.4": "5.3",
    "5.5": "5.4",
    "5.6": "5.5",
    "5.7": "5.6",
    "5.8": "5.7",
    "5.10": "5.8",
    "5.11": "5.9",
    "5.12": "5.10",
    "5.13": "5.11",
    "5.14": "5.12",
    "5.15": "5.13",
    "5.16": "5.14",
    "5.17": "5.15",
    "5.18": "5.16",
    "5.19": "5.17",
    "5.20": "5.18",
}


FIGURE_NUMBER_PATTERN = re.compile(
    r"(?<![\d.])("
    + "|".join(re.escape(number) for number in sorted(FIGURE_RENUMBERING, key=len, reverse=True))
    + r")(?!\d)"
)


def renumber_figure_text(text: str) -> str:
    return FIGURE_NUMBER_PATTERN.sub(lambda match: FIGURE_RENUMBERING[match.group(1)], text)


def renumber_figures_and_references(body_children: list[ET.Element]) -> int:
    changed = 0
    figure_words = ("рис.", "рис ", "рисунок", "Рис.", "Рис ", "Рисунок")
    for child in body_children:
        if child.tag != q(W, "p"):
            continue
        current = element_text(child)
        if not any(word in current for word in figure_words):
            continue
        updated = renumber_figure_text(current)
        updated = updated.replace("рис 2.6", "рис. 2.6")
        if updated != current:
            set_element_text(child, updated)
            changed += 1
    return changed


def table_rows(table: ET.Element) -> list[list[ET.Element]]:
    return [
        row.findall(q(W, "tc"))
        for row in table.findall(q(W, "tr"))
    ]


def find_table_after_caption(body_children: list[ET.Element], caption: str) -> ET.Element:
    for index, child in enumerate(body_children):
        if child.tag == q(W, "p") and element_text(child).startswith(caption):
            for next_child in body_children[index + 1 :]:
                if next_child.tag == q(W, "tbl"):
                    return next_child
                if element_text(next_child).startswith("Таблица "):
                    break
    raise RuntimeError(f"Table not found after caption: {caption}")


def find_relationship_target(rels_root: ET.Element, rel_id: str) -> str:
    for rel in rels_root.findall(q(REL, "Relationship")):
        if rel.get("Id") == rel_id:
            target = rel.get("Target")
            if not target:
                raise RuntimeError(f"Image relationship {rel_id} has no target")
            return target
    raise RuntimeError(f"Image relationship not found: {rel_id}")


def image_size_emu(path: Path) -> tuple[int, int]:
    with Image.open(path) as image:
        width, height = image.size
    scale = min(MAX_WIDTH_EMU / width, MAX_HEIGHT_EMU / height)
    return int(width * scale), int(height * scale)


def resize_drawing(paragraph: ET.Element, image_path: Path) -> None:
    cx, cy = image_size_emu(image_path)
    for extent in paragraph.findall(f".//{q(WP, 'extent')}"):
        extent.set("cx", str(cx))
        extent.set("cy", str(cy))
    for extent in paragraph.findall(f".//{q(A, 'ext')}"):
        extent.set("cx", str(cx))
        extent.set("cy", str(cy))


def blip_rel_id(paragraph: ET.Element) -> str | None:
    blip = paragraph.find(f".//{q(A, 'blip')}")
    if blip is None:
        return None
    return blip.get(q(R, "embed"))


def read_metrics_case(results_root: Path, case_name: str) -> dict:
    with (results_root / case_name / "metrics.json").open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_triplet(results_root: Path) -> TripletStats:
    summary = results_root / "triplet_summary.csv"
    if not summary.exists():
        raise FileNotFoundError(summary)

    case_rows: dict[str, dict[str, str]] = {}
    parameters: dict[str, float] = {}
    mode = "cases"
    header: list[str] | None = None

    with summary.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle):
            if not row:
                continue
            if row[0] == "case":
                header = row
                mode = "cases"
                continue
            if row[0] == "parameter":
                mode = "parameters"
                continue
            if mode == "cases":
                if header is None:
                    raise RuntimeError(f"Bad summary header in {summary}")
                case_rows[row[0]] = dict(zip(header, row))
            else:
                parameters[row[0]] = float(row[1])

    def case_stats(case_name: str) -> CaseStats:
        row = case_rows[case_name]
        patch_length: float | None = None
        metrics_path = results_root / case_name / "metrics.json"
        if metrics_path.exists():
            metrics = read_metrics_case(results_root, case_name)
            contact = metrics.get("contact", {})
            if "contact_patch_length" in contact:
                patch_length = float(contact["contact_patch_length"])
        return CaseStats(
            load_steps=int(row["load_steps"]),
            converged_load_steps=int(row["converged_load_steps"]),
            nonlinear_iterations=int(row["nonlinear_iterations"]),
            max_penetration=float(row["max_penetration"]),
            total_time_seconds=float(row["total_time_seconds"]),
            contact_force_norm=float(row["contact_force_norm"]),
            contact_patch_length=patch_length,
        )

    return TripletStats(
        penalty=case_stats("penalty_contact"),
        augmented=case_stats("augmented_lagrangian_contact"),
        surrogate=case_stats("no_contact_surrogate"),
        max_pressure=parameters["surrogate_max_pressure"],
        half_angle_rad=parameters["surrogate_contact_half_angle_rad"],
        half_angle_deg=parameters["surrogate_contact_half_angle_deg"],
        half_width=parameters["surrogate_contact_half_width"],
    )


def fmt_k_n(force_n: float) -> str:
    return f"{force_n / 1000.0:.6f}"


def fmt_time(seconds: float) -> str:
    return f"{seconds:.2f}"


def fmt_len(value: float) -> str:
    return f"{value:.3f}"


def fmt_width(value: float) -> str:
    return f"{value:.4f}"


def fmt_pressure(value: float) -> str:
    return f"{value:.5f}"


def fmt_angle(value: float) -> str:
    return f"{value:.4f}"


def fmt_sci(value: float, digits: int = 3) -> str:
    if value == 0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / (10 ** exponent)
    return f"{mantissa:.{digits}f}·10^{{{exponent}}}"


def fmt_ratio(value: float, digits: int = 2) -> str:
    return f"{value:.{digits}f}"


def patch_tables(body_children: list[ET.Element], main: TripletStats, coarse: TripletStats) -> None:
    load_label = "плоскость 10 мм, зазор = 2 мм"

    material_table = body_children[336]
    if material_table.tag != q(W, "tbl"):
        raise RuntimeError("Expected material table at body index 336")
    rows_32 = table_rows(material_table)
    set_element_text(rows_32[2][0], "ν")
    set_element_text(rows_32[3][0], "G")

    table_43 = find_table_after_caption(
        body_children,
        "Таблица 4.3 — Сравнение штрафного метода и метода расширенного Лагранжа",
    )
    rows_43 = table_rows(table_43)
    data_43 = [
        ("грубая сетка", "Штрафной метод", coarse.penalty),
        ("грубая сетка", "метод расширенного Лагранжа", coarse.augmented),
        ("основная расчетная сетка", "Штрафной метод", main.penalty),
        ("основная расчетная сетка", "метод расширенного Лагранжа", main.augmented),
    ]
    for row_index, (mesh, method, stats) in enumerate(data_43, start=1):
        values = [
            mesh,
            method,
            load_label,
            fmt_k_n(stats.contact_force_norm),
            fmt_sci(stats.max_penetration),
            str(stats.nonlinear_iterations),
            f"{stats.converged_load_steps}/{stats.load_steps}",
            fmt_time(stats.total_time_seconds),
        ]
        for column_index, value in enumerate(values):
            set_element_text(rows_43[row_index][column_index], value)

    table_44 = find_table_after_caption(
        body_children,
        "Таблица 4.4 — Длина пятна контакта для двух контактных методов",
    )
    rows_44 = table_rows(table_44)
    if main.penalty.contact_patch_length is None or main.augmented.contact_patch_length is None:
        raise RuntimeError("Contact patch length is missing from main metrics")
    set_element_text(rows_44[1][1], fmt_len(main.penalty.contact_patch_length))
    set_element_text(rows_44[2][1], fmt_len(main.augmented.contact_patch_length))

    table_51 = find_table_after_caption(
        body_children,
        "Таблица 5.1 — Параметры восстановленной параболической нагрузки",
    )
    rows_51 = table_rows(table_51)
    data_51 = [
        ("грубая сетка", coarse),
        ("основная расчетная сетка", main),
    ]
    for row_index, (mesh, stats) in enumerate(data_51, start=1):
        values = [
            mesh,
            fmt_pressure(stats.max_pressure),
            fmt_angle(stats.half_angle_deg),
            fmt_width(stats.half_width),
            fmt_time(stats.surrogate.total_time_seconds),
        ]
        for column_index, value in enumerate(values):
            set_element_text(rows_51[row_index][column_index], value)


def patch_paragraphs(body_children: list[ET.Element], main: TripletStats) -> None:
    pen = main.penalty
    aug = main.augmented
    sur = main.surrogate
    if pen.contact_patch_length is None or aug.contact_patch_length is None:
        raise RuntimeError("Contact patch length is missing from main metrics")

    pen_k_n = fmt_k_n(pen.contact_force_norm)
    aug_k_n = fmt_k_n(aug.contact_force_norm)
    pen_g = fmt_sci(pen.max_penetration)
    aug_g = fmt_sci(aug.max_penetration)
    pen_l = fmt_len(pen.contact_patch_length)
    aug_l = fmt_len(aug.contact_patch_length)
    delta_r = abs(aug.contact_force_norm - pen.contact_force_norm) / aug.contact_force_norm * 100.0
    g_ratio = pen.max_penetration / aug.max_penetration
    time_ratio = aug.total_time_seconds / pen.total_time_seconds
    length_delta = abs(pen.contact_patch_length - aug.contact_patch_length) / aug.contact_patch_length * 100.0
    sur_pen_ratio = pen.total_time_seconds / sur.total_time_seconds
    sur_aug_ratio = aug.total_time_seconds / sur.total_time_seconds

    replacements = {
        331: (
            "Параметры гиперупругой модели были получены из линейно-упругих характеристик E и ν "
            "по стандартным соотношениям изотропной теории упругости [1, 2]:"
        ),
        332: "G=E/(2(1+ν)),(3.3)",
        333: "K=E/(3(1-2ν)).(3.4)",
        338: (
            "Значение коэффициента Пуассона ν=0.48 соответствует почти несжимаемому материалу. "
            "При этом материал не считается строго несжимаемым, так как объемная деформация допускается, "
            "но сильно штрафуется за счет большого значения K. Такой вариант постановки является компромиссом "
            "между физической осмысленностью описания резиноподобного материала и устойчивостью численного решения."
        ),
        435: (
            "Однако в гиперупругой контактной задаче свойства системы могут существенно ухудшаться. "
            "Во-первых, почти несжимаемый материал с коэффициентом Пуассона ν=0.48 приводит к высокой объемной "
            "жесткости. Во-вторых, контактные ограничения добавляют локально жесткие связи на части границы. "
            "В-третьих, изменение активной контактной зоны может вызывать резкое изменение касательной матрицы "
            "между итерациями."
        ),
        444: (
            "Третье ограничение связано с конечным элементом. Используется четырехузловой билинейный элемент Q4 "
            "в формулировке только через перемещения и с полным интегрированием. Для почти несжимаемого материала "
            "при ν=0.48 такая постановка может приводить к эффекту объемного запирания, то есть к искусственному "
            "завышению жесткости модели [1, 5]. В настоящей работе этот риск учитывался при анализе результатов "
            "и частично компенсировался сгущением сетки, однако специальные методы борьбы с запиранием, например "
            "выборочное пониженное интегрирование или смешанная постановка перемещения–давление, не применялись."
        ),
        611: (
            "Таким образом, основным итогом главы является не выбор одного универсального метода, а установление "
            "их рациональных областей применения: штрафной метод — для устойчивого инженерного расчета, метод "
            "расширенного Лагранжа — для уточнения контактного ограничения и получения базового распределения давления."
        ),
        612: (
            "Полученное методом расширенного Лагранжа распределение давления далее используется как исходная "
            "информация для бесконтактной аппроксимации. Это позволяет проверить более узкую идею: можно ли после "
            "одного явного контактного расчета заменить найденное контактное действие заданной нагрузкой и тем "
            "самым удешевить последующие расчеты напряженно-деформированного состояния."
        ),
        614: (
            "В этой главе контактная постановка не заменяется произвольной внешней нагрузкой. Эквивалентная нагрузка "
            "строится по результатам уже решенной явной контактной задачи, поэтому бесконтактная аппроксимация "
            "рассматривается как зависимая расчетная схема, а не как самостоятельная альтернатива контактному алгоритму."
        ),
        615: (
            "После решения явной контактной задачи возникает инженерный вопрос: можно ли использовать найденное "
            "распределение давления как заранее заданную поверхностную нагрузку и получить приемлемое описание "
            "общего напряженно-деформированного состояния без повторного решения контактной подзадачи. Такая "
            "аппроксимация не определяет ни пятно контакта, ни давление из условий непроникания, но может быть "
            "полезна как вычислительное упрощение для уже известного контактного состояния."
        ),
        616: (
            "В данной главе рассматривается возвращение к бесконтактной постановке, но с нагрузкой, восстановленной "
            "по явному контакту. В качестве источника давления используется расчет методом расширенного Лагранжа, "
            "так как в предыдущей главе было показано, что он существенно точнее выполняет условие непроникания "
            "при практически той же суммарной реакции [3]."
        ),
        617: (
            "Цель главы состоит в оценке применимости такой аппроксимации: показать порядок выигрыша по времени, "
            "сопоставить получаемые поля с явным контактным расчетом и сформулировать ограничения, при которых "
            "восстановленная нагрузка может использоваться в инженерном анализе."
        ),
        619: (
            "В явной контактной задаче нормальное давление является неизвестной величиной, возникающей в результате "
            "решения. Оно определяется условиями непроникания, текущей геометрией деформируемого тела, положением "
            "жесткой плоскости и активной областью контакта [3]. В бесконтактной аппроксимации эта логика меняется: "
            "давление считается уже известным и прикладывается к внешнему контуру шины как заданная поверхностная нагрузка."
        ),
        624: (
            "где fc — контактный вклад, определяемый алгоритмом контакта, а fp — заданная эквивалентная нагрузка, "
            "восстановленная по контактному давлению. Главное различие состоит в том, что fc пересчитывается по "
            "текущему зазору и активной контактной зоне, тогда как fp фиксируется до начала бесконтактного расчета."
        ),
        625: (
            "Схема перехода от явного контакта к эквивалентной бесконтактной нагрузке показана на рис. 5.1. "
            "На этой схеме контактный вклад, зависящий от зазора и активной области, заменяется заранее заданной "
            "поверхностной нагрузкой на найденном участке внешнего контура."
        ),
        556: (
            "Расчеты выполнялись для двух сеток с локальным сгущением в зоне контакта. "
            "Для обеих сеток задавалось одинаковое нагружение: начальный зазор 2 мм "
            "и вертикальное перемещение жесткой опорной плоскости на 10 мм при закреплении "
            "внутреннего контура. Поэтому полученные различия характеризуют влияние контактного "
            "метода и, в меньшей степени, плотности сетки."
        ),
        557: "Сводные результаты приведены в таблице 4.3.",
        561: (
            "Из таблицы видно, что суммарная нормальная реакция практически не зависит от выбранного "
            f"контактного метода. На основной расчетной сетке штрафной метод дает RN={pen_k_n} кН, "
            f"а метод расширенного Лагранжа — RN={aug_k_n} кН. Относительное различие можно оценить как"
        ),
        564: f"δR≈|{aug_k_n}-{pen_k_n}|/{aug_k_n}⋅100%≈{delta_r:.2f}%.(4.26)",
        566: (
            "Принципиальное различие проявляется в выполнении условия непроникания. "
            f"На основной расчетной сетке максимальная пенетрация для штрафного метода равна {pen_g} мм, "
            f"а для метода расширенного Лагранжа — {aug_g} мм. Отношение этих величин составляет"
        ),
        567: f"gpenpen/gpenAL≈{pen_g}/{aug_g}≈{fmt_sci(g_ratio, 2)}.(4.27)",
        571: (
            "Однако повышение точности контактного ограничения достигается за счет вычислительной стоимости. "
            f"На основной сетке штрафной метод потребовал {pen.nonlinear_iterations} нелинейных итераций "
            f"и {fmt_time(pen.total_time_seconds)} с расчетного времени. Метод расширенного Лагранжа потребовал "
            f"{aug.nonlinear_iterations} нелинейных итераций и {fmt_time(aug.total_time_seconds)} с. "
            "Отношение времен равно"
        ),
        572: (
            f"tAL/tpen={fmt_time(aug.total_time_seconds)}/{fmt_time(pen.total_time_seconds)}"
            f"≈{fmt_ratio(time_ratio)}.(4.28)"
        ),
        573: f"Таким образом, метод расширенного Лагранжа оказался примерно в {time_ratio:.1f} раза дороже по времени. Сравнение расчетного времени показано на рис. 4.5.",
        577: (
            "Также заметно, что метод расширенного Лагранжа требует большего числа шагов нагружения. "
            f"На основной сетке штрафной метод завершил расчет за {pen.load_steps} шагов, тогда как метод "
            f"расширенного Лагранжа потребовал {aug.load_steps} шагов. Это связано с большей чувствительностью "
            "метода к изменению активной контактной зоны и необходимостью устойчивого обновления множителей Лагранжа."
        ),
        584: (
            "Полученные результаты показывают, что оба контактных алгоритма приводят к близкому интегральному "
            "силовому ответу. Это означает, что для рассматриваемого перемещения жесткой плоскости оба метода "
            "в целом одинаково передают нагрузку от шины к жесткой плоскости. На основной расчетной сетке "
            f"суммарная реакция отличается примерно на {delta_r:.2f} %, то есть менее чем на 0.2 %, что для "
            "инженерного расчета можно считать практически совпадающим результатом."
        ),
        587: f"gpen,maxpen={pen_g} мм,",
        589: f"gpen,maxAL={aug_g} мм.",
        594: (
            "Оба метода дали практически одинаковую длину контактной зоны. Это связано с тем, что при выбранном "
            "шаге сетки активная область определяется одними и теми же граничными фасетками, а различие между "
            "методами проявляется прежде всего в величине остаточной пенетрации."
        ),
        598: f"δL≈|{pen_l}-{aug_l}|/{aug_l}⋅100%≈{length_delta:.2f}%.(4.30)",
        599: (
            "Следовательно, по длине пятна контакта методы совпадают в пределах округления. "
            "Различие между ними в рассматриваемой постановке связано главным образом с различной точностью "
            "выполнения контактного ограничения."
        ),
        600: (
            "С вычислительной точки зрения штрафной метод оказался заметно более экономичным. Он потребовал меньше "
            "шагов нагружения и меньше времени расчета. Метод расширенного Лагранжа потребовал большего числа шагов "
            "и внешних итераций, поскольку помимо равновесия по перемещениям необходимо добиться стабилизации "
            f"контактных множителей. На основной сетке расчет методом расширенного Лагранжа оказался примерно "
            f"в {time_ratio:.1f} раза дольше."
        ),
        607: (
            "Метод расширенного Лагранжа обеспечил значительно более точное выполнение контактного ограничения. "
            f"На основной сетке максимальная пенетрация уменьшилась с {pen_g} мм для штрафного метода до "
            f"{aug_g} мм. При этом суммарная нормальная реакция практически не изменилась. Это означает, что метод "
            "расширенного Лагранжа улучшает локальное выполнение условия непроникания, не нарушая общий силовой баланс задачи."
        ),
        638: (
            "Для основной расчетной сетки при перемещении жесткой опорной плоскости на 10 мм и начальном "
            "зазоре 2 мм суммарная нормальная реакция составила"
        ),
        639: f"RN={pen_k_n} кН.",
        641: f"gpen,max={pen_g} мм.",
        643: f"Lc={pen_l} мм.",
        645: (
            "Поле перемещений позволяет оценить общую форму деформирования сектора под действием перемещения "
            "опорной плоскости. Наибольшие перемещения возникают в зоне взаимодействия с опорной плоскостью, "
            "а по мере удаления от наружного контура их величина уменьшается. Соответствующее поле перемещений "
            "показано на рис. 5.2. По этому рисунку можно оценить общую форму осадки сектора и локализацию "
            "перемещений вблизи опорной поверхности."
        ),
        668: f"RN={aug_k_n} кН,",
        670: f"gpen,max={aug_g} мм.",
        672: f"Lc={aug_l} мм.",
        707: (
            f"где R — внешний радиус шины, а φc выражается в радианах. Для R=300 мм и "
            f"φc={fmt_angle(main.half_angle_deg)}°={main.half_angle_rad:.6f} рад получаем"
        ),
        708: f"a≈300·{main.half_angle_rad:.6f}≈{main.half_width:.1f} мм,",
        709: f"что согласуется с восстановленным значением {fmt_width(main.half_width)} мм.",
        714: (
            "Для дальнейшего сопоставления используется основная расчетная сетка. На ней максимальное давление "
            f"составило pmax={fmt_pressure(main.max_pressure)} МПа, полуугол пятна контакта — "
            f"{fmt_angle(main.half_angle_deg)}°, а полуширина пятна контакта — {fmt_width(main.half_width)} мм."
        ),
        715: (
            "Восстановленная нагрузка не является независимой физической моделью контакта. Она наследует форму, "
            "положение и интенсивность от уже решенной контактной задачи. Поэтому ее смысл состоит не в определении "
            "контактного состояния, а в проверке того, насколько заданное давление способно заменить найденное "
            "контактное действие при последующем расчете полей перемещений, напряжений и энергии деформации."
        ),
        717: (
            "После восстановления параметров давления решается гиперупругая задача без явного контакта на той же "
            "геометрии и с теми же параметрами материала. В этой постановке жесткая плоскость уже не входит в систему "
            "ограничений: внутренний контур закрепляется, а действие опоры заменяется нормальной поверхностной "
            "нагрузкой на внешней дуге сектора."
        ),
        722: (
            "где Γp — участок внешнего контура, на котором прикладывается восстановленное давление, N — матрица "
            "функций формы на граничной фасетке, p(ξ) — заданное параболическое давление, n — внешняя нормаль к "
            "поверхности. Схема приложения восстановленной нагрузки показана на рис. 5.14. В отличие от явной "
            "контактной задачи, здесь нет проверки зазора и перестроения активного множества: нагруженная дуга "
            "выбирается заранее по результатам контактного расчета."
        ),
        727: (
            "В этой системе отсутствуют неизвестные контактные реакции, множители Лагранжа и проверка активного "
            "множества. Нелинейность сохраняется только за счет гиперупругой модели материала и конечно-деформационной "
            "кинематики. Поэтому расчет становится проще, но одновременно теряет способность самостоятельно находить "
            "границу контакта и перераспределять давление при изменении деформированной конфигурации."
        ),
        729: f"pmax={fmt_pressure(main.max_pressure)} МПа,",
        730: f"φc={fmt_angle(main.half_angle_deg)}°={main.half_angle_rad:.6f} рад,(5.8)",
        731: f"a={fmt_width(main.half_width)} мм.",
        732: (
            "Поверхностное давление прикладывалось только в пределах дуги |ξ|≤a, а вне этой области принималось "
            "равным нулю. Таким образом, расчетная задача без явного контакта воспроизводит заданное интегральное "
            "действие найденного контакта, но не решает задачу контакта как часть нелинейной системы."
        ),
        733: (
            "Главное вычислительное преимущество такой постановки состоит в сокращении времени расчета. "
            f"Для основной расчетной сетки расчетная задача без явного контакта с восстановленной нагрузкой была "
            f"решена за {fmt_time(sur.total_time_seconds)} с. Для сравнения, явный контакт штрафным методом на той "
            f"же сетке потребовал {fmt_time(pen.total_time_seconds)} с, а метод расширенного Лагранжа — "
            f"{fmt_time(aug.total_time_seconds)} с. Отношение времен составляет"
        ),
        734: (
            f"tpen/ts≈{fmt_time(pen.total_time_seconds)}/{fmt_time(sur.total_time_seconds)}"
            f"≈{fmt_ratio(sur_pen_ratio)},(5.9)"
        ),
        735: (
            f"tAL/ts≈{fmt_time(aug.total_time_seconds)}/{fmt_time(sur.total_time_seconds)}"
            f"≈{fmt_ratio(sur_aug_ratio)}.(5.10)"
        ),
        736: (
            "Следовательно, бесконтактная аппроксимация действительно снижает вычислительную стоимость, особенно "
            "по сравнению с методом расширенного Лагранжа. При этом выигрыш по отношению к штрафному методу оказался "
            "умеренным: в рассматриваемом расчете время уменьшилось примерно в 1.8 раза."
        ),
        737: (
            "Это ускорение достигается за счет потери самостоятельности модели. Если изменить геометрию, материал, "
            "положение плоскости, величину перемещения или начальный зазор, ранее восстановленное давление уже не "
            "обязано соответствовать новому контактному состоянию. Поэтому такая аппроксимация применима только тогда, "
            "когда давление заранее известно или может быть надежно получено из отдельного явного контактного расчета."
        ),
        738: (
            "С инженерной точки зрения расчет без явного контакта с восстановленной нагрузкой удобен не как замена "
            "контактной постановки, а как быстрый последующий расчет при фиксированном контактном сценарии. Он может "
            "использоваться для оценки влияния уже найденного давления на поля напряжений и деформаций, но не для "
            "самостоятельного определения контактных характеристик."
        ),
        740: (
            "Сравнение явного контакта и бесконтактной аппроксимации должно выполняться не по контактным величинам, "
            "а по полям напряженно-деформированного состояния внутри тела шины. Это принципиально важно: в расчетной "
            "задаче без явного контакта нет активного множества, нормального зазора и контактной реакции как неизвестных "
            "задачи. Поэтому проверяется не совпадение контактного алгоритма, а то, насколько заданная нагрузка передает "
            "механический эффект найденного контактного давления."
        ),
        745: (
            "В качестве основных сравниваемых величин выбраны поле перемещений, эквивалентные напряжения по Мизесу, "
            "плотность энергии деформации, силовые характеристики, а также профили напряжений по контуру и по "
            "радиальному сечению. Такой набор позволяет отдельно оценить общую деформационную картину и локальные "
            "отличия в зоне приложения восстановленной нагрузки."
        ),
        746: (
            "Поле перемещений является первой характеристикой для проверки. Бесконтактная нагрузка должна в первую "
            "очередь воспроизводить общий характер осадки и распределение перемещений в теле шины, но от нее не "
            "следует ожидать полного совпадения в окрестности краев пятна контакта. Сравнение полей перемещений "
            "показано на рис. 5.15."
        ),
        749: (
            "Следующим шагом сравниваются силовые характеристики на внешнем контуре. В явном контакте они формируются "
            "контактным алгоритмом, а в расчетной задаче без явного контакта задаются через эквивалентное давление. "
            "Сравнение модулей контактных и эквивалентных нагрузочных сил приведено на рис. 5.16. Этот график позволяет "
            "оценить не только суммарный уровень нагрузки, но и отличие формы распределения по дуге."
        ),
        753: (
            "Эквивалентные напряжения по Мизесу используются как обобщенная характеристика интенсивности напряженного "
            "состояния. Хотя для резиноподобного материала эта величина не является критерием разрушения в строгом "
            "смысле, она удобна для сопоставления локализации напряжений между различными расчетными постановками "
            "[1, 5]. Соответствующее сравнение показано на рис. 5.17. По нему видно, насколько восстановленная "
            "нагрузка передает общий уровень напряженного состояния и где возникают локальные расхождения."
        ),
        756: (
            "Дополнительно анализируется плотность энергии деформации. Эта величина удобна для гиперупругой постановки, "
            "поскольку непосредственно связана с выбранным потенциалом материала W [4, 5]. Сравнение на рис. 5.18 "
            "показывает, насколько близко две постановки воспроизводят накопленную упругую энергию в теле шины; "
            "локальные отличия при этом ожидаемы в зоне приложения давления."
        ),
        760: (
            "Для более детального анализа используются графики профилей напряжений. Профили по внешнему контуру "
            "позволяют оценить расхождения непосредственно вблизи зоны нагружения, а профили по радиальному сечению "
            "показывают, как отличие поверхностной нагрузки передается вглубь резинового массива. Эти результаты "
            "приведены на рис. 5.19 и 5.20."
        ),
        767: (
            "Основное отличие бесконтактной аппроксимации состоит в том, что она фиксирует нагрузку заранее и поэтому "
            "не реагирует на изменение активной зоны в ходе решения. Из-за этого она может передавать интегральный "
            "уровень нагрузки и общий характер деформирования, но локально расходиться с явным контактом, особенно "
            "вблизи краев пятна контакта и на внешнем контуре [3]."
        ),
        768: (
            "Таким образом, бесконтактная аппроксимация с восстановленной нагрузкой может использоваться для быстрой "
            "оценки напряженно-деформированного состояния при уже известном контактном давлении. Ее нельзя рассматривать "
            "как самостоятельный контактный метод: она не определяет активную область, не контролирует зазор и не "
            "обновляет давление при изменении условий нагружения."
        ),
        777: (
            "Наконец, бесконтактная аппроксимация с восстановленной нагрузкой не является самостоятельной физической "
            "моделью контакта. Она не определяет пятно контакта и давление из условий непроникания, а использует данные, "
            "уже полученные из явного контактного расчета. Поэтому ее следует применять только как инженерное упрощение "
            "для известного контактного состояния и не переносить без повторной калибровки на другие режимы нагружения."
        ),
        778: (
            "Таким образом, разработанная модель предназначена для сопоставления контактных алгоритмов и оценки "
            "ограниченной возможности замены явного контакта восстановленной нагрузкой. Она дает механически "
            "содержательные выводы в рамках принятой двумерной постановки, но не предназначена для полного описания "
            "реальной шины в условиях качения, трения и трехмерного контактного взаимодействия."
        ),
        780: (
            "В главах 4 и 5 было выполнено сопоставление двух явных контактных алгоритмов и бесконтактной аппроксимации "
            "с восстановленной нагрузкой. Полученные результаты позволяют разделить расчетные подходы по их роли: "
            "контактные методы определяют контактное состояние, а бесконтактная постановка использует уже найденное "
            "давление как заданное внешнее воздействие."
        ),
        781: (
            "Штрафной метод показал себя как наиболее простая и устойчивая схема для получения контактного решения. "
            "Он обеспечивает близкий к методу расширенного Лагранжа интегральный силовой ответ, но допускает конечную "
            "пенетрацию. Поэтому его рационально использовать для предварительных и серийных расчетов, где важны "
            "устойчивость и умеренная вычислительная стоимость."
        ),
        782: (
            "Метод расширенного Лагранжа обеспечивает существенно более точное выполнение условия непроникания. "
            "На основной расчетной сетке при практически совпадающей суммарной реакции он уменьшает максимальную "
            "пенетрацию примерно на шесть порядков по сравнению со штрафным методом. Это делает его более подходящим "
            "для проверочных расчетов и восстановления контактного давления."
        ),
        783: (
            "Бесконтактная аппроксимация с восстановленной нагрузкой не является самостоятельной контактной моделью, "
            "так как не определяет активную область контакта и распределение давления. После калибровки по явному "
            "контактному расчету она может сократить время последующего анализа, но ее результаты следует трактовать "
            "как приближение к известному контактному состоянию, а не как новое решение контактной задачи."
        ),
        784: (
            "Таким образом, основной результат расчетной части состоит в построении согласованной схемы применения "
            "трех подходов: штрафной метод — для устойчивого контактного расчета, метод расширенного Лагранжа — для "
            "более точного определения контактного давления, бесконтактная аппроксимация с восстановленной нагрузкой — "
            "для экономичного последующего анализа при уже найденном и зафиксированном контактном давлении."
        ),
        786: (
            "В работе решена задача построения и исследования конечно-элементной модели контактного взаимодействия "
            "массивной шины с жесткой плоской опорной поверхностью. Расчетная область была принята в виде двумерного "
            "кольцевого сектора, моделирующего поперечное сечение шины, а расчет выполнялся в постановке плоской "
            "деформации. Такая схема позволила сосредоточиться на нормальном контактном взаимодействии и сопоставлении "
            "численных методов его учета."
        ),
        788: (
            "Вторая задача была связана с формированием расчетной схемы контактной задачи. В работе была задана "
            "геометрия кольцевого сектора, выделена потенциальная контактная область на внешнем контуре, построена "
            "сетка с локальным сгущением в зоне контакта и введена жесткая плоская опорная поверхность. В расчетной "
            "части использовалась схема с закрепленным внутренним контуром, начальным зазором 2 мм и перемещением "
            "жесткой плоскости на 10 мм. Контакт рассматривался как гладкий односторонний нормальный контакт без трения."
        ),
        789: (
            "Третья задача состояла в переходе от линейно-упругой модели к конечно-деформационной гиперупругой "
            "постановке. Для описания резиноподобного материала использовалась сжимаемая неогуковская модель "
            "в почти несжимаемом режиме. Параметры материала задавались через модуль Юнга и коэффициент Пуассона "
            "с последующим переходом к модулю сдвига, модулю объемного сжатия и параметру λ. Значение ν=0.48 "
            "привело к высокой объемной жесткости, что было учтено при анализе обусловленности и устойчивости "
            "численного решения."
        ),
        791: (
            "Пятая задача состояла в сопоставлении штрафного метода и метода расширенного Лагранжа. На основной "
            f"расчетной сетке штрафной метод дал суммарную нормальную реакцию {pen_k_n} кН при максимальной "
            f"пенетрации {pen_g} мм. Метод расширенного Лагранжа дал близкую суммарную реакцию {aug_k_n} кН, "
            f"но уменьшил максимальную пенетрацию до {aug_g} мм. Следовательно, интегральный силовой ответ двух "
            "методов практически совпал, тогда как точность выполнения локального условия непроникания различалась "
            "примерно на шесть порядков."
        ),
        792: (
            "При этом метод расширенного Лагранжа оказался существенно дороже по времени. "
            f"На основной расчетной сетке штрафной расчет был выполнен за {fmt_time(pen.total_time_seconds)} с, "
            f"а расчет методом расширенного Лагранжа — за {fmt_time(aug.total_time_seconds)} с. Увеличение времени "
            "связано с внешними итерациями по множителям, стабилизацией обновления контактного давления и более "
            "осторожным пошаговым нагружением."
        ),
        793: (
            "Шестая задача заключалась в построении бесконтактной аппроксимации контактного нагружения. Для этого "
            "по результату метода расширенного Лагранжа была восстановлена эквивалентная параболическая нагрузка "
            f"на внешнем контуре шины. Для основной расчетной сетки были получены параметры: максимальное давление "
            f"pmax={fmt_pressure(main.max_pressure)} МПа, полуугол пятна контакта {fmt_angle(main.half_angle_deg)}° "
            f"и полуширина пятна контакта {fmt_width(main.half_width)} мм."
        ),
        794: (
            "Седьмая задача состояла в оценке применимости этой аппроксимации. Гиперупругая расчетная задача без "
            f"явного контакта с восстановленной нагрузкой была решена за {fmt_time(sur.total_time_seconds)} с, что "
            f"примерно в {fmt_ratio(sur_pen_ratio)} раза быстрее штрафного контактного расчета и примерно в "
            f"{fmt_ratio(sur_aug_ratio)} раза быстрее расчета методом расширенного Лагранжа. Поэтому бесконтактная "
            "аппроксимация может быть полезна для быстрого повторного анализа напряженно-деформированного состояния, "
            "если контактное давление уже известно из предварительного явного контактного расчета."
        ),
        795: (
            "В то же время такая аппроксимация не является самостоятельной заменой контактной задачи. Она не определяет "
            "активную область контакта и распределение давления из условий непроникания, а использует эти данные как "
            "исходные. Поэтому ее рациональная область применения ограничена задачами, где контактное состояние уже "
            "предварительно найдено или надежно задано."
        ),
        797: (
            "Таким образом, цель работы достигнута. Сформирована конечно-элементная постановка контактной задачи "
            "массивной шины с жесткой плоской опорой, выполнено сопоставление штрафного метода и метода расширенного "
            "Лагранжа, а также показана возможность построения быстрой, но зависимой от предварительного контактного "
            "расчета бесконтактной аппроксимации. Полученные результаты позволяют разделить области применения "
            "расчетных подходов: штрафной метод целесообразен для устойчивых инженерных расчетов, метод расширенного "
            "Лагранжа — для более точного контактного решения и восстановления давления, а бесконтактная аппроксимация — "
            "для экономичного последующего анализа уже известного контактного состояния."
        ),
        798: "7. СПИСОК ИСПОЛЬЗОВАННЫХ ИСТОЧНИКОВ",
    }

    for index, value in replacements.items():
        if index >= len(body_children):
            raise RuntimeError(f"Body index is out of range: {index}")
        set_element_text(body_children[index], value)


def patch_images(
    body_children: list[ET.Element],
    rels_root: ET.Element,
) -> dict[str, Path]:
    caption_to_figure = {
        "4.2": "4.2",
        "4.3": "4.3",
        "4.4": "4.4",
        "4.5": "4.5",
        "4.6": "4.6",
        "4.7": "4.7",
        "5.1": "5.1",
        "5.2": "5.2",
        "5.4": "5.4",
        "5.5": "5.5",
        "5.6": "5.6",
        "5.7": "5.7",
        "5.8": "5.8",
        "5.10": "5.10",
        "5.11": "5.11",
        "5.12": "5.12",
        "5.13": "5.13",
        "5.14": "5.15",
        "5.15": "5.16",
        "5.16": "5.17",
        "5.17": "5.18",
        "5.18": "5.19",
        "5.19": "5.20",
        "5.20": "5.21",
    }
    caption_prefix = "Рисунок "
    replacements: dict[str, Path] = {}

    for index, child in enumerate(body_children):
        if child.tag != q(W, "p"):
            continue
        caption = element_text(child)
        if not caption.startswith(caption_prefix):
            continue
        parts = caption.split()
        if len(parts) < 2:
            continue
        number = parts[1]
        figure_number = caption_to_figure.get(number)
        if figure_number is None:
            continue
        chapter, local_number = figure_number.split(".")
        image_path = FIGURES / f"fig_{int(chapter)}_{int(local_number):02d}.png"
        if not image_path.exists():
            raise FileNotFoundError(image_path)

        image_paragraph = None
        for previous in reversed(body_children[:index]):
            rel_id = blip_rel_id(previous)
            if rel_id:
                image_paragraph = previous
                break
        if image_paragraph is None:
            raise RuntimeError(f"Image paragraph not found before caption {number}")
        rel_id = blip_rel_id(image_paragraph)
        if rel_id is None:
            raise RuntimeError(f"Image relationship not found before caption {number}")
        target = find_relationship_target(rels_root, rel_id)
        replacements[f"word/{target}"] = image_path
        resize_drawing(image_paragraph, image_path)

    if len(replacements) != len(caption_to_figure):
        raise RuntimeError(
            f"Expected {len(caption_to_figure)} image replacements, got {len(replacements)}"
        )
    return replacements


def count_captions(body_children: list[ET.Element], prefix: str) -> int:
    return sum(
        1
        for child in body_children
        if child.tag == q(W, "p") and element_text(child).startswith(prefix)
    )


def count_sources(body_children: list[ET.Element]) -> int:
    in_sources = False
    total = 0
    for child in body_children:
        text = element_text(child)
        if text == "СПИСОК ЛИТЕРАТУРЫ" or text == "7. СПИСОК ИСПОЛЬЗОВАННЫХ ИСТОЧНИКОВ":
            in_sources = True
            continue
        if in_sources and re.match(r"^\d+\.\s+", text):
            total += 1
    return total


def update_abstract_counts(body_children: list[ET.Element]) -> tuple[int, int, int, int]:
    if len(body_children) <= 2:
        raise RuntimeError("Document body is unexpectedly short")
    current = element_text(body_children[2])
    page_match = re.search(r"содержит\s+(\d+)\s+страниц", current)
    pages = int(page_match.group(1)) if page_match else 99
    figures = count_captions(body_children, "Рисунок ")
    tables = count_captions(body_children, "Таблица ")
    sources = count_sources(body_children)
    set_element_text(
        body_children[2],
        (
            "Выпускная квалификационная работа содержит "
            f"{pages} страниц, {figures} рисунка, {tables} таблиц, "
            f"{sources} использованных источников."
        ),
    )
    return pages, figures, tables, sources


def update_docx() -> tuple[int, int]:
    main = read_triplet(MAIN_RESULTS)
    coarse = read_triplet(COARSE_RESULTS)

    with zipfile.ZipFile(SOURCE_DOCX, "r") as zin:
        document_bytes = zin.read("word/document.xml")
        register_namespaces_from_xml(document_bytes)
        document_root = ET.fromstring(document_bytes)
        rels_root = ET.fromstring(zin.read("word/_rels/document.xml.rels"))
        body = document_root.find(q(W, "body"))
        if body is None:
            raise RuntimeError("word/document.xml has no w:body")
        body_children = list(body)

        patch_tables(body_children, main, coarse)
        patch_paragraphs(body_children, main)
        media_replacements = patch_images(body_children, rels_root)
        collapsed_ref_fields = collapse_all_ref_fields(body_children)
        renumbered_references = renumber_figures_and_references(body_children)
        abstract_counts = update_abstract_counts(body_children)
        if document_root.get(q(MC, "Ignorable")):
            document_root.set(q(MC, "Ignorable"), "w14 wp14")

        document_xml = ET.tostring(document_root, encoding="utf-8", xml_declaration=True)
        rels_xml = ET.tostring(rels_root, encoding="utf-8", xml_declaration=True)

        tmp = OUTPUT_DOCX.with_suffix(".tmp.docx")
        if tmp.exists():
            tmp.unlink()
        with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename == "word/document.xml":
                    continue
                if item.filename == "word/_rels/document.xml.rels":
                    continue
                if item.filename in media_replacements:
                    zout.writestr(item, media_replacements[item.filename].read_bytes())
                else:
                    zout.writestr(item, zin.read(item.filename))
            zout.writestr("word/document.xml", document_xml)
            zout.writestr("word/_rels/document.xml.rels", rels_xml)

    if OUTPUT_DOCX.exists():
        OUTPUT_DOCX.unlink()
    shutil.move(str(tmp), str(OUTPUT_DOCX))
    return len(media_replacements), collapsed_ref_fields, renumbered_references, abstract_counts


def main() -> None:
    image_count, collapsed_ref_fields, renumbered_references, abstract_counts = update_docx()
    with zipfile.ZipFile(OUTPUT_DOCX, "r") as archive:
        bad_file = archive.testzip()
    if bad_file is not None:
        raise RuntimeError(f"Bad file in generated docx: {bad_file}")
    print(f"source={SOURCE_DOCX}")
    print(f"output={OUTPUT_DOCX}")
    print(f"updated_images={image_count}")
    print(f"collapsed_ref_fields={collapsed_ref_fields}")
    print(f"renumbered_figure_paragraphs={renumbered_references}")
    print(
        "abstract_counts="
        f"pages:{abstract_counts[0]},figures:{abstract_counts[1]},"
        f"tables:{abstract_counts[2]},sources:{abstract_counts[3]}"
    )
    print("zip_check=ok")


if __name__ == "__main__":
    main()
