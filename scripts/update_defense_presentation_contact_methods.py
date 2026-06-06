#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import os
import re
import shutil
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_DIR = REPO_ROOT / "docs"
ASSET_DIR = DOCS_DIR / "presentation_assets"
PRESENTATION_GLOB = "MironovMV_RK5-81B_*.pptx"
NORMAL_GAP_SOURCE_IMAGE = ASSET_DIR / "normal_gap_literature.png"
NORMAL_GAP_PRESENTATION_MEDIA = "ppt/media/image-normal-gap-literature.png"

P_NS = "http://schemas.openxmlformats.org/presentationml/2006/main"
A_NS = "http://schemas.openxmlformats.org/drawingml/2006/main"
R_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CT_NS = "http://schemas.openxmlformats.org/package/2006/content-types"
SLIDE_FIGSIZE = (14.2, 4.85)

NS = {
    "p": P_NS,
    "a": A_NS,
    "r": R_NS,
    "rel": PKG_REL_NS,
    "ct": CT_NS,
}

SLIDE_CONTENT_TYPE = "application/vnd.openxmlformats-officedocument.presentationml.slide+xml"
SLIDE_REL_TYPE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide"
IMAGE_REL_TYPE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image"
LAYOUT_REL_TYPE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout"


def qname(namespace: str, name: str) -> str:
    return f"{{{namespace}}}{name}"


def register_namespaces() -> None:
    ET.register_namespace("p", P_NS)
    ET.register_namespace("a", A_NS)
    ET.register_namespace("r", R_NS)
    ET.register_namespace("", PKG_REL_NS)
    ET.register_namespace("ct", CT_NS)


def presentation_path() -> Path:
    candidates = sorted(DOCS_DIR.glob(PRESENTATION_GLOB), key=lambda path: path.stat().st_mtime, reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No presentation matching {PRESENTATION_GLOB} in {DOCS_DIR}")
    return candidates[0]


def configure_matplotlib() -> None:
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#222222",
            "axes.linewidth": 1.0,
            "axes.titlesize": 18,
            "axes.labelsize": 14,
            "axes.titleweight": "semibold",
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman"],
            "font.size": 12,
            "legend.fontsize": 10.5,
            "legend.frameon": True,
            "legend.facecolor": "white",
            "legend.edgecolor": "#bbbbbb",
            "legend.framealpha": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "savefig.dpi": 220,
            "axes.formatter.use_mathtext": True,
            "mathtext.fontset": "stix",
        }
    )


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def mesh_label(raw: str) -> str:
    labels = {
        "contact_focused_very_coarse": "очень грубая",
        "contact_focused_coarse": "грубая",
        "contact_focused_medium": "средняя",
        "contact_focused_dense": "плотная",
        "contact_focused_heavy": "тяжелая",
    }
    return labels.get(raw, raw)


def method_label(raw: str) -> str:
    if raw == "penalty":
        return "Штрафной метод"
    if raw == "augmented_lagrangian":
        return "Модифицированный метод Лагранжа"
    return raw


def method_style(raw: str) -> dict[str, str]:
    if raw == "penalty":
        return {"linestyle": "-", "marker": "o"}
    return {"linestyle": "--", "marker": "s"}


def build_parameter_study_figure(output_path: Path) -> None:
    rows = read_csv_rows(REPO_ROOT / "results" / "main_scale_contact_method_comparison_mesh_sensitivity" / "study_summary.csv")
    success_rows = [row for row in rows if row.get("success", "").lower() == "true"]
    meshes = [
        "contact_focused_very_coarse",
        "contact_focused_coarse",
        "contact_focused_medium",
        "contact_focused_dense",
        "contact_focused_heavy",
    ]
    colors = {
        "contact_focused_very_coarse": "#a61c1c",
        "contact_focused_coarse": "#1f4e79",
        "contact_focused_medium": "#b45f06",
        "contact_focused_dense": "#2f6b3b",
        "contact_focused_heavy": "#0b6e75",
    }

    fig, axes = plt.subplots(1, 2, figsize=SLIDE_FIGSIZE)
    fields = [
        ("max_penetration", "Максимальное проникание, мм"),
        ("total_time_seconds", "Полное время расчета, с"),
    ]

    for mesh in meshes:
        for method in ["penalty", "augmented_lagrangian"]:
            group = sorted(
                [
                    row
                    for row in success_rows
                    if row["mesh"] == mesh and row["contact_method"] == method
                ],
                key=lambda row: float(row["contact_parameter_value"]),
            )
            if not group:
                continue
            x_values = np.array([float(row["contact_parameter_value"]) for row in group])
            style = method_style(method)
            for axis, (field_name, ylabel) in zip(axes, fields):
                y_values = np.array([max(float(row[field_name]), 1.0e-10) for row in group])
                axis.plot(
                    x_values,
                    y_values,
                    color=colors[mesh],
                    linestyle=style["linestyle"],
                    marker=style["marker"],
                    linewidth=1.85,
                    markersize=4.8,
                    markerfacecolor="white",
                    markeredgewidth=1.0,
                    alpha=0.95,
                )
                axis.set_xscale("log")
                axis.set_yscale("log")
                axis.set_xlabel("Контактный параметр")
                axis.set_ylabel(ylabel)
                axis.grid(True, which="major", linestyle="--", linewidth=0.65, color="#aeb6bf", alpha=0.9)
                axis.grid(True, which="minor", linestyle=":", linewidth=0.45, color="#d5dbe3", alpha=0.85)

    axes[0].set_title("Контроль проникания")
    axes[1].set_title("Вычислительная цена")

    mesh_handles = [
        plt.Line2D([0], [0], color=colors[mesh], linewidth=2.0, label=mesh_label(mesh))
        for mesh in meshes
    ]
    method_handles = [
        plt.Line2D(
            [0],
            [0],
            color="#444444",
            linewidth=2.0,
            linestyle=method_style(method)["linestyle"],
            marker=method_style(method)["marker"],
            markerfacecolor="white",
            label=method_label(method),
        )
        for method in ["penalty", "augmented_lagrangian"]
    ]
    legend1 = fig.legend(
        handles=mesh_handles,
        loc="upper center",
        bbox_to_anchor=(0.31, 0.90),
        ncol=3,
        title="Сетка",
    )
    fig.add_artist(legend1)
    fig.legend(
        handles=method_handles,
        loc="upper center",
        bbox_to_anchor=(0.77, 0.90),
        ncol=2,
        title="Метод контакта",
    )
    fig.suptitle("Серия расчетов: 5 сеток × 10 контактных параметров × 2 метода", y=0.99, fontsize=18, fontweight="semibold")
    fig.subplots_adjust(left=0.075, right=0.975, bottom=0.16, top=0.62, wspace=0.25)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def parse_triplet_summary(path: Path) -> dict[str, dict[str, float | str]]:
    rows: dict[str, dict[str, float | str]] = {}
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            if not row.get("case"):
                break
            rows[row["case"]] = row
    return rows


def build_final_contact_figure(output_path: Path) -> None:
    root = REPO_ROOT / "results" / "main_scale_hyperelastic_reference_triplet_coarse_moving_plane"
    summary = parse_triplet_summary(root / "triplet_summary.csv")
    pressure_rows = {
        "penalty_contact": read_csv_rows(root / "penalty_contact" / "outer_boundary_contact_pressure.csv"),
        "augmented_lagrangian_contact": read_csv_rows(root / "augmented_lagrangian_contact" / "outer_boundary_contact_pressure.csv"),
    }

    fig = plt.figure(figsize=SLIDE_FIGSIZE)
    grid = fig.add_gridspec(1, 3, width_ratios=[2.2, 1.0, 1.0], wspace=0.35)
    pressure_axis = fig.add_subplot(grid[0, 0])
    penetration_axis = fig.add_subplot(grid[0, 1])
    time_axis = fig.add_subplot(grid[0, 2])

    colors = {
        "penalty_contact": "#1f4e79",
        "augmented_lagrangian_contact": "#a61c1c",
    }
    labels = {
        "penalty_contact": "Штрафной метод",
        "augmented_lagrangian_contact": "Модифицированный метод Лагранжа",
    }
    for case_name, rows in pressure_rows.items():
        x_values = np.array([float(row["arc_coordinate"]) for row in rows])
        y_values = np.array([float(row["average_pressure_mpa"]) for row in rows])
        pressure_axis.plot(
            x_values,
            y_values,
            color=colors[case_name],
            linewidth=2.3,
            marker="o" if case_name == "penalty_contact" else "s",
            markersize=4.0,
            markevery=max(1, len(rows) // 16),
            markerfacecolor="white",
            label=labels[case_name],
        )

    pressure_axis.set_title("Давление по пятну контакта")
    pressure_axis.set_xlabel("Дуговая координата, мм")
    pressure_axis.set_ylabel("Контактное давление, МПа")
    pressure_axis.grid(True, which="major", linestyle="--", linewidth=0.65, color="#aeb6bf", alpha=0.9)
    pressure_axis.grid(True, which="minor", linestyle=":", linewidth=0.45, color="#d5dbe3", alpha=0.85)
    pressure_axis.legend(loc="best")

    case_order = ["penalty_contact", "augmented_lagrangian_contact"]
    short_labels = ["Штрафной\nметод", "Модифицированный\nметод Лагранжа"]
    max_penetration = [max(float(summary[case]["max_penetration"]), 1.0e-10) for case in case_order]
    total_time = [float(summary[case]["total_time_seconds"]) for case in case_order]

    penetration_axis.bar(
        short_labels,
        max_penetration,
        color=[colors[case] for case in case_order],
        width=0.62,
    )
    penetration_axis.set_yscale("log")
    penetration_axis.set_title("Проникание")
    penetration_axis.set_ylabel("Максимум, мм")
    penetration_axis.grid(True, which="major", axis="y", linestyle="--", linewidth=0.65, color="#aeb6bf", alpha=0.9)

    time_axis.bar(
        short_labels,
        total_time,
        color=[colors[case] for case in case_order],
        width=0.62,
    )
    time_axis.set_title("Время")
    time_axis.set_ylabel("Полное время, с")
    time_axis.grid(True, which="major", axis="y", linestyle="--", linewidth=0.65, color="#aeb6bf", alpha=0.9)

    for axis in [penetration_axis, time_axis]:
        axis.tick_params(axis="x", rotation=0, labelsize=10)

    fig.suptitle("Финальная контактная постановка с подвижной плоскостью", y=0.98, fontsize=18, fontweight="semibold")
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.19, top=0.74, wspace=0.35)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def read_cropped_image(path: Path, *, tolerance: float = 0.985, padding: int = 10) -> np.ndarray:
    image = plt.imread(path)
    rgb = image[:, :, :3] if image.ndim == 3 else image
    mask = np.any(rgb < tolerance, axis=2) if rgb.ndim == 3 else rgb < tolerance
    if not np.any(mask):
        return image

    rows, cols = np.nonzero(mask)
    row_min = max(int(rows.min()) - padding, 0)
    row_max = min(int(rows.max()) + padding + 1, image.shape[0])
    col_min = max(int(cols.min()) - padding, 0)
    col_max = min(int(cols.max()) + padding + 1, image.shape[1])
    return image[row_min:row_max, col_min:col_max]


def build_image_pair_figure(
    output_path: Path,
    left_image: Path,
    right_image: Path,
    *,
    left_title: str,
    right_title: str,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=SLIDE_FIGSIZE)
    for axis, image_path, title in [
        (axes[0], left_image, left_title),
        (axes[1], right_image, right_title),
    ]:
        axis.imshow(read_cropped_image(image_path))
        axis.axis("off")
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.04, top=0.96, wspace=0.03)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def build_single_image_figure(output_path: Path, image_path: Path) -> None:
    fig, axis = plt.subplots(1, 1, figsize=SLIDE_FIGSIZE)
    axis.imshow(read_cropped_image(image_path))
    axis.axis("off")
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.04, top=0.96)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def build_contour_stress_profile_figure(output_path: Path, csv_path: Path) -> None:
    rows = read_csv_rows(csv_path)
    contour_labels = {
        "inner": "Внутренний контур",
        "mid": "Средний слой",
        "outer": "Внешний контур",
    }
    contour_styles = {
        "inner": {"color": "#1f4e79", "marker": "o"},
        "mid": {"color": "#b45f06", "marker": "s"},
        "outer": {"color": "#2f6b3b", "marker": "^"},
    }
    fields = [
        ("sigma_rr_mpa", "σ_rr, МПа"),
        ("sigma_tt_mpa", "σ_θθ, МПа"),
        ("tau_rt_mpa", "τ_rθ, МПа"),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(10.8, 5.8), sharex=True)
    for axis, (field_name, ylabel) in zip(axes, fields):
        for contour_name in ["inner", "mid", "outer"]:
            contour_rows = sorted(
                [row for row in rows if row["contour"] == contour_name],
                key=lambda row: float(row["relative_angle_deg"]),
            )
            if not contour_rows:
                continue
            x_values = np.array([float(row["relative_angle_deg"]) for row in contour_rows])
            y_values = np.array([float(row[field_name]) for row in contour_rows])
            style = contour_styles[contour_name]
            axis.plot(
                x_values,
                y_values,
                color=style["color"],
                linewidth=1.9,
                marker=style["marker"],
                markersize=3.2,
                markevery=max(1, len(x_values) // 18),
                markerfacecolor="white",
                markeredgewidth=0.9,
                label=contour_labels[contour_name],
            )
        axis.set_ylabel(ylabel)
        axis.grid(True, which="major", linestyle="--", linewidth=0.6, color="#aeb6bf", alpha=0.9)
        axis.grid(True, which="minor", linestyle=":", linewidth=0.4, color="#d5dbe3", alpha=0.85)
        axis.minorticks_on()

    axes[-1].set_xlabel("Угол относительно центра контакта, град")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.98),
        ncol=3,
        frameon=True,
        framealpha=1.0,
        edgecolor="#777777",
    )
    fig.subplots_adjust(left=0.105, right=0.975, bottom=0.12, top=0.83, hspace=0.32)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def add_formula_card(
    axis: plt.Axes,
    *,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    lines: list[str],
    facecolor: str,
    edgecolor: str,
) -> None:
    card = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=1.0,
        edgecolor="#4a5568",
        facecolor="#fbfbfb",
        transform=axis.transAxes,
    )
    axis.add_patch(card)
    axis.text(
        x + 0.025,
        y + height - 0.055,
        title,
        transform=axis.transAxes,
        fontsize=15,
        fontweight="semibold",
        color="#111111",
        va="top",
    )

    cursor_y = y + height - 0.135
    for line in lines:
        is_formula = line.startswith("$")
        axis.text(
            x + 0.025,
            cursor_y,
            line,
            transform=axis.transAxes,
            fontsize=13.4 if is_formula else 10.6,
            color="#111111" if is_formula else "#222222",
            va="top",
        )
        cursor_y -= 0.082 if is_formula else 0.057


def build_formula_slide_figure(
    output_path: Path,
    cards: list[dict[str, object]],
    *,
    note: str,
) -> None:
    fig = plt.figure(figsize=SLIDE_FIGSIZE)
    axis = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    axis.axis("off")

    for card in cards:
        add_formula_card(axis, **card)

    axis.text(
        0.035,
        0.055,
        note,
        transform=axis.transAxes,
        fontsize=12.4,
        color="#222222",
        va="bottom",
    )
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def build_formula_assets(asset_dir: Path) -> dict[str, Path]:
    outputs = {
        "neo_hookean": asset_dir / "formula_neo_hookean.png",
        "penalty": asset_dir / "formula_penalty_method.png",
        "augmented_lagrangian": asset_dir / "formula_augmented_lagrangian_method.png",
        "inner_outer_loop": asset_dir / "formula_inner_outer_contact_loop.png",
    }

    common_card = {
        "width": 0.29,
        "height": 0.78,
    }
    x_positions = [0.045, 0.355, 0.665]
    y = 0.145

    build_formula_slide_figure(
        outputs["neo_hookean"],
        [
            {
                **common_card,
                "x": x_positions[0],
                "y": y,
                "title": "Кинематика",
                "facecolor": "#edf7f6",
                "edgecolor": "#7fb7b5",
                "lines": [
                    r"$\mathbf{F}=\mathbf{I}+\partial\mathbf{u}/\partial\mathbf{X}$",
                    r"$J=\det\mathbf{F}>0$",
                    r"$\mathbf{C}=\mathbf{F}^{T}\mathbf{F},\quad I_1=\mathrm{tr}\,\mathbf{C}$",
                    "Плоская деформация:",
                    r"$F_{33}=1,\quad \varepsilon_{zz}=0$",
                    "Интегрирование ведется в исходной конфигурации.",
                ],
            },
            {
                **common_card,
                "x": x_positions[1],
                "y": y,
                "title": "Энергия материала",
                "facecolor": "#f3f0ff",
                "edgecolor": "#b197fc",
                "lines": [
                    r"$W=\frac{\mu}{2}(I_1-3-2\ln J)+\frac{\lambda}{2}(\ln J)^2$",
                    r"$\lambda=K-\frac{2}{3}\mu$",
                    r"$\mu=\frac{E}{2(1+\nu)},\quad K=\frac{E}{3(1-2\nu)}$",
                    "Финальные параметры:",
                    "E = 11,84 МПа,   ν = 0,48",
                    "μ = 4,0 МПа,   K = 98,67 МПа",
                ],
            },
            {
                **common_card,
                "x": x_positions[2],
                "y": y,
                "title": "Вычислительная форма",
                "facecolor": "#fff4e6",
                "edgecolor": "#f0ad4e",
                "lines": [
                    "Используется второе напряжение Пиола:",
                    r"$\mathbf{S}=\mu(\mathbf{I}-\mathbf{C}^{-1})+\lambda\ln J\,\mathbf{C}^{-1}$",
                    r"$\mathbf{P}=\mathbf{F}\mathbf{S}$",
                    r"$\mathbf{f}_{int}^{e}=\sum_q \mathbf{P}_q\nabla_0 N\,J_0w_qh$",
                    "Q4-элемент полностью интегрирован.",
                    "Почти несжимаемость может давать объемное запирание.",
                ],
            },
        ],
        note="Неогуковская модель выбрана как первый физически согласованный шаг для резиноподобного материала при известных E и ν.",
    )

    build_formula_slide_figure(
        outputs["penalty"],
        [
            {
                **common_card,
                "x": x_positions[0],
                "y": y,
                "title": "Вариационная постановка",
                "facecolor": "#edf7f6",
                "edgecolor": "#7fb7b5",
                "lines": [
                    r"$\Pi(\mathbf{u})=\int_{\Omega_0}W(\mathbf{F})\,d\Omega-W_{ext}$",
                    r"$\Pi_{\varepsilon}=\Pi+\frac{1}{2}\int_{\Gamma_c}\varepsilon_n\langle-g_n\rangle_+^2\,d\Gamma$",
                    r"$\delta\Pi_{\varepsilon}=0$",
                    "Условие непроникновения заменяется упругой реакцией.",
                ],
            },
            {
                **common_card,
                "x": x_positions[1],
                "y": y,
                "title": "Давление контакта",
                "facecolor": "#f3f0ff",
                "edgecolor": "#b197fc",
                "lines": [
                    r"$g_n=\mathbf{n}\cdot(\mathbf{x}+\mathbf{u}-\mathbf{x}_{pl})$",
                    r"$\delta_n=\langle-g_n\rangle_+$",
                    r"$p_n=\varepsilon_n\delta_n$",
                    r"$\mathbf{t}_c=p_n\mathbf{n}$",
                    "Если зазор положительный,",
                    "контактная реакция равна нулю.",
                ],
            },
            {
                **common_card,
                "x": x_positions[2],
                "y": y,
                "title": "КЭ-вклад",
                "facecolor": "#fff4e6",
                "edgecolor": "#f0ad4e",
                "lines": [
                    r"$\mathbf{f}_c^e=\sum_q \mathbf{N}^{T}\mathbf{t}_c\,J_{\Gamma}w_qh$",
                    r"$\mathbf{K}_c^e=\sum_q\varepsilon_n\,\mathbf{N}^{T}(\mathbf{n}\otimes\mathbf{n})\mathbf{N}\,J_{\Gamma}w_qh$",
                    r"$\mathbf{R}=\mathbf{f}_{int}+\mathbf{f}_c-\mathbf{f}_{ext}=0$",
                    "Рост εn уменьшает проникание,",
                    "но ухудшает обусловленность.",
                ],
            },
        ],
        note="Штрафной метод удобен как инженерная рабочая постановка: он устойчив, но допускает конечное проникание в зоне контакта.",
    )

    build_formula_slide_figure(
        outputs["augmented_lagrangian"],
        [
            {
                **common_card,
                "x": x_positions[0],
                "y": y,
                "title": "Вариационная идея",
                "facecolor": "#edf7f6",
                "edgecolor": "#7fb7b5",
                "lines": [
                    r"$\Pi_{AL}=\Pi+\int_{\Gamma_c}\lambda_n\langle-g_n\rangle_+\,d\Gamma$",
                    r"$\quad+\frac{1}{2}\int_{\Gamma_c}\rho\langle-g_n\rangle_+^2\,d\Gamma$",
                    r"$\delta\Pi_{AL}=0$",
                    "Множитель λn приближает",
                    "неизвестное контактное давление.",
                ],
            },
            {
                **common_card,
                "x": x_positions[1],
                "y": y,
                "title": "Давление в расчете",
                "facecolor": "#f3f0ff",
                "edgecolor": "#b197fc",
                "lines": [
                    r"$p_n^{k+1}=\langle\lambda_n^k-\rho\,g_n(\mathbf{u}^{k+1})\rangle_+$",
                    r"$\mathbf{t}_c=p_n^{k+1}\mathbf{n}$",
                    r"$\lambda_{trial}^{k+1}=p_n^{k+1}$",
                    r"$\lambda_n^{k+1}=\lambda_n^k+\alpha(\lambda_{trial}^{k+1}-\lambda_n^k)$",
                    "В финальных расчетах",
                    "использовалась релаксация множителей.",
                ],
            },
            {
                **common_card,
                "x": x_positions[2],
                "y": y,
                "title": "Итерационный контроль",
                "facecolor": "#fff4e6",
                "edgecolor": "#f0ad4e",
                "lines": [
                    r"$\mathbf{K}_c^e=\sum_q\rho\,\mathbf{N}^{T}(\mathbf{n}\otimes\mathbf{n})\mathbf{N}\,J_{\Gamma}w_qh$",
                    r"$\max(\delta_n)\leq tol_{\delta}$",
                    r"$\|\Delta\lambda\|/\|\lambda\|\leq tol_{\lambda}$",
                    "Активная область уточняется",
                    "после нелинейного шага.",
                    "Цена — больше итераций;",
                    "выигрыш — почти нулевое проникание.",
                ],
            },
        ],
        note="Модифицированный метод Лагранжа использован как уточняющая контактная постановка для контроля непроникновения и восстановления давления.",
    )

    build_formula_slide_figure(
        outputs["inner_outer_loop"],
        [
            {
                **common_card,
                "x": x_positions[0],
                "y": y,
                "title": "Разделение задачи",
                "facecolor": "#edf7f6",
                "edgecolor": "#7fb7b5",
                "lines": [
                    "По Вриггерсу контактную задачу",
                    "удобно разделять на два уровня:",
                    "1. равновесие деформируемого тела;",
                    "2. выполнение контактных условий.",
                    r"$g_n\geq0,\quad p_n\geq0,\quad p_ng_n=0$",
                    "Активная область заранее неизвестна.",
                ],
            },
            {
                **common_card,
                "x": x_positions[1],
                "y": y,
                "title": "Внутренний цикл",
                "facecolor": "#f3f0ff",
                "edgecolor": "#b197fc",
                "lines": [
                    "При фиксированных множителях",
                    "и текущей активной области решается:",
                    r"$\mathbf{R}(\mathbf{u},\lambda^k)=\mathbf{f}_{int}+\mathbf{f}_c-\mathbf{f}_{ext}=0$",
                    r"$\mathbf{K}_T\Delta\mathbf{u}=-\mathbf{R}$",
                    r"$\mathbf{u}_{i+1}=\mathbf{u}_i+\alpha_i\Delta\mathbf{u}$",
                    r"$\|\mathbf{R}\|\leq tol_R$",
                    "Это Ньютоновский поиск равновесия шины.",
                ],
            },
            {
                **common_card,
                "x": x_positions[2],
                "y": y,
                "title": "Внешний цикл",
                "facecolor": "#fff4e6",
                "edgecolor": "#f0ad4e",
                "lines": [
                    "После найденного равновесия",
                    "уточняются контактные величины:",
                    r"$p_n^{k+1}=\langle\lambda_n^k-\rho g_n(\mathbf{u}^{*})\rangle_+$",
                    r"$A^{k+1}=\{q:\ p_n^{k+1}>0\}$",
                    r"$\lambda_n^{k+1}=\lambda_n^k+\omega(p_n^{k+1}-\lambda_n^k)$",
                    r"$\max\delta_n\leq tol_{\delta},\quad \|\Delta\lambda\|/\|\lambda\|\leq tol_{\lambda}$",
                    "Цикл останавливается, когда контакт стабилен.",
                ],
            },
        ],
        note="Внутренний цикл ищет равновесную форму шины; внешний цикл проверяет непроникновение и корректирует контактное давление.",
    )

    return outputs


def build_lagrange_result_assets(asset_dir: Path) -> dict[str, Path]:
    al_dir = REPO_ROOT / "results" / "main_scale_hyperelastic_reference_triplet_coarse_moving_plane" / "augmented_lagrangian_contact"
    surrogate_dir = REPO_ROOT / "results" / "main_scale_hyperelastic_reference_triplet_coarse_moving_plane" / "no_contact_surrogate"
    outputs = {
        "al_displacement_jacobian": asset_dir / "lagrange_displacement_jacobian.png",
        "al_stress_fields": asset_dir / "lagrange_stress_fields.png",
        "al_stress_profiles": asset_dir / "lagrange_stress_profiles.png",
        "surrogate_supplement": asset_dir / "surrogate_supplement_comparison.png",
    }

    build_image_pair_figure(
        outputs["al_displacement_jacobian"],
        al_dir / "displacement_magnitude.png",
        al_dir / "jacobian_determinant.png",
        left_title="Поле перемещений",
        right_title="Объемное изменение det(F)",
    )
    build_single_image_figure(
        outputs["al_stress_fields"],
        al_dir / "von_mises_stress.png",
    )
    build_contour_stress_profile_figure(
        outputs["al_stress_profiles"],
        al_dir / "ring_contour_stress_profiles.csv",
    )
    build_image_pair_figure(
        outputs["surrogate_supplement"],
        al_dir / "von_mises_stress.png",
        surrogate_dir / "von_mises_stress.png",
        left_title="Контактная постановка",
        right_title="Бесконтактная аппроксимация",
    )
    return outputs


def replace_text(root: ET.Element, replacements: dict[str, str]) -> None:
    for text_node in root.findall(".//a:t", NS):
        if text_node.text is None:
            continue
        text = text_node.text
        for old, new in replacements.items():
            text = text.replace(old, new)
        text_node.text = text


def set_shape_text_by_name(root: ET.Element, shape_name: str, text: str) -> None:
    for shape in root.findall(".//p:sp", NS):
        c_nv_pr = shape.find(".//p:cNvPr", NS)
        if c_nv_pr is None or c_nv_pr.attrib.get("name") != shape_name:
            continue
        text_nodes = shape.findall(".//a:t", NS)
        if not text_nodes:
            return
        text_nodes[0].text = text
        for node in text_nodes[1:]:
            node.text = ""
        return


def next_shape_id(root: ET.Element) -> int:
    max_id = 0
    for c_nv_pr in root.findall(".//p:cNvPr", NS):
        try:
            max_id = max(max_id, int(c_nv_pr.attrib.get("id", "0")))
        except ValueError:
            continue
    return max_id + 1


def remove_shape_by_name(root: ET.Element, name: str) -> None:
    sp_tree = root.find(".//p:cSld/p:spTree", NS)
    if sp_tree is None:
        return
    for child in list(sp_tree):
        c_nv_pr = child.find(".//p:cNvPr", NS)
        if c_nv_pr is not None and c_nv_pr.attrib.get("name") == name:
            sp_tree.remove(child)


def add_text_box(
    root: ET.Element,
    *,
    shape_id: int,
    name: str,
    x: int,
    y: int,
    cx: int,
    cy: int,
    paragraphs: list[str],
    font_size_pt: int = 18,
) -> None:
    sp_tree = root.find(".//p:cSld/p:spTree", NS)
    if sp_tree is None:
        raise RuntimeError("Could not find slide shape tree")

    shape = ET.SubElement(sp_tree, qname(P_NS, "sp"))
    nv_sp_pr = ET.SubElement(shape, qname(P_NS, "nvSpPr"))
    ET.SubElement(nv_sp_pr, qname(P_NS, "cNvPr"), {"id": str(shape_id), "name": name})
    ET.SubElement(nv_sp_pr, qname(P_NS, "cNvSpPr"), {"txBox": "1"})
    ET.SubElement(nv_sp_pr, qname(P_NS, "nvPr"))

    sp_pr = ET.SubElement(shape, qname(P_NS, "spPr"))
    xfrm = ET.SubElement(sp_pr, qname(A_NS, "xfrm"))
    ET.SubElement(xfrm, qname(A_NS, "off"), {"x": str(x), "y": str(y)})
    ET.SubElement(xfrm, qname(A_NS, "ext"), {"cx": str(cx), "cy": str(cy)})
    prst_geom = ET.SubElement(sp_pr, qname(A_NS, "prstGeom"), {"prst": "rect"})
    ET.SubElement(prst_geom, qname(A_NS, "avLst"))

    tx_body = ET.SubElement(shape, qname(P_NS, "txBody"))
    ET.SubElement(tx_body, qname(A_NS, "bodyPr"), {
        "wrap": "square",
        "rtlCol": "0",
        "anchor": "t",
    })
    ET.SubElement(tx_body, qname(A_NS, "lstStyle"))
    for paragraph in paragraphs:
        p = ET.SubElement(tx_body, qname(A_NS, "p"))
        p_pr = ET.SubElement(p, qname(A_NS, "pPr"), {
            "algn": "just",
        })
        ET.SubElement(p_pr, qname(A_NS, "spcAft")).append(
            ET.Element(qname(A_NS, "spcPts"), {"val": "850"})
        )
        run = ET.SubElement(p, qname(A_NS, "r"))
        r_pr = ET.SubElement(run, qname(A_NS, "rPr"), {
            "lang": "ru-RU",
            "sz": str(font_size_pt * 100),
        })
        ET.SubElement(r_pr, qname(A_NS, "latin"), {"typeface": "Times New Roman"})
        ET.SubElement(r_pr, qname(A_NS, "cs"), {"typeface": "Times New Roman"})
        text_node = ET.SubElement(run, qname(A_NS, "t"))
        text_node.text = paragraph


def prepare_normal_gap_slide(root: ET.Element) -> None:
    set_slide_image(root, "rId1")
    remove_shape_by_name(root, "Normal Gap Source")
    add_text_box(
        root,
        shape_id=next_shape_id(root),
        name="Normal Gap Source",
        x=502920,
        y=5372100,
        cx=6583680,
        cy=275000,
        paragraphs=[
            "Источник рисунка: S. Yamaguchi, Y. Sugawara, M. Takeda, Discover Applied Sciences, 2021, Fig. 2, CC BY 4.0."
        ],
        font_size_pt=9,
    )


def make_conclusion_slide(root: ET.Element) -> None:
    set_shape_text_by_name(root, "Text 0", "Заключение")
    sp_tree = root.find(".//p:cSld/p:spTree", NS)
    if sp_tree is None:
        raise RuntimeError("Could not find slide shape tree")

    keep_shapes = {"Text 0", "Text 3", "Text 4", "Shape 1", "Shape 2"}
    for child in list(sp_tree):
        if child.tag != qname(P_NS, "sp"):
            continue
        c_nv_pr = child.find(".//p:cNvPr", NS)
        shape_name = c_nv_pr.attrib.get("name", "") if c_nv_pr is not None else ""
        if shape_name not in keep_shapes:
            sp_tree.remove(child)

    conclusion_paragraphs = [
        "В работе рассмотрена контактная задача массивной шины с жесткой плоской опорной поверхностью в конечно-деформационной постановке. Материал шины задан сжимаемой неогуковской моделью в почти несжимаемом режиме; контакт считается гладким, без трения, с неизвестной заранее активной областью.",
        "Сравнение штрафного метода и модифицированного метода Лагранжа показало, что оба подхода дают близкую интегральную реакцию и близкую форму пятна контакта. Основное различие проявляется в качестве выполнения условия непроникновения: штрафной метод дешевле и устойчивее как инженерный расчет, а модифицированный метод Лагранжа существенно снижает проникание ценой большего числа итераций.",
        "Для инженерной интерпретации контакта наиболее информативны не только суммарная реакция, но и распределение контактного давления по дуге, активная зона и поля напряжений в шине. Бесконтактная аппроксимация восстановленной нагрузкой подтверждает близость общей картины деформированно-напряженного состояния, но не заменяет контактную постановку.",
        "Полученная модель является базовой расчетной схемой для дальнейшего развития: учета трения, смешанной u/p-постановки для почти несжимаемости, более сложных резиноподобных материалов и перехода к пространственной модели шины.",
    ]
    add_text_box(
        root,
        shape_id=next_shape_id(root),
        name="Conclusion Text",
        x=731520,
        y=1005840,
        cx=10668000,
        cy=5166360,
        paragraphs=conclusion_paragraphs,
        font_size_pt=18,
    )


def set_footer_number(root: ET.Element, number: int) -> None:
    for shape in root.findall(".//p:sp", NS):
        off = shape.find(".//a:xfrm/a:off", NS)
        if off is None:
            continue
        x = int(off.attrib.get("x", "0"))
        y = int(off.attrib.get("y", "0"))
        if x < 10_000_000 or y < 6_000_000:
            continue
        text_nodes = shape.findall(".//a:t", NS)
        if text_nodes:
            text_nodes[0].text = str(number)
            for node in text_nodes[1:]:
                node.text = ""


def set_slide_image(root: ET.Element, rel_id: str = "rId1") -> None:
    blip = root.find(".//p:pic//a:blip", NS)
    if blip is not None:
        blip.attrib[qname(R_NS, "embed")] = rel_id


def set_picture_geometry(root: ET.Element, *, x: int, y: int, cx: int, cy: int) -> None:
    picture = root.find(".//p:pic", NS)
    if picture is None:
        return
    off = picture.find(".//a:xfrm/a:off", NS)
    ext = picture.find(".//a:xfrm/a:ext", NS)
    if off is not None:
        off.attrib["x"] = str(x)
        off.attrib["y"] = str(y)
    if ext is not None:
        ext.attrib["cx"] = str(cx)
        ext.attrib["cy"] = str(cy)


def xml_bytes(root: ET.Element) -> bytes:
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def rels_bytes(items: list[tuple[str, str, str]]) -> bytes:
    root = ET.Element(qname(PKG_REL_NS, "Relationships"))
    for rel_id, rel_type, target in items:
        ET.SubElement(
            root,
            qname(PKG_REL_NS, "Relationship"),
            {"Id": rel_id, "Type": rel_type, "Target": target},
        )
    return xml_bytes(root)


def replace_relationship_target(xml_data: bytes, rel_id: str, target: str) -> bytes:
    root = ET.fromstring(xml_data)
    for rel in root:
        if rel.attrib.get("Id") == rel_id:
            rel.attrib["Target"] = target
            break
    return xml_bytes(root)


def next_relationship_id(existing_ids: set[str]) -> str:
    max_num = 0
    for rel_id in existing_ids:
        match = re.fullmatch(r"rId(\d+)", rel_id)
        if match:
            max_num = max(max_num, int(match.group(1)))
    candidate = max_num + 1
    while f"rId{candidate}" in existing_ids:
        candidate += 1
    existing_ids.add(f"rId{candidate}")
    return f"rId{candidate}"


def update_app_slide_count(xml_data: bytes, count: int) -> bytes:
    root = ET.fromstring(xml_data)
    for child in root:
        if child.tag.endswith("Slides"):
            child.text = str(count)
            return xml_bytes(root)
    return xml_data


def update_presentation(
    pptx_path: Path,
    parameter_study_image: Path,
    final_contact_image: Path,
    formula_assets: dict[str, Path],
    lagrange_assets: dict[str, Path],
) -> Path:
    register_namespaces()
    replacements = {
        "Расширенный метод Лагранжа": "Модифицированный метод Лагранжа",
        "Расширенный Лагранж": "Модифицированный метод Лагранжа",
        "Расш. Лагранж": "Модифицированный метод Лагранжа",
        "Метод расширенного Лагранжа": "Модифицированный метод Лагранжа",
        "метод расширенного Лагранжа": "модифицированный метод Лагранжа",
        "метода расширенного Лагранжа": "модифицированного метода Лагранжа",
        "расширенный Лагранж": "модифицированный метод Лагранжа",
        "метод РЛ": "модифицированный метод Лагранжа",
        "метода РЛ": "модифицированного метода Лагранжа",
        "РЛ": "модифицированного метода Лагранжа",
        "ускорение к модифицированного метода Лагранжа": "ускорение суррогатной постановки",
        "давление используется для калибровки": "давление задает локальное нагружение",
        "Модифицированный метод Лагранжа выбран как база для восстановления эквивалентной нагрузки, потому что точнее выполняет непроникновение.": "Модифицированный метод Лагранжа выбран как уточняющая постановка, потому что точнее выполняет непроникновение.",
        "μ | 3,0 МПа": "μ | 4,0 МПа",
        "3,0 МПа": "4,0 МПа",
    }

    new_slide_specs = [
        {
            "file": "ppt/slides/slide18.xml",
            "rels": "ppt/slides/_rels/slide18.xml.rels",
            "media": "ppt/media/image-formula-neo-hookean.png",
            "source_image": formula_assets["neo_hookean"],
            "title": "Неогуковская модель материала",
            "caption": "Материал описан через энергию конечных деформаций; параметры заданы модулем сдвига μ и объемным модулем K.",
        },
        {
            "file": "ppt/slides/slide19.xml",
            "rels": "ppt/slides/_rels/slide19.xml.rels",
            "media": "ppt/media/image-formula-penalty-method.png",
            "source_image": formula_assets["penalty"],
            "title": "Штрафной метод контактного взаимодействия",
            "caption": "Контактная реакция пропорциональна прониканию; метод прост и устойчив, но допускает конечное нарушение условия непроникновения.",
        },
        {
            "file": "ppt/slides/slide20.xml",
            "rels": "ppt/slides/_rels/slide20.xml.rels",
            "media": "ppt/media/image-formula-augmented-lagrangian.png",
            "source_image": formula_assets["augmented_lagrangian"],
            "title": "Модифицированный метод Лагранжа",
            "caption": "Контактное давление уточняется через множители Лагранжа и штрафную стабилизацию, поэтому проникание подавляется значительно сильнее.",
        },
        {
            "file": "ppt/slides/slide21.xml",
            "rels": "ppt/slides/_rels/slide21.xml.rels",
            "media": "ppt/media/image-formula-inner-outer-loop.png",
            "source_image": formula_assets["inner_outer_loop"],
            "title": "Внешний и внутренний циклы контактного решения",
            "caption": "Внутренний цикл ищет равновесие тела, внешний цикл уточняет активную область контакта и контактные множители.",
        },
        {
            "file": "ppt/slides/slide12.xml",
            "rels": "ppt/slides/_rels/slide12.xml.rels",
            "media": "ppt/media/image-contact-parameter-study.png",
            "source_image": parameter_study_image,
            "title": "Чувствительность методов к параметру и сетке",
            "caption": "Серия расчетов показывает компромисс: штрафной метод проще и устойчивее, модифицированный метод Лагранжа лучше подавляет проникание, но требует более дорогой итерационной настройки.",
        },
        {
            "file": "ppt/slides/slide13.xml",
            "rels": "ppt/slides/_rels/slide13.xml.rels",
            "media": "ppt/media/image-final-contact-comparison.png",
            "source_image": final_contact_image,
            "title": "Финальное сравнение контактных постановок",
            "caption": "Интегральная реакция и пятно контакта близки; главное отличие проявляется в контроле проникания и вычислительной цене метода.",
        },
        {
            "file": "ppt/slides/slide14.xml",
            "rels": "ppt/slides/_rels/slide14.xml.rels",
            "media": "ppt/media/image-lagrange-displacement-jacobian.png",
            "source_image": lagrange_assets["al_displacement_jacobian"],
            "title": "Деформированное состояние: модифицированный метод Лагранжа",
            "caption": "Расчет показывает конечные деформации шины в зоне контакта; det(F) остается близким к единице, что соответствует почти несжимаемому материалу.",
        },
        {
            "file": "ppt/slides/slide15.xml",
            "rels": "ppt/slides/_rels/slide15.xml.rels",
            "media": "ppt/media/image-lagrange-stress-fields.png",
            "source_image": lagrange_assets["al_stress_fields"],
            "title": "Напряженное состояние: модифицированный метод Лагранжа",
            "caption": "Максимальные напряжения локализуются в зоне контакта и переходных участках кольцевого сектора; картина соответствует изгибу и локальному сжатию шины.",
        },
        {
            "file": "ppt/slides/slide16.xml",
            "rels": "ppt/slides/_rels/slide16.xml.rels",
            "media": "ppt/media/image-lagrange-stress-profiles.png",
            "source_image": lagrange_assets["al_stress_profiles"],
            "title": "Профили напряжений в сечениях шины",
            "caption": "Профили по внутреннему, среднему и внешнему контурам показывают, как компоненты напряжений меняются вдоль зоны контакта.",
        },
        {
            "file": "ppt/slides/slide17.xml",
            "rels": "ppt/slides/_rels/slide17.xml.rels",
            "media": "ppt/media/image-surrogate-supplement.png",
            "source_image": lagrange_assets["surrogate_supplement"],
            "title": "Бесконтактная аппроксимация как дополнительная проверка",
            "caption": "Восстановленная нагрузка дает похожую общую картину деформированно-напряженного состояния, но не контролирует зазор и не заменяет контактную постановку.",
        },
    ]
    generated_slide_files = {spec["file"] for spec in new_slide_specs}
    generated_slide_targets = {spec["file"].replace("ppt/", "") for spec in new_slide_specs}

    temp_path = DOCS_DIR / "_presentation_update_tmp.pptx"
    if temp_path.exists():
        temp_path.unlink()
    try:
        with ZipFile(pptx_path, "r") as zin:
            original_files = {name: zin.read(name) for name in zin.namelist()}

        slide_names = sorted(
            [name for name in original_files if re.fullmatch(r"ppt/slides/slide\d+\.xml", name)],
            key=lambda name: int(re.search(r"slide(\d+)\.xml", name).group(1)),
        )
        excluded_original_slides = {
            "ppt/slides/slide3.xml",
            "ppt/slides/slide5.xml",
            "ppt/slides/slide7.xml",
            "ppt/slides/slide8.xml",
            "ppt/slides/slide9.xml",
            "ppt/slides/slide10.xml",
        }
        slide_names = [
            name
            for name in slide_names
            if name not in generated_slide_files and name not in excluded_original_slides
        ]
        original_by_name = {Path(name).name: name for name in slide_names}
        formula_slide_files = [spec["file"] for spec in new_slide_specs[:4]]
        comparison_slide_files = [spec["file"] for spec in new_slide_specs[4:6]]
        result_slide_files = [spec["file"] for spec in new_slide_specs[6:]]

        final_order = []
        for slide_name in ["slide1.xml", "slide2.xml", "slide4.xml"]:
            if slide_name in original_by_name:
                final_order.append(original_by_name[slide_name])
        final_order += formula_slide_files
        final_order += comparison_slide_files
        if "slide6.xml" in original_by_name:
            final_order.append(original_by_name["slide6.xml"])
        final_order += result_slide_files
        if "slide11.xml" in original_by_name:
            final_order.append(original_by_name["slide11.xml"])

        modified_files: dict[str, bytes] = {}

        visible_numbers = {slide_name: index for index, slide_name in enumerate(final_order, start=1)}
        for visible_index, slide_name in enumerate(final_order, start=1):
            if slide_name in original_files:
                root = ET.fromstring(original_files[slide_name])
                replace_text(root, replacements)
                set_footer_number(root, visible_index)
                if Path(slide_name).name == "slide4.xml":
                    prepare_normal_gap_slide(root)
                    slide_rels = "ppt/slides/_rels/slide4.xml.rels"
                    modified_files[slide_rels] = replace_relationship_target(
                        original_files[slide_rels],
                        "rId1",
                        f"../media/{Path(NORMAL_GAP_PRESENTATION_MEDIA).name}",
                    )
                    modified_files[NORMAL_GAP_PRESENTATION_MEDIA] = NORMAL_GAP_SOURCE_IMAGE.read_bytes()
                if Path(slide_name).name == "slide11.xml":
                    make_conclusion_slide(root)
                modified_files[slide_name] = xml_bytes(root)

        template_root = ET.fromstring(original_files["ppt/slides/slide8.xml"])
        for spec in new_slide_specs:
            root = ET.fromstring(ET.tostring(template_root, encoding="utf-8"))
            replace_text(root, replacements)
            set_shape_text_by_name(root, "Text 0", spec["title"])
            set_shape_text_by_name(root, "Text 6", spec["caption"])
            set_footer_number(root, visible_numbers[spec["file"]])
            set_slide_image(root, "rId1")
            if spec["file"] == "ppt/slides/slide16.xml":
                remove_shape_by_name(root, "Shape 5")
                remove_shape_by_name(root, "Text 6")
                set_picture_geometry(
                    root,
                    x=1516380,
                    y=914400,
                    cx=9144000,
                    cy=4908000,
                )
            modified_files[spec["file"]] = xml_bytes(root)
            modified_files[spec["rels"]] = rels_bytes(
                [
                    ("rId1", IMAGE_REL_TYPE, f"../media/{Path(spec['media']).name}"),
                    ("rId2", LAYOUT_REL_TYPE, "../slideLayouts/slideLayout1.xml"),
                ]
            )
            modified_files[spec["media"]] = Path(spec["source_image"]).read_bytes()

        presentation_root = ET.fromstring(original_files["ppt/presentation.xml"])
        presentation_rels_root = ET.fromstring(original_files["ppt/_rels/presentation.xml.rels"])
        original_presentation_root = ET.fromstring(original_files["ppt/presentation.xml"])
        original_rels_root = ET.fromstring(original_files["ppt/_rels/presentation.xml.rels"])
        original_rel_target_by_id = {
            rel.attrib["Id"]: rel.attrib.get("Target")
            for rel in original_rels_root
            if rel.attrib.get("Type") == SLIDE_REL_TYPE
        }
        original_slide_id_by_target: dict[str, int] = {}
        original_slide_id_list = original_presentation_root.find("p:sldIdLst", NS)
        if original_slide_id_list is not None:
            for element in original_slide_id_list:
                rel_id = element.attrib.get(qname(R_NS, "id"))
                target = original_rel_target_by_id.get(rel_id)
                if target is not None and target not in generated_slide_targets:
                    original_slide_id_by_target[target] = int(element.attrib["id"])

        for rel in list(presentation_rels_root):
            if rel.attrib.get("Type") == SLIDE_REL_TYPE and rel.attrib.get("Target") in generated_slide_targets:
                presentation_rels_root.remove(rel)

        existing_rel_ids = {rel.attrib["Id"] for rel in presentation_rels_root}
        slide_rel_by_target = {
            rel.attrib.get("Target"): rel.attrib.get("Id")
            for rel in presentation_rels_root
            if rel.attrib.get("Type") == SLIDE_REL_TYPE
        }

        for spec in new_slide_specs:
            target = spec["file"].replace("ppt/", "")
            rel_id = next_relationship_id(existing_rel_ids)
            ET.SubElement(
                presentation_rels_root,
                qname(PKG_REL_NS, "Relationship"),
                {"Id": rel_id, "Type": SLIDE_REL_TYPE, "Target": target},
            )
            slide_rel_by_target[target] = rel_id

        slide_id_list = presentation_root.find("p:sldIdLst", NS)
        if slide_id_list is None:
            raise RuntimeError("Could not find p:sldIdLst in presentation.xml")
        for child in list(slide_id_list):
            slide_id_list.remove(child)
        existing_slide_ids = list(original_slide_id_by_target.values())
        next_slide_id = max(existing_slide_ids) + 1
        for slide_name in final_order:
            target = slide_name.replace("ppt/", "")
            rel_id = slide_rel_by_target[target]
            slide_id = next_slide_id if slide_name in generated_slide_files else original_slide_id_by_target[target]
            if slide_name in generated_slide_files:
                next_slide_id += 1
            ET.SubElement(
                slide_id_list,
                qname(P_NS, "sldId"),
                {"id": str(slide_id), qname(R_NS, "id"): rel_id},
            )

        modified_files["ppt/presentation.xml"] = xml_bytes(presentation_root)
        modified_files["ppt/_rels/presentation.xml.rels"] = xml_bytes(presentation_rels_root)

        content_types_root = ET.fromstring(original_files["[Content_Types].xml"])
        existing_overrides = {
            override.attrib.get("PartName")
            for override in content_types_root.findall("ct:Override", NS)
        }
        for spec in new_slide_specs:
            part_name = "/" + spec["file"]
            if part_name not in existing_overrides:
                ET.SubElement(
                    content_types_root,
                    qname(CT_NS, "Override"),
                    {"PartName": part_name, "ContentType": SLIDE_CONTENT_TYPE},
                )
        modified_files["[Content_Types].xml"] = xml_bytes(content_types_root)

        if "docProps/app.xml" in original_files:
            modified_files["docProps/app.xml"] = update_app_slide_count(
                original_files["docProps/app.xml"],
                len(final_order),
            )

        skip_files = set(modified_files)
        with ZipFile(temp_path, "w", ZIP_DEFLATED) as zout:
            for name, data in original_files.items():
                if name in skip_files:
                    continue
                zout.writestr(name, data)
            for name, data in modified_files.items():
                zout.writestr(name, data)

        try:
            shutil.copy2(temp_path, pptx_path)
            return pptx_path
        except PermissionError:
            fallback_path = pptx_path.with_name(f"{pptx_path.stem}_обновлено{pptx_path.suffix}")
            shutil.copy2(temp_path, fallback_path)
            return fallback_path
    finally:
        if temp_path.exists():
            temp_path.unlink()


def main() -> int:
    configure_matplotlib()
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    parameter_study_image = ASSET_DIR / "contact_method_parameter_study.png"
    final_contact_image = ASSET_DIR / "final_contact_method_comparison.png"
    build_parameter_study_figure(parameter_study_image)
    build_final_contact_figure(final_contact_image)
    formula_assets = build_formula_assets(ASSET_DIR)
    lagrange_assets = build_lagrange_result_assets(ASSET_DIR)
    updated_presentation_path = update_presentation(
        presentation_path(),
        parameter_study_image,
        final_contact_image,
        formula_assets,
        lagrange_assets,
    )
    print(f"Updated presentation: {updated_presentation_path}")
    print(f"Saved assets: {parameter_study_image}, {final_contact_image}")
    for asset_name, asset_path in formula_assets.items():
        print(f"Saved {asset_name}: {asset_path}")
    for asset_name, asset_path in lagrange_assets.items():
        print(f"Saved {asset_name}: {asset_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
