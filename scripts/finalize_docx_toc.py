#!/usr/bin/env python3
from __future__ import annotations

import io
import re
import shutil
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "docs" / "MironovMV_RK5-81B_moving_plane_results.docx"
LAYOUT = ROOT / "docs" / "layout_pages.tsv"

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
MC = "http://schemas.openxmlformats.org/markup-compatibility/2006"


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


def parse_layout() -> tuple[int, dict[str, int]]:
    total_pages = 0
    pages: dict[str, int] = {}
    heading_pattern = re.compile(r"^(?:ВВЕДЕНИЕ|[2-7]\.\s|[2-5]\.\d+\s)")
    for line in LAYOUT.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        left, _, right = line.partition("\t")
        if left == "TOTAL":
            total_pages = int(right)
            continue
        if not left.isdigit() or not heading_pattern.match(right):
            continue
        # The first block is the old static TOC itself; actual body headings occur later
        # and overwrite those entries in the map.
        pages[right] = int(left)
    if not total_pages:
        raise RuntimeError("TOTAL page count not found in layout_pages.tsv")
    return total_pages, pages


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
        if text == "7. СПИСОК ИСПОЛЬЗОВАННЫХ ИСТОЧНИКОВ":
            in_sources = True
            continue
        if in_sources and re.match(r"^\d+\.\s+", text):
            total += 1
    return total


def toc_title(text: str) -> str | None:
    match = re.match(r"^(.*?)(\d+)$", text)
    if not match:
        return None
    title = match.group(1).strip()
    if title == "СОДЕРЖАНИЕ":
        return None
    return title


def update_docx() -> tuple[int, int, int, int, int]:
    total_pages, page_map = parse_layout()
    with zipfile.ZipFile(DOCX, "r") as zin:
        document_bytes = zin.read("word/document.xml")
        register_namespaces_from_xml(document_bytes)
        root = ET.fromstring(document_bytes)
        body = root.find(q(W, "body"))
        if body is None:
            raise RuntimeError("word/document.xml has no w:body")
        body_children = list(body)

        updated_toc = 0
        in_toc = False
        for child in body_children:
            text = element_text(child)
            if text == "СОДЕРЖАНИЕ":
                in_toc = True
                continue
            if in_toc and text == "ВВЕДЕНИЕ":
                break
            if not in_toc:
                continue
            title = toc_title(text)
            if title is None:
                continue
            page = page_map.get(title)
            if page is None:
                raise RuntimeError(f"TOC title has no page: {title}")
            set_element_text(child, f"{title}{page}")
            updated_toc += 1

        figures = count_captions(body_children, "Рисунок ")
        tables = count_captions(body_children, "Таблица ")
        sources = count_sources(body_children)
        set_element_text(
            body_children[2],
            (
                "Выпускная квалификационная работа содержит "
                f"{total_pages} страниц, {figures} рисунка, {tables} таблиц, "
                f"{sources} использованных источников."
            ),
        )
        if root.get(q(MC, "Ignorable")):
            root.set(q(MC, "Ignorable"), "w14 wp14")

        document_xml = ET.tostring(root, encoding="utf-8", xml_declaration=True)
        tmp = DOCX.with_suffix(".tmp.docx")
        if tmp.exists():
            tmp.unlink()
        with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename == "word/document.xml":
                    continue
                zout.writestr(item, zin.read(item.filename))
            zout.writestr("word/document.xml", document_xml)

    shutil.move(str(tmp), str(DOCX))
    return updated_toc, total_pages, figures, tables, sources


def main() -> None:
    updated_toc, pages, figures, tables, sources = update_docx()
    with zipfile.ZipFile(DOCX, "r") as archive:
        bad_file = archive.testzip()
    if bad_file is not None:
        raise RuntimeError(f"Bad file in generated docx: {bad_file}")
    print(f"updated_toc_entries={updated_toc}")
    print(f"abstract_counts=pages:{pages},figures:{figures},tables:{tables},sources:{sources}")
    print("zip_check=ok")


if __name__ == "__main__":
    main()
