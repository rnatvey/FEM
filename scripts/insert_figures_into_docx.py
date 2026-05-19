#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
FIGURES = DOCS / "figures_numbered"
MANIFEST = FIGURES / "manifest.csv"

EMU_PER_CM = 360_000
MAX_WIDTH_EMU = int(15.6 * EMU_PER_CM)
MAX_HEIGHT_EMU = int(10.6 * EMU_PER_CM)

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
PIC = "http://schemas.openxmlformats.org/drawingml/2006/picture"
REL = "http://schemas.openxmlformats.org/package/2006/relationships"
CT = "http://schemas.openxmlformats.org/package/2006/content-types"


for prefix, uri in {
    "w": W,
    "r": R,
    "wp": WP,
    "a": A,
    "pic": PIC,
    "rel": REL,
    "ct": CT,
}.items():
    ET.register_namespace(prefix, uri)


def q(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


@dataclass(frozen=True)
class Figure:
    number: str
    image: Path
    caption: str


def load_manifest() -> dict[str, Figure]:
    figures: dict[str, Figure] = {}
    with MANIFEST.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            number = row["number"]
            image = FIGURES / row["file"]
            if not image.exists():
                raise FileNotFoundError(image)
            figures[number] = Figure(number=number, image=image, caption=row["caption"])
    return figures


def paragraph_text(paragraph: ET.Element) -> str:
    return "".join((node.text or "") for node in paragraph.findall(f".//{q(W, 't')}")).strip()


def element_text(element: ET.Element) -> str:
    return "".join((node.text or "") for node in element.findall(f".//{q(W, 't')}")).strip()


def remove_placeholder_tables(body: ET.Element) -> int:
    placeholder = "\u041c\u0435\u0441\u0442\u043e \u0434\u043b\u044f \u0440\u0438\u0441\u0443\u043d\u043a\u0430"
    removed = 0
    for child in list(body):
        if child.tag == q(W, "tbl") and placeholder in element_text(child):
            body.remove(child)
            removed += 1
    return removed


def image_size_emu(path: Path) -> tuple[int, int]:
    with Image.open(path) as img:
        width, height = img.size
    scale = min(MAX_WIDTH_EMU / width, MAX_HEIGHT_EMU / height)
    return int(width * scale), int(height * scale)


def image_paragraph(rel_id: str, file_name: str, number: str, cx: int, cy: int, docpr_id: int) -> ET.Element:
    p = ET.Element(q(W, "p"))
    p_pr = ET.SubElement(p, q(W, "pPr"))
    jc = ET.SubElement(p_pr, q(W, "jc"))
    jc.set(q(W, "val"), "center")
    spacing = ET.SubElement(p_pr, q(W, "spacing"))
    spacing.set(q(W, "before"), "120")
    spacing.set(q(W, "after"), "120")

    run = ET.SubElement(p, q(W, "r"))
    drawing = ET.SubElement(run, q(W, "drawing"))
    inline = ET.SubElement(drawing, q(WP, "inline"))
    inline.set("distT", "0")
    inline.set("distB", "0")
    inline.set("distL", "0")
    inline.set("distR", "0")

    extent = ET.SubElement(inline, q(WP, "extent"))
    extent.set("cx", str(cx))
    extent.set("cy", str(cy))

    effect_extent = ET.SubElement(inline, q(WP, "effectExtent"))
    for attr in ("l", "t", "r", "b"):
        effect_extent.set(attr, "0")

    doc_pr = ET.SubElement(inline, q(WP, "docPr"))
    doc_pr.set("id", str(docpr_id))
    doc_pr.set("name", f"Рисунок {number}")

    c_nv = ET.SubElement(inline, q(WP, "cNvGraphicFramePr"))
    locks = ET.SubElement(c_nv, q(A, "graphicFrameLocks"))
    locks.set("noChangeAspect", "1")

    graphic = ET.SubElement(inline, q(A, "graphic"))
    graphic_data = ET.SubElement(graphic, q(A, "graphicData"))
    graphic_data.set("uri", PIC)
    pic = ET.SubElement(graphic_data, q(PIC, "pic"))

    nv_pic_pr = ET.SubElement(pic, q(PIC, "nvPicPr"))
    c_nv_pr = ET.SubElement(nv_pic_pr, q(PIC, "cNvPr"))
    c_nv_pr.set("id", "0")
    c_nv_pr.set("name", file_name)
    c_nv_pic_pr = ET.SubElement(nv_pic_pr, q(PIC, "cNvPicPr"))
    pic_locks = ET.SubElement(c_nv_pic_pr, q(A, "picLocks"))
    pic_locks.set("noChangeAspect", "1")
    pic_locks.set("noChangeArrowheads", "1")

    blip_fill = ET.SubElement(pic, q(PIC, "blipFill"))
    blip = ET.SubElement(blip_fill, q(A, "blip"))
    blip.set(q(R, "embed"), rel_id)
    stretch = ET.SubElement(blip_fill, q(A, "stretch"))
    ET.SubElement(stretch, q(A, "fillRect"))

    sp_pr = ET.SubElement(pic, q(PIC, "spPr"))
    xfrm = ET.SubElement(sp_pr, q(A, "xfrm"))
    off = ET.SubElement(xfrm, q(A, "off"))
    off.set("x", "0")
    off.set("y", "0")
    ext = ET.SubElement(xfrm, q(A, "ext"))
    ext.set("cx", str(cx))
    ext.set("cy", str(cy))
    prst_geom = ET.SubElement(sp_pr, q(A, "prstGeom"))
    prst_geom.set("prst", "rect")
    ET.SubElement(prst_geom, q(A, "avLst"))

    return p


def next_relationship_id(rels_root: ET.Element) -> int:
    max_id = 0
    for rel in rels_root.findall(q(REL, "Relationship")):
        rel_id = rel.get("Id", "")
        match = re.fullmatch(r"rId(\d+)", rel_id)
        if match:
            max_id = max(max_id, int(match.group(1)))
    return max_id + 1


def next_docpr_id(document_root: ET.Element) -> int:
    max_id = 0
    for doc_pr in document_root.findall(f".//{q(WP, 'docPr')}"):
        value = doc_pr.get("id")
        if value and value.isdigit():
            max_id = max(max_id, int(value))
    return max_id + 1


def ensure_png_content_type(content_types_root: ET.Element) -> None:
    for item in content_types_root.findall(q(CT, "Default")):
        if item.get("Extension") == "png":
            return
    default = ET.SubElement(content_types_root, q(CT, "Default"))
    default.set("Extension", "png")
    default.set("ContentType", "image/png")


def find_source_docx() -> Path:
    candidates = [
        path
        for path in DOCS.glob("*.docx")
        if not path.name.startswith("~$")
        and "_с_рисунками" not in path.stem
        and "_with_figures" not in path.stem
    ]
    if not candidates:
        raise FileNotFoundError("No source .docx found in docs")
    valid: list[Path] = []
    for path in candidates:
        try:
            with zipfile.ZipFile(path) as archive:
                archive.getinfo("word/document.xml")
        except (KeyError, zipfile.BadZipFile):
            continue
        valid.append(path)
    if not valid:
        raise FileNotFoundError("No readable source .docx found in docs")
    return max(valid, key=lambda path: path.stat().st_mtime)


def build_output_path(source: Path) -> Path:
    return source.with_name(f"{source.stem}_с_рисунками{source.suffix}")


def insert_figures(source: Path, output: Path) -> tuple[list[str], int]:
    figures = load_manifest()
    caption_prefix = "\u0420\u0438\u0441\u0443\u043d\u043e\u043a"
    inserted: list[str] = []

    with zipfile.ZipFile(source, "r") as zin:
        document_root = ET.fromstring(zin.read("word/document.xml"))
        rels_root = ET.fromstring(zin.read("word/_rels/document.xml.rels"))
        content_types_root = ET.fromstring(zin.read("[Content_Types].xml"))
        body = document_root.find(q(W, "body"))
        if body is None:
            raise RuntimeError("word/document.xml has no w:body")
        removed_placeholders = remove_placeholder_tables(body)

        rel_counter = next_relationship_id(rels_root)
        docpr_counter = next_docpr_id(document_root)
        image_rels: dict[str, tuple[str, Path]] = {}

        body_children = list(body)
        for index, child in reversed(list(enumerate(body_children))):
            if child.tag != q(W, "p"):
                continue
            text = paragraph_text(child)
            if not text.startswith(caption_prefix):
                continue
            match = re.search(r"(\d+\.\d+)", text)
            if not match:
                continue
            number = match.group(1)
            figure = figures.get(number)
            if figure is None:
                raise KeyError(f"No figure for caption {number}")

            rel_id = f"rId{rel_counter}"
            rel_counter += 1
            media_name = f"fig_{number.replace('.', '_')}.png"
            image_rels[media_name] = (rel_id, figure.image)

            rel = ET.SubElement(rels_root, q(REL, "Relationship"))
            rel.set("Id", rel_id)
            rel.set("Type", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image")
            rel.set("Target", f"media/{media_name}")

            cx, cy = image_size_emu(figure.image)
            body.insert(index, image_paragraph(rel_id, media_name, number, cx, cy, docpr_counter))
            docpr_counter += 1
            inserted.append(number)

        inserted.reverse()
        missing = sorted(set(figures) - set(inserted), key=lambda item: tuple(map(int, item.split("."))))
        if missing:
            raise RuntimeError(f"Captions not found for figures: {', '.join(missing)}")

        ensure_png_content_type(content_types_root)

        document_xml = ET.tostring(document_root, encoding="utf-8", xml_declaration=True)
        rels_xml = ET.tostring(rels_root, encoding="utf-8", xml_declaration=True)
        content_types_xml = ET.tostring(content_types_root, encoding="utf-8", xml_declaration=True)

        tmp = output.with_suffix(".tmp.docx")
        if tmp.exists():
            tmp.unlink()
        with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            skip = {"word/document.xml", "word/_rels/document.xml.rels", "[Content_Types].xml"}
            for item in zin.infolist():
                if item.filename in skip:
                    continue
                zout.writestr(item, zin.read(item.filename))
            zout.writestr("word/document.xml", document_xml)
            zout.writestr("word/_rels/document.xml.rels", rels_xml)
            zout.writestr("[Content_Types].xml", content_types_xml)
            for media_name, (_, image_path) in image_rels.items():
                zout.writestr(f"word/media/{media_name}", image_path.read_bytes())

    if output.exists():
        output.unlink()
    shutil.move(str(tmp), str(output))
    return inserted, removed_placeholders


def main() -> None:
    source = find_source_docx()
    output = build_output_path(source)
    inserted, removed_placeholders = insert_figures(source, output)
    print(f"source={source}")
    print(f"output={output}")
    print(f"inserted={len(inserted)}")
    print(f"removed_placeholders={removed_placeholders}")


if __name__ == "__main__":
    main()
