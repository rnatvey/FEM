#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "main_scale_hyperelastic_reference_triplet_coarse_symmetric_anchor"
OUT = ROOT / "docs" / "figures_numbered"

# GOST-like drafting palette: mostly monochrome, with muted accents only where
# they separate boundary conditions, contact state, or comparison curves.
INK = "#202020"
BLUE = "#2f5f8f"
RED = "#8f2f2f"
GREEN = "#3f6f4f"
ORANGE = "#8a6630"
GRAY = "#5c5c5c"
LIGHT_GRAY = "#f6f7f8"
MID_GRAY = "#bfc5cc"
TEAL = "#4f7378"

OBJECT_LW = 1.05
ACCENT_LW = 1.55
THIN_LW = 0.45
DIM_LW = 0.75
ARROW_SCALE = 8

R_INNER_MM = 250.0
R_OUTER_MM = 300.0
THETA1_DEG = 210.0
THETA2_DEG = 330.0


@dataclass(frozen=True)
class FigureTask:
    number: str
    caption: str
    builder: Callable[[Path], str]


MANIFEST: list[dict[str, str]] = []


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "normal",
            "axes.labelsize": 9,
            "axes.edgecolor": INK,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "grid.color": "#d8dde3",
            "grid.linewidth": 0.45,
            "grid.linestyle": "--",
            "legend.fontsize": 8,
            "legend.frameon": False,
            "mathtext.fontset": "dejavusans",
            "savefig.dpi": 300,
            "svg.fonttype": "none",
        }
    )


def file_name(number: str) -> str:
    chapter, item = number.split(".")
    return f"fig_{chapter}_{int(item):02d}.png"


def save(fig: plt.Figure, path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return "generated schematic/chart; editable SVG saved"


def copy_image(src: Path) -> Callable[[Path], str]:
    def builder(dst: Path) -> str:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        return str(src.relative_to(ROOT))

    return builder


def read_numeric_csv(path: Path) -> list[dict[str, float | str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows: list[dict[str, float | str]] = []
        for row in reader:
            parsed: dict[str, float | str] = {}
            for key, value in row.items():
                if value is None:
                    parsed[key] = ""
                    continue
                try:
                    parsed[key] = float(value)
                except ValueError:
                    parsed[key] = value
            rows.append(parsed)
        return rows


def read_triplet_summary() -> tuple[dict[str, dict[str, float | str]], dict[str, float]]:
    cases: dict[str, dict[str, float | str]] = {}
    params: dict[str, float] = {}
    with (RESULTS / "triplet_summary.csv").open("r", encoding="utf-8-sig", newline="") as handle:
        raw_rows = list(csv.reader(handle))

    section = "cases"
    headers: list[str] = []
    for row in raw_rows:
        if not row:
            section = "params"
            headers = []
            continue
        if not headers:
            headers = row
            continue
        if section == "cases":
            item: dict[str, float | str] = {}
            for key, value in zip(headers, row):
                try:
                    item[key] = float(value)
                except ValueError:
                    item[key] = value
            cases[str(item["case"])] = item
        else:
            if len(row) >= 2:
                try:
                    params[row[0]] = float(row[1])
                except ValueError:
                    pass
    return cases, params


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = GRAY, lw: float = DIM_LW, **kwargs) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(arrowstyle="->", color=color, lw=lw, shrinkA=0, shrinkB=0, mutation_scale=ARROW_SCALE, **kwargs),
    )


def double_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = GRAY, lw: float = DIM_LW) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(arrowstyle="<->", color=color, lw=lw, shrinkA=0, shrinkB=0, mutation_scale=ARROW_SCALE),
    )


def setup_plain(ax: plt.Axes, xlim: tuple[float, float] = (0, 10), ylim: tuple[float, float] = (0, 6)) -> None:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def box(ax: plt.Axes, xy: tuple[float, float], width: float, height: float, text: str, fc: str = "white", ec: str = BLUE) -> None:
    patch = patches.Rectangle(
        xy,
        width,
        height,
        facecolor=fc,
        edgecolor=ec,
        linewidth=OBJECT_LW,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=8.7)


def ring_sector(ax: plt.Axes, center: tuple[float, float] = (0, 0), r_outer: float = 3.0, r_inner: float = 2.5,
                theta1: float = 210, theta2: float = 330, fc: str = LIGHT_GRAY, ec: str = INK, lw: float = OBJECT_LW) -> None:
    wedge = patches.Wedge(center, r_outer, theta1, theta2, width=r_outer - r_inner, facecolor=fc, edgecolor=ec, linewidth=lw)
    ax.add_patch(wedge)


def outer_arc_points(theta1: float, theta2: float, radius: float = 3.0, center: tuple[float, float] = (0, 0), n: int = 80) -> tuple[np.ndarray, np.ndarray]:
    th = np.deg2rad(np.linspace(theta1, theta2, n))
    return center[0] + radius * np.cos(th), center[1] + radius * np.sin(th)


def ring_arc_mm(theta1: float, theta2: float, radius: float = R_OUTER_MM, n: int = 160) -> tuple[np.ndarray, np.ndarray]:
    th = np.deg2rad(np.linspace(theta1, theta2, n))
    return radius * np.cos(th), radius * np.sin(th)


def add_ring_mesh(ax: plt.Axes, center: tuple[float, float] = (0, 0), theta1: float = 210, theta2: float = 330) -> None:
    for r in np.linspace(2.5, 3.0, 5):
        x, y = outer_arc_points(theta1, theta2, r, center, 100)
        ax.plot(x, y, color=MID_GRAY, lw=THIN_LW)
    for th in np.linspace(theta1, theta2, 13):
        rad = math.radians(th)
        ax.plot(
            [center[0] + 2.5 * math.cos(rad), center[0] + 3.0 * math.cos(rad)],
            [center[1] + 2.5 * math.sin(rad), center[1] + 3.0 * math.sin(rad)],
            color=MID_GRAY,
            lw=THIN_LW,
        )


def draw_results_ring_mesh(ax: plt.Axes) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-170, 170)
    ax.set_ylim(-315, -205)
    ax.axis("off")
    for radius in np.linspace(R_INNER_MM, R_OUTER_MM, 6):
        x, y = ring_arc_mm(THETA1_DEG, THETA2_DEG, radius, 180)
        ax.plot(x, y, color=MID_GRAY, lw=THIN_LW)
    for angle in np.linspace(THETA1_DEG, THETA2_DEG, 17):
        xi, yi = ring_arc_mm(angle, angle, R_INNER_MM, 1)
        xo, yo = ring_arc_mm(angle, angle, R_OUTER_MM, 1)
        ax.plot([xi[0], xo[0]], [yi[0], yo[0]], color=MID_GRAY, lw=THIN_LW)
    x, y = ring_arc_mm(THETA1_DEG, THETA2_DEG, R_INNER_MM, 180)
    ax.plot(x, y, color=INK, lw=OBJECT_LW)
    x, y = ring_arc_mm(THETA1_DEG, THETA2_DEG, R_OUTER_MM, 180)
    ax.plot(x, y, color=INK, lw=OBJECT_LW)
    for angle in (THETA1_DEG, THETA2_DEG):
        xi, yi = ring_arc_mm(angle, angle, R_INNER_MM, 1)
        xo, yo = ring_arc_mm(angle, angle, R_OUTER_MM, 1)
        ax.plot([xi[0], xo[0]], [yi[0], yo[0]], color=INK, lw=OBJECT_LW)
    ax.plot([-160, 160], [-302.0, -302.0], color=INK, lw=OBJECT_LW)


def plot_contact_profile(
    dst: Path,
    case_name: str,
    columns: list[str],
    labels: list[str],
    ylabel: str,
    colors: list[str] | None = None,
) -> str:
    rows = read_numeric_csv(RESULTS / case_name / "contact_patch_profiles.csv")
    x = np.array([float(row["relative_angle_deg"]) for row in rows])
    active = np.array([float(row.get("active", 0.0)) > 0.0 for row in rows])
    colors = colors or [RED, BLUE, GREEN]
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    if np.any(active):
        ax.axvspan(float(x[active].min()), float(x[active].max()), color=LIGHT_GRAY, zorder=0)
    for column, label, color in zip(columns, labels, colors):
        y = np.array([float(row[column]) for row in rows])
        ax.plot(x, y, color=color, lw=OBJECT_LW, marker="o", ms=2.8, label=label)
    ax.axhline(0, color=GRAY, lw=THIN_LW)
    ax.set_xlabel("Относительный угол по внешней дуге, град")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    if len(columns) > 1:
        ax.legend(loc="best")
    return save(fig, dst)


def contact_force_boundary(case_name: str) -> Callable[[Path], str]:
    def builder(dst: Path) -> str:
        return plot_contact_profile(
            dst,
            case_name,
            ["integrated_normal_force_kn"],
            ["нормальная контактная сила"],
            "Сила на фасетке, кН",
            [RED],
        )

    return builder


def penetration_boundary(case_name: str) -> Callable[[Path], str]:
    def builder(dst: Path) -> str:
        return plot_contact_profile(
            dst,
            case_name,
            ["average_penetration", "maximum_penetration"],
            ["средняя пенетрация", "максимальная пенетрация"],
            "Пенетрация, мм",
            [BLUE, RED],
        )

    return builder


def contact_pressure_profile(case_name: str) -> Callable[[Path], str]:
    def builder(dst: Path) -> str:
        return plot_contact_profile(
            dst,
            case_name,
            ["average_pressure_mpa"],
            ["контактное давление"],
            "Контактное давление, МПа",
            [GREEN],
        )

    return builder


def figure_2_1(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    setup_plain(ax, (-0.6, 7.2), (-0.4, 5.5))
    pts = np.array([[1.0, 0.8], [2.2, 0.25], [4.7, 0.55], [6.2, 1.7], [5.75, 3.7], [4.1, 4.8], [1.65, 4.25], [0.45, 2.5]])
    ax.add_patch(patches.Polygon(pts, closed=True, facecolor=LIGHT_GRAY, edgecolor=INK, lw=OBJECT_LW))
    ax.plot(pts[[6, 7, 0], 0], pts[[6, 7, 0], 1], color=BLUE, lw=ACCENT_LW)
    ax.plot(pts[[2, 3, 4], 0], pts[[2, 3, 4], 1], color=RED, lw=ACCENT_LW)
    ax.text(3.3, 2.55, "Ω", fontsize=24, ha="center", color=GRAY)
    ax.text(0.75, 3.65, "Γᵤ", fontsize=15, color=BLUE, fontweight="bold")
    ax.text(6.08, 3.0, "Γₜ", fontsize=15, color=RED, fontweight="bold")
    for p in [(0.9, 1.05), (0.75, 1.85), (0.85, 2.65), (1.15, 3.45)]:
        ax.plot([p[0] - 0.25, p[0] + 0.25], [p[1] - 0.25, p[1] + 0.25], color=BLUE, lw=DIM_LW)
        ax.plot([p[0] - 0.25, p[0] + 0.25], [p[1] + 0.25, p[1] - 0.25], color=BLUE, lw=DIM_LW)
    for start, end in [((5.85, 1.6), (6.55, 1.35)), ((6.0, 2.25), (6.75, 2.25)), ((5.75, 3.0), (6.45, 3.25))]:
        arrow(ax, start, end, RED)
    ax.text(5.25, 0.8, "поверхностная\nнагрузка t", color=RED, ha="center", fontsize=9)
    ax.text(1.75, 0.25, "заданные\nперемещения u", color=BLUE, ha="center", fontsize=9)
    return save(fig, dst)


def figure_2_2(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7, 4.8))
    setup_plain(ax, (-0.6, 7.0), (-0.5, 5.1))
    pts = np.array([[1.1, 0.8], [5.6, 0.55], [6.0, 3.8], [0.75, 4.2]])
    ax.add_patch(patches.Polygon(pts, closed=True, facecolor=LIGHT_GRAY, edgecolor=INK, lw=OBJECT_LW))
    for i, (x, y) in enumerate(pts, start=1):
        ax.scatter(x, y, s=50, color=BLUE, zorder=3)
        ax.text(x + 0.12, y + 0.12, f"{i}", color=BLUE, fontweight="bold")
        arrow(ax, (x, y), (x + 0.45, y + 0.32), ORANGE, lw=1.2)
        ax.text(x + 0.48, y + 0.35, f"u{i}, v{i}", color=ORANGE, fontsize=8)
    px, py = 3.45, 2.55
    ax.scatter(px, py, s=42, color=RED, zorder=4)
    ax.text(px + 0.2, py + 0.1, "u(x,y)", color=RED, fontsize=10)
    for x, y in pts:
        ax.plot([x, px], [y, py], color=MID_GRAY, lw=0.8, ls=":")
    ax.text(3.45, 4.55, "u(x,y) = Σ Nᵢ(x,y) uᵢ", ha="center", fontsize=13)
    return save(fig, dst)


def figure_2_3(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    setup_plain(ax, (-1.9, 2.1), (-1.7, 1.9))
    square = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    ax.add_patch(patches.Polygon(square, closed=True, facecolor=LIGHT_GRAY, edgecolor=INK, lw=OBJECT_LW))
    ax.axhline(0, color=GRAY, lw=0.9)
    ax.axvline(0, color=GRAY, lw=0.9)
    arrow(ax, (-1.5, 0), (1.55, 0), GRAY, lw=1)
    arrow(ax, (0, -1.35), (0, 1.45), GRAY, lw=1)
    ax.text(1.65, -0.1, "ξ", fontsize=14)
    ax.text(0.1, 1.5, "η", fontsize=14)
    labels = ["1 (-1,-1)", "2 (1,-1)", "3 (1,1)", "4 (-1,1)"]
    offsets = [(-0.35, -0.3), (0.1, -0.3), (0.1, 0.12), (-0.75, 0.12)]
    for (x, y), label, (dx, dy) in zip(square, labels, offsets):
        ax.scatter(x, y, s=56, color=BLUE, zorder=3)
        ax.text(x + dx, y + dy, label, fontsize=9, color=BLUE)
    gp = 1 / math.sqrt(3)
    for x in [-gp, gp]:
        for y in [-gp, gp]:
            ax.scatter(x, y, marker="x", s=45, color=RED, lw=1.0)
    return save(fig, dst)


def figure_2_4(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    setup_plain(ax, (0, 10), (0, 6))
    matrix_text = (
        "B = [  N₁,x    0     N₂,x    0     N₃,x    0     N₄,x    0  ]\n"
        "    [   0     N₁,y    0     N₂,y    0     N₃,y    0     N₄,y ]\n"
        "    [  N₁,y   N₁,x   N₂,y   N₂,x   N₃,y   N₃,x   N₄,y   N₄,x ]"
    )
    ax.text(0.8, 3.2, matrix_text, family="DejaVu Sans Mono", fontsize=12, va="center")
    ax.text(5.0, 5.1, "ε = B uₑ", ha="center", fontsize=18, color=BLUE, fontweight="bold")
    ax.text(2.1, 1.0, "εₓ", fontsize=11, color=GRAY)
    ax.text(2.1, 0.65, "εᵧ", fontsize=11, color=GRAY)
    ax.text(2.1, 0.3, "γₓᵧ", fontsize=11, color=GRAY)
    ax.text(6.6, 0.55, "uᵢ, vᵢ", ha="center", fontsize=9, color=GRAY)
    return save(fig, dst)


def figure_2_5(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    setup_plain(ax, (0, 11), (0, 6))
    box(ax, (0.4, 3.6), 2.1, 1.1, "геометрия\nэлемента", "#f8fbff")
    box(ax, (0.4, 1.8), 2.1, 1.1, "матрица\nD", "#fffaf2", ORANGE)
    box(ax, (3.4, 2.7), 2.3, 1.2, "B(ξᵢ,ηᵢ),\ndet J", "#f8fbff")
    box(ax, (6.7, 2.7), 3.6, 1.2, "Kₑ = Σ Bᵀ D B detJ wᵢwⱼ", "#f7fff8", GREEN)
    arrow(ax, (2.5, 4.15), (3.4, 3.55), BLUE)
    arrow(ax, (2.5, 2.35), (3.4, 3.0), ORANGE)
    arrow(ax, (5.7, 3.3), (6.7, 3.3), GREEN)
    for x in [4.35, 4.75]:
        for y in [1.0, 1.4]:
            ax.scatter(x, y, marker="x", color=RED, s=50)
    ax.add_patch(patches.Rectangle((4.05, 0.7), 1.0, 1.0, fill=False, edgecolor=BLUE, lw=1.2))
    return save(fig, dst)


def figure_2_6(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    setup_plain(ax, (-1.55, 1.75), (-1.45, 1.65))
    ax.add_patch(patches.Rectangle((-1, -1), 2, 2, facecolor=LIGHT_GRAY, edgecolor=INK, lw=OBJECT_LW))
    ax.axhline(0, color=MID_GRAY, lw=THIN_LW)
    ax.axvline(0, color=MID_GRAY, lw=THIN_LW)
    gp = 1 / math.sqrt(3)
    for x in [-gp, gp]:
        ax.axvline(x, color="#dde3ea", lw=THIN_LW, ls="--")
        ax.text(x, -1.22, "±1/√3" if x > 0 else "-1/√3", ha="center", fontsize=8, color=GRAY)
    for y in [-gp, gp]:
        ax.axhline(y, color="#dde3ea", lw=THIN_LW, ls="--")
    for x in [-gp, gp]:
        for y in [-gp, gp]:
            ax.scatter(x, y, marker="x", s=55, color=RED, lw=1.0)
    arrow(ax, (-1.25, 0), (1.35, 0), GRAY)
    arrow(ax, (0, -1.25), (0, 1.35), GRAY)
    ax.text(1.43, -0.05, "ξ", fontsize=13)
    ax.text(0.05, 1.4, "η", fontsize=13)
    return save(fig, dst)


def figure_2_7(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    setup_plain(ax, (0, 10), (0, 6))
    for x0, label in [(0.7, "Kₑ¹"), (0.7, "Kₑ²")]:
        pass
    def small_matrix(x: float, y: float, label: str, color: str) -> None:
        ax.text(x + 0.55, y + 1.15, label, ha="center", color=color, fontweight="bold")
        for i in range(3):
            for j in range(3):
                ax.add_patch(patches.Rectangle((x + 0.35 * j, y + 0.35 * i), 0.35, 0.35, facecolor="#f8fbff", edgecolor=color, lw=0.8))

    small_matrix(0.8, 3.0, "Kₑ¹", BLUE)
    small_matrix(0.8, 1.2, "Kₑ²", ORANGE)
    ax.text(2.8, 3.45, "локальные\nстепени свободы", ha="center", fontsize=9)
    ax.text(2.8, 1.65, "карта\nсборки", ha="center", fontsize=9)
    arrow(ax, (2.0, 3.45), (3.8, 3.45), BLUE)
    arrow(ax, (2.0, 1.65), (3.8, 2.35), ORANGE)
    ax.text(6.6, 5.05, "глобальная матрица K", ha="center", fontsize=12, color=GRAY)
    for i in range(7):
        for j in range(7):
            face = "white"
            if 1 <= i <= 3 and 1 <= j <= 3:
                face = "#dcebf7"
            if 3 <= i <= 5 and 2 <= j <= 4:
                face = "#fdebd3"
            ax.add_patch(patches.Rectangle((4.4 + 0.42 * j, 1.2 + 0.42 * i), 0.42, 0.42, facecolor=face, edgecolor=MID_GRAY, lw=0.6))
    return save(fig, dst)


def figure_2_8(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    setup_plain(ax, (-3.7, 4.2), (-3.5, 1.5))
    ring_sector(ax, (0, 0), 3.0, 2.45, 205, 335)
    add_ring_mesh(ax, (0, 0), 205, 335)
    for y in np.linspace(-2.6, -1.9, 6):
        ax.plot([-3.1, -2.7], [y - 0.18, y + 0.18], color=BLUE, lw=1)
    ax.text(-3.15, -1.65, "закрепления", color=BLUE, fontsize=9, ha="center")
    for x in np.linspace(-0.9, 0.9, 5):
        arrow(ax, (x, -3.05), (x, -2.7), RED)
    ax.text(0, -3.3, "внешняя нагрузка / перемещение", color=RED, fontsize=9, ha="center")
    box(ax, (1.7, -0.2), 2.0, 0.9, "K u = F", "#ffffff", GRAY)
    box(ax, (1.7, -1.45), 2.0, 0.9, "K_ff u_f = F_f", "#f7fff8", GREEN)
    arrow(ax, (2.7, -0.2), (2.7, -0.55), GREEN)
    return save(fig, dst)


def figure_2_9(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    setup_plain(ax, (-0.5, 8.2), (-0.4, 5.2))
    pts = np.array([[0.8, 0.8], [4.4, 0.6], [4.8, 3.8], [0.6, 4.2]])
    ax.add_patch(patches.Polygon(pts, closed=True, facecolor=LIGHT_GRAY, edgecolor=INK, lw=OBJECT_LW))
    for x, y in pts:
        ax.scatter(x, y, s=42, color=BLUE)
    gp = np.array([[2.0, 1.7], [3.35, 1.6], [3.45, 2.85], [2.05, 3.0]])
    for i, (x, y) in enumerate(gp, 1):
        ax.scatter(x, y, marker="x", s=55, color=RED, lw=1.0)
        ax.text(x + 0.1, y + 0.1, f"σᵍ{i}", color=RED, fontsize=9)
    box(ax, (5.55, 2.5), 2.1, 0.9, "σ в точках\nГаусса", "#fff8f8", RED)
    box(ax, (5.55, 1.0), 2.1, 0.9, "экстраполяция\nи усреднение", "#f8fbff", BLUE)
    arrow(ax, (4.8, 2.5), (5.55, 2.95), RED)
    arrow(ax, (6.6, 2.5), (6.6, 1.9), BLUE)
    return save(fig, dst)


def figure_2_10(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    setup_plain(ax, (-3.2, 5.8), (-2.8, 3.0))
    ring_sector(ax, (0, 0), 2.5, 1.9, 25, 105)
    p = (1.25, 2.15)
    ax.scatter(*p, s=42, color=RED)
    arrow(ax, (0, 0), p, BLUE)
    er = np.array(p) / np.linalg.norm(p)
    et = np.array([-er[1], er[0]])
    arrow(ax, p, tuple(np.array(p) + 0.8 * er), BLUE)
    arrow(ax, p, tuple(np.array(p) + 0.8 * et), ORANGE)
    ax.text(p[0] + 0.75 * er[0] + 0.1, p[1] + 0.75 * er[1], "eᵣ", color=BLUE, fontsize=11)
    ax.text(p[0] + 0.85 * et[0], p[1] + 0.85 * et[1], "eθ", color=ORANGE, fontsize=11)
    arrow(ax, (-2.8, -2.2), (-1.6, -2.2), GRAY)
    arrow(ax, (-2.8, -2.2), (-2.8, -1.0), GRAY)
    ax.text(-1.5, -2.28, "x")
    ax.text(-2.95, -0.9, "y")
    ax.text(2.8, 1.5, "σ_rr = σ_xx cos²φ + σ_yy sin²φ\n       + 2τ_xy sinφ cosφ", fontsize=11)
    ax.text(2.8, 0.25, "σ_θθ, τ_rθ — аналогично\nчерез поворот осей", fontsize=10, color=GRAY)
    ax.text(0.3, 0.45, "φ", fontsize=12, color=GRAY)
    return save(fig, dst)


def figure_2_12(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    setup_plain(ax, (-3.8, 3.8), (-3.35, 0.9))
    ring_sector(ax, (0, 0), 3.0, 2.5, 205, 335, "#f8fbff")
    add_ring_mesh(ax, (0, 0), 205, 335)
    for theta, label, dx in [(260, "край пятна\nконтакта", -1.0), (280, "край пятна\nконтакта", 1.0), (270, "зона\nмаксимального\nсжатия", 0.0)]:
        x, y = outer_arc_points(theta - 3, theta + 3, 3.02, (0, 0), 20)
        ax.plot(x, y, color=RED if theta == 270 else ORANGE, lw=ACCENT_LW, solid_capstyle="butt")
        rad = math.radians(theta)
        lx, ly = 3.32 * math.cos(rad) + dx, 3.32 * math.sin(rad) - 0.2
        ax.text(lx, ly, label, ha="center", fontsize=8, color=RED if theta == 270 else ORANGE)
    ax.plot([-3.4, 3.4], [-3.05, -3.05], color=GRAY, lw=1.2)
    ax.text(2.45, -3.25, "жесткая плоскость", color=GRAY, fontsize=8)
    return save(fig, dst)


def figure_3_3(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    setup_plain(ax, (-3.8, 3.8), (-3.35, 0.9))
    ring_sector(ax, (0, 0), 3.0, 2.5, 205, 335, "#f8fbff")
    add_ring_mesh(ax, (0, 0), 205, 335)
    x, y = outer_arc_points(225, 315, 3.02, (0, 0), 100)
    ax.plot(x, y, color=BLUE, lw=ACCENT_LW, label="потенциальная область")
    x, y = outer_arc_points(260, 280, 3.08, (0, 0), 60)
    ax.plot(x, y, color=RED, lw=ACCENT_LW, solid_capstyle="butt", label="активная область")
    ax.plot([-3.4, 3.4], [-3.05, -3.05], color=INK, lw=OBJECT_LW)
    ax.text(-1.7, -2.45, "Γ_c^pot", color=BLUE, fontsize=12)
    ax.text(0.45, -3.0, "Γ_c^act", color=RED, fontsize=12)
    ax.legend(loc="upper center", ncol=2, frameon=False)
    return save(fig, dst)


def figure_3_1(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    setup_plain(ax, (-4.35, 4.35), (-3.9, 0.95))
    ring_sector(ax, (0, 0), 3.0, 2.5, 240, 300)
    add_ring_mesh(ax, (0, 0), 240, 300)
    plane_y = -3.35
    outer_low_y = -3.0
    ax.plot([-3.85, 3.85], [plane_y, plane_y], color=INK, lw=OBJECT_LW)
    for x0 in np.arange(-3.75, 3.9, 0.35):
        ax.plot([x0 - 0.18, x0 + 0.18], [plane_y - 0.16, plane_y], color="#8a8a8a", lw=THIN_LW)

    for angle in (240, 300):
        ax.plot(
            [0, 3.32 * math.cos(math.radians(angle))],
            [0, 3.32 * math.sin(math.radians(angle))],
            color="#9da3aa",
            lw=DIM_LW,
            ls=(0, (4, 4)),
        )
    ax.scatter([0], [0], s=12, color=INK, zorder=4)
    ax.text(0.08, 0.08, "O", ha="left", va="bottom", fontsize=10, color=INK)

    arrow(ax, (0, 0), (2.5 * math.cos(math.radians(248)), 2.5 * math.sin(math.radians(248))), BLUE)
    arrow(ax, (0, 0), (3.0 * math.cos(math.radians(292)), 3.0 * math.sin(math.radians(292))), BLUE)
    ax.annotate("Rᵢ = 250 мм", xy=(2.5 * math.cos(math.radians(248)), 2.5 * math.sin(math.radians(248))), xytext=(-3.55, -1.15),
                ha="left", va="center", fontsize=9, color=BLUE,
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=DIM_LW, shrinkA=4, shrinkB=4))
    ax.annotate("Rₒ = 300 мм", xy=(3.0 * math.cos(math.radians(300)), 3.0 * math.sin(math.radians(300))), xytext=(2.78, -0.66),
                ha="left", va="center", fontsize=9, color=BLUE,
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=DIM_LW, shrinkA=4, shrinkB=4))

    ax.add_patch(patches.Arc((0, 0), 1.35, 1.35, theta1=240, theta2=300, color=BLUE, lw=DIM_LW))
    ax.text(0, -0.84, "φ = 60°", color=BLUE, fontsize=10, ha="center", va="center")

    x_gap = 1.05
    double_arrow(ax, (x_gap, plane_y), (x_gap, outer_low_y), RED)
    ax.plot([0.18, x_gap + 0.14], [outer_low_y, outer_low_y], color=RED, lw=THIN_LW, ls=(0, (3, 3)))
    ax.plot([0.18, x_gap + 0.14], [plane_y, plane_y], color=RED, lw=THIN_LW, ls=(0, (3, 3)))
    ax.text(1.23, -3.18, "g₀ = 2 мм", color=RED, fontsize=9, va="center", ha="left")
    for x in (-0.42, 0.0, 0.42):
        arrow(ax, (x, -2.31), (x, -2.78), RED)
    ax.annotate("заданное\nперемещение u_y", xy=(0.42, -2.47), xytext=(2.32, -2.52),
                ha="left", va="center", fontsize=9, color=RED,
                arrowprops=dict(arrowstyle="-", color=RED, lw=DIM_LW, shrinkA=4, shrinkB=4))
    return save(fig, dst)


def figure_3_4(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    setup_plain(ax, (0, 10), (0, 5.8))
    X = np.array([[1.0, 1.1], [3.4, 1.0], [3.2, 3.4], [0.8, 3.2]])
    x = np.array([[6.0, 0.9], [8.75, 1.25], [8.15, 4.0], [5.65, 3.25]])
    ax.add_patch(patches.Polygon(X, closed=True, facecolor=LIGHT_GRAY, edgecolor=BLUE, lw=OBJECT_LW))
    ax.add_patch(patches.Polygon(x, closed=True, facecolor="#fffafa", edgecolor=RED, lw=OBJECT_LW))
    ax.text(2.1, 4.0, "начальная\nконфигурация X", ha="center", color=BLUE)
    ax.text(7.15, 4.7, "текущая\nконфигурация x", ha="center", color=RED)
    arrow(ax, (3.7, 2.35), (5.35, 2.35), GRAY, lw=OBJECT_LW)
    ax.text(4.5, 2.65, "x = X + u", ha="center", fontsize=11)
    arrow(ax, (1.8, 1.7), (2.55, 2.55), BLUE)
    arrow(ax, (6.55, 1.7), (7.45, 2.95), RED)
    ax.text(1.65, 1.45, "dX", color=BLUE)
    ax.text(6.35, 1.45, "dx", color=RED)
    ax.text(4.75, 0.6, "F = ∂x / ∂X", ha="center", fontsize=10, color=GRAY)
    return save(fig, dst)


def figure_3_5(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    setup_plain(ax, (-1, 8), (-0.5, 5))
    ax.plot([0.3, 7.5], [0.9, 0.9], color=INK, lw=OBJECT_LW)
    ax.text(6.1, 0.55, "жесткая плоскость", color=GRAY, fontsize=9)
    xs = np.linspace(1.1, 6.8, 120)
    ys = 2.35 + 0.25 * np.sin((xs - 1.1) * np.pi / 5.7)
    ax.plot(xs, ys, color=BLUE, lw=ACCENT_LW)
    xp = 4.15
    yp = float(2.35 + 0.25 * np.sin((xp - 1.1) * np.pi / 5.7))
    ax.scatter(xp, yp, color=RED, s=44)
    double_arrow(ax, (xp, 0.9), (xp, yp), RED)
    ax.text(xp + 0.15, (yp + 0.9) / 2, "gₙ", color=RED, fontsize=13)
    arrow(ax, (xp - 1.0, 0.9), (xp - 1.0, 1.85), GRAY)
    ax.text(xp - 0.8, 1.55, "n", fontsize=12, color=GRAY)
    ax.text(2.2, 3.15, "текущее положение\nточки поверхности", color=BLUE, ha="center")
    return save(fig, dst)


def figure_3_6(dst: Path) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.8))
    states = [("отрыв", "gₙ > 0\npₙ = 0", 2.5, False), ("касание", "gₙ = 0\npₙ = 0", 1.1, False), ("активный контакт", "gₙ = 0\npₙ > 0", 1.1, True)]
    for ax, (title, text, ybody, pressure) in zip(axes, states):
        setup_plain(ax, (0, 4), (0, 4))
        ax.plot([0.5, 3.5], [1, 1], color=INK, lw=OBJECT_LW)
        ax.plot([0.8, 3.2], [ybody, ybody], color=BLUE, lw=ACCENT_LW)
        ax.scatter(2.0, ybody, color=RED, s=35)
        if pressure:
            arrow(ax, (2.0, 0.95), (2.0, 1.75), RED, lw=OBJECT_LW)
            ax.text(2.15, 1.45, "pₙ", color=RED)
        elif ybody > 1.2:
            double_arrow(ax, (2.0, 1.0), (2.0, ybody), RED)
        ax.text(2, 3.35, title, ha="center", fontweight="bold")
        ax.text(2, 0.2, text, ha="center", fontsize=10)
    fig.tight_layout()
    return save(fig, dst)


def figure_3_7(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    setup_plain(ax, (-3.8, 3.8), (-3.35, 0.9))
    ring_sector(ax, (0, 0), 3.0, 2.5, 210, 330, "#f8fbff")
    add_ring_mesh(ax, (0, 0), 210, 330)
    for a1, a2 in zip(np.linspace(230, 310, 9)[:-1], np.linspace(230, 310, 9)[1:]):
        x, y = outer_arc_points(a1, a2, 3.05, (0, 0), 8)
        ax.plot(x, y, color=RED, lw=ACCENT_LW, solid_capstyle="butt")
    mid = math.radians(245)
    label_point = (3.05 * math.cos(mid), 3.05 * math.sin(mid))
    ax.annotate(
        "контактная фасетка",
        xy=label_point,
        xytext=(-2.65, -2.82),
        ha="center",
        fontsize=9,
        arrowprops=dict(arrowstyle="->", color=RED, lw=1.1),
        color=RED,
    )
    ax.text(0, -3.12, "потенциальная контактная граница Γc", ha="center", fontsize=9, color=BLUE)
    return save(fig, dst)


def figure_3_8(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    setup_plain(ax, (0, 8), (0, 4.8))
    ax.plot([0.7, 7.4], [0.8, 0.8], color=INK, lw=OBJECT_LW)
    ax.text(6.0, 0.45, "жесткая плоскость", color=GRAY, fontsize=9)
    p1, p2 = np.array([1.4, 2.35]), np.array([6.6, 2.0])
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=BLUE, lw=ACCENT_LW)
    ax.scatter([p1[0], p2[0]], [p1[1], p2[1]], color=BLUE, s=35)
    for t, active in [(0.28, False), (0.72, True)]:
        p = p1 * (1 - t) + p2 * t
        ax.scatter(p[0], p[1], marker="x", s=55, color=RED, lw=1.0)
        double_arrow(ax, (p[0], 0.8), (p[0], p[1]), RED, lw=1.2)
        ax.text(p[0] + 0.12, (p[1] + 0.8) / 2, "gₙ", color=RED, fontsize=10)
        ax.text(p[0], p[1] + 0.35, "активна" if active else "неактивна", ha="center", fontsize=8, color=GREEN if active else GRAY)
    return save(fig, dst)


def figure_3_9(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    setup_plain(ax, (0, 12), (0, 6))
    box(ax, (0.5, 3.6), 2.0, 0.9, "uᵏ", "#f8fbff")
    box(ax, (3.2, 3.6), 2.3, 0.9, "сборка R(uᵏ), Kᵏ", "#f8fbff")
    box(ax, (6.3, 3.6), 2.2, 0.9, "Kᵏ Δu = -R", "#fffaf2", ORANGE)
    box(ax, (9.2, 3.6), 2.0, 0.9, "uᵏ⁺¹ = uᵏ + αΔu", "#f7fff8", GREEN)
    for x1, x2 in [(2.5, 3.2), (5.5, 6.3), (8.5, 9.2)]:
        arrow(ax, (x1, 4.05), (x2, 4.05), GRAY)
    box(ax, (4.7, 1.3), 2.5, 0.9, "проверка\nсходимости", "#ffffff", GRAY)
    arrow(ax, (10.2, 3.6), (6.0, 2.2), GRAY)
    arrow(ax, (4.7, 1.75), (1.5, 3.6), GRAY)
    ax.text(2.6, 2.2, "если нет", fontsize=9, color=GRAY)
    arrow(ax, (7.2, 1.75), (10.1, 1.75), GREEN)
    ax.text(9.0, 2.0, "если да: следующий\nшаг нагрузки", color=GREEN, fontsize=9)
    return save(fig, dst)


def figure_3_10(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("номер шага")
    ax.set_ylabel("параметр нагрузки λ")
    ax.grid(True, ls="--", lw=0.5, color="#d4dae2")
    x = np.array([0, 1.3, 2.6, 3.9, 5.2, 6.5, 7.8, 9.1])
    y = np.array([0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.88, 1.0])
    ax.plot(x[:5], y[:5], marker="o", color=BLUE, lw=OBJECT_LW, label="успешные шаги")
    ax.plot([5.2, 7.8], [0.60, 0.88], color=RED, lw=OBJECT_LW, ls="--", label="неудачный крупный шаг")
    ax.plot([5.2, 6.5, 7.8, 9.1], [0.60, 0.75, 0.88, 1.0], marker="o", color=GREEN, lw=OBJECT_LW, label="деление шага")
    ax.annotate("шаг делится", xy=(6.5, 0.75), xytext=(6.0, 0.45), arrowprops=dict(arrowstyle="->", color=RED), color=RED)
    ax.legend(frameon=True, loc="lower right")
    return save(fig, dst)


def figure_3_11(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    setup_plain(ax, (0, 12), (0, 6))
    box(ax, (0.6, 3.8), 2.0, 0.9, "итерация\nНьютона", "#f8fbff")
    box(ax, (3.5, 3.8), 2.5, 0.9, "касательная\nсистема", "#f8fbff")
    box(ax, (7.0, 3.8), 2.7, 0.9, "линейный\nрешатель", "#fffaf2", ORANGE)
    box(ax, (7.0, 2.1), 2.7, 0.9, "CG + IC\nили LDLT", "#ffffff", ORANGE)
    box(ax, (10.4, 3.8), 1.3, 0.9, "Δu", "#f7fff8", GREEN)
    arrow(ax, (2.6, 4.25), (3.5, 4.25), GRAY)
    arrow(ax, (6.0, 4.25), (7.0, 4.25), GRAY)
    arrow(ax, (8.35, 3.8), (8.35, 3.0), ORANGE)
    arrow(ax, (9.7, 4.25), (10.4, 4.25), GREEN)
    return save(fig, dst)


def figure_4_1(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    setup_plain(ax, (0, 12), (0, 6))
    box(ax, (0.6, 2.6), 2.3, 1.0, "условие\nнепроникания\n gₙ ≥ 0", "#f8fbff")
    box(ax, (4.1, 3.9), 2.7, 1.0, "штрафной\nметод", "#fff8f8", RED)
    box(ax, (4.1, 1.1), 2.7, 1.0, "расширенный\nЛагранж", "#f7fff8", GREEN)
    box(ax, (8.0, 3.9), 3.0, 1.0, "pₙ = εₙ ⟨-gₙ⟩", "#ffffff", RED)
    box(ax, (8.0, 1.1), 3.0, 1.0, "pₙ = λ + ρ⟨-gₙ⟩", "#ffffff", GREEN)
    arrow(ax, (2.9, 3.1), (4.1, 4.4), RED)
    arrow(ax, (2.9, 3.1), (4.1, 1.6), GREEN)
    arrow(ax, (6.8, 4.4), (8.0, 4.4), RED)
    arrow(ax, (6.8, 1.6), (8.0, 1.6), GREEN)
    return save(fig, dst)


def figure_4_2(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    setup_plain(ax, (0, 12), (0, 6))
    box(ax, (0.7, 3.8), 2.0, 0.9, "λᵐ", "#f7fff8", GREEN)
    box(ax, (3.4, 3.8), 2.7, 0.9, "решение\nравновесия uᵐ", "#f8fbff")
    box(ax, (6.9, 3.8), 2.0, 0.9, "зазоры gₙ", "#fffaf2", ORANGE)
    box(ax, (9.8, 3.8), 1.8, 0.9, "λᵐ⁺¹", "#f7fff8", GREEN)
    for x1, x2 in [(2.7, 3.4), (6.1, 6.9), (8.9, 9.8)]:
        arrow(ax, (x1, 4.25), (x2, 4.25), GRAY)
    ax.text(6.15, 2.35, "λᵐ⁺¹ = max(0, λᵐ + ρ ⟨-gₙ⟩)", ha="center", fontsize=13)
    arrow(ax, (10.7, 3.8), (1.7, 3.8), GREEN, lw=1.0, connectionstyle="arc3,rad=-0.25")
    return save(fig, dst)


def figure_4_3(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    setup_plain(ax, (0, 12), (0, 6))
    ax.add_patch(patches.Rectangle((0.6, 0.7), 10.8, 4.7, fill=False, edgecolor=GREEN, lw=OBJECT_LW))
    ax.text(1.0, 5.15, "внешний цикл AL: m = 0, 1, ...", color=GREEN, fontsize=11, fontweight="bold")
    ax.add_patch(patches.Rectangle((2.0, 1.5), 7.8, 2.7, fill=False, edgecolor=BLUE, lw=OBJECT_LW, ls="--"))
    ax.text(2.35, 3.9, "внутренний цикл Ньютона: k = 0, 1, ...", color=BLUE, fontsize=10)
    box(ax, (2.5, 2.4), 1.7, 0.8, "R,K", "#f8fbff")
    box(ax, (5.0, 2.4), 1.7, 0.8, "Δu", "#fffaf2", ORANGE)
    box(ax, (7.5, 2.4), 1.7, 0.8, "u ← u+Δu", "#f7fff8", GREEN)
    arrow(ax, (4.2, 2.8), (5.0, 2.8), GRAY)
    arrow(ax, (6.7, 2.8), (7.5, 2.8), GRAY)
    box(ax, (4.6, 0.9), 3.1, 0.8, "обновление λ\nпосле внутренней сходимости", "#f7fff8", GREEN)
    arrow(ax, (8.35, 2.4), (6.15, 1.7), GREEN)
    return save(fig, dst)


def figure_4_4(dst: Path) -> str:
    cases, _ = read_triplet_summary()
    labels = ["Штрафной\nметод", "Расширенный\nЛагранж"]
    values = [float(cases["penalty_contact"]["max_penetration"]), float(cases["augmented_lagrangian_contact"]["max_penetration"])]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.bar(labels, values, color=[RED, GREEN], width=0.55)
    ax.set_yscale("log")
    ax.set_ylim(min(values) * 0.35, max(values) * 4.0)
    ax.set_ylabel("Максимальная пенетрация, мм")
    ax.grid(True, axis="y", which="both", ls="--", lw=0.5, color="#d4dae2")
    for i, value in enumerate(values):
        ax.text(i, value * 1.45, f"{value:.3e}", ha="center", fontsize=9)
    return save(fig, dst)


def figure_4_5(dst: Path) -> str:
    cases, _ = read_triplet_summary()
    labels = ["Штрафной\nметод", "Расширенный\nЛагранж"]
    values = [float(cases["penalty_contact"]["total_time_seconds"]), float(cases["augmented_lagrangian_contact"]["total_time_seconds"])]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.bar(labels, values, color=[RED, GREEN], width=0.55)
    ax.set_ylim(0, max(values) * 1.18)
    ax.set_ylabel("Время расчета, с")
    ax.grid(True, axis="y", ls="--", lw=0.5, color="#d4dae2")
    for i, value in enumerate(values):
        ax.text(i, value + max(values) * 0.035, f"{value:.2f}", ha="center", fontsize=9)
    return save(fig, dst)


def figure_4_6(dst: Path) -> str:
    penalty = read_numeric_csv(RESULTS / "penalty_contact" / "contact_patch_profiles.csv")
    al = read_numeric_csv(RESULTS / "augmented_lagrangian_contact" / "contact_patch_profiles.csv")
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for rows, label, color in [(penalty, "штрафной метод", RED), (al, "расширенный Лагранж", GREEN)]:
        x = [float(r["relative_angle_deg"]) for r in rows]
        y = [float(r["average_pressure_mpa"]) for r in rows]
        ax.plot(x, y, lw=OBJECT_LW, color=color, label=label)
    ax.set_xlabel("Относительный угол по внешней дуге, град")
    ax.set_ylabel("Среднее контактное давление, МПа")
    ax.grid(True, ls="--", lw=0.5, color="#d4dae2")
    ax.legend(frameon=True)
    return save(fig, dst)


def get_font(size: int) -> ImageFont.ImageFont:
    for name in ["arial.ttf", "DejaVuSans.ttf", "calibri.ttf"]:
        try:
            return ImageFont.truetype(name, size)
        except OSError:
            continue
    return ImageFont.load_default()


def trim_image(img: Image.Image) -> Image.Image:
    bg = Image.new(img.mode, img.size, "white")
    diff = ImageChops.difference(img.convert("RGB"), bg.convert("RGB"))
    bbox = diff.getbbox()
    if bbox:
        return img.crop(bbox)
    return img


def compose_side_by_side(sources: list[Path], labels: list[str], dst: Path, source_note: str) -> str:
    images = [trim_image(Image.open(src).convert("RGB")) for src in sources]
    target_h = 780
    resized = []
    for img in images:
        scale = target_h / img.height
        resized.append(img.resize((int(img.width * scale), target_h), Image.Resampling.LANCZOS))
    label_h = 80
    gap = 36
    margin = 32
    width = sum(img.width for img in resized) + gap * (len(resized) - 1) + 2 * margin
    height = target_h + label_h + 2 * margin
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    font = get_font(34)
    x = margin
    for img, label in zip(resized, labels):
        draw.text((x + img.width / 2, margin + 10), label, fill=(34, 34, 34), font=font, anchor="ma")
        canvas.paste(img, (x, margin + label_h))
        x += img.width + gap
    dst.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(dst, dpi=(220, 220))
    return source_note


def draw_active_arc_panel(ax: plt.Axes, rows: list[dict[str, float | str]], title: str) -> None:
    draw_results_ring_mesh(ax)
    active_length = 0.0
    for row in rows:
        if float(row.get("active_length", 0.0)) <= 0.0 and float(row.get("active", 0.0)) <= 0.0:
            continue
        x_mid = float(row["reference_midpoint_x"])
        y_mid = float(row["reference_midpoint_y"])
        radius = float(np.hypot(x_mid, y_mid))
        theta = math.degrees(math.atan2(y_mid, x_mid))
        length = max(float(row.get("active_length", 0.0)), 0.35 * float(row.get("facet_length", 0.0)))
        half_angle = math.degrees(length / max(radius, 1.0)) / 2.0
        x, y = ring_arc_mm(theta - half_angle, theta + half_angle, radius + 2.2, 10)
        ax.plot(x, y, color=RED, lw=ACCENT_LW, solid_capstyle="butt")
        active_length += float(row.get("active_length", 0.0))
    if title:
        ax.text(0, -214, title, ha="center", va="center", fontsize=10, color=GRAY, fontweight="bold")
    ax.text(0, -309.5, f"активная длина: {active_length:.2f} мм", ha="center", fontsize=8.5, color=RED)


def figure_4_7(dst: Path) -> str:
    penalty = read_numeric_csv(RESULTS / "penalty_contact" / "contact_facets.csv")
    al = read_numeric_csv(RESULTS / "augmented_lagrangian_contact" / "contact_facets.csv")
    fig, axes = plt.subplots(2, 1, figsize=(7.8, 6.2))
    draw_active_arc_panel(axes[0], penalty, "Штрафной метод")
    draw_active_arc_panel(axes[1], al, "Расширенный Лагранж")
    fig.subplots_adjust(hspace=0.08)
    return save(fig, dst)


def figure_5_12(dst: Path) -> str:
    al = read_numeric_csv(RESULTS / "augmented_lagrangian_contact" / "contact_facets.csv")
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    draw_active_arc_panel(ax, al, "")
    return save(fig, dst)


def figure_5_1(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(9, 4.8))
    setup_plain(ax, (-4.2, 7.2), (-3.3, 1.3))
    ring_sector(ax, (-2.0, 0.0), 2.1, 1.65, 215, 325, "#f8fbff")
    ax.plot([-4.1, 0.1], [-2.25, -2.25], color=INK, lw=OBJECT_LW)
    for x in np.linspace(-2.6, -1.4, 5):
        arrow(ax, (x, -2.3), (x, -1.85), RED)
    ax.text(-2.0, 0.75, "явный контакт:\nдавление найдено\nиз условий gₙ ≥ 0", ha="center", fontsize=9)
    arrow(ax, (0.5, -1.0), (2.0, -1.0), GRAY, lw=OBJECT_LW)
    ax.text(1.25, -0.65, "замена", ha="center", fontsize=10)
    ring_sector(ax, (4.1, 0.0), 2.1, 1.65, 215, 325, "#f8fbff")
    for a in np.linspace(260, 280, 7):
        rad = math.radians(a)
        start = (4.1 + 2.25 * math.cos(rad), 2.25 * math.sin(rad))
        end = (4.1 + 1.85 * math.cos(rad), 1.85 * math.sin(rad))
        arrow(ax, start, end, RED)
    ax.text(4.1, 0.75, "бесконтактная задача:\nзаданная нагрузка p(s)", ha="center", fontsize=9)
    return save(fig, dst)


def figure_5_14(dst: Path) -> str:
    _, params = read_triplet_summary()
    half_width = params["surrogate_contact_half_width"]
    p0 = params["surrogate_max_pressure"]
    al = read_numeric_csv(RESULTS / "augmented_lagrangian_contact" / "contact_patch_profiles.csv")
    s = np.array([float(r["arc_coordinate"]) for r in al])
    p = np.array([float(r["average_pressure_mpa"]) for r in al])
    s0 = s[np.argmax(p)]
    x = s - s0
    xs = np.linspace(-half_width * 1.08, half_width * 1.08, 300)
    parabola = np.maximum(0.0, p0 * (1.0 - (xs / half_width) ** 2))
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(x, p, "o", color=GREEN, ms=4, label="давление AL")
    ax.plot(xs, parabola, color=RED, lw=OBJECT_LW, label="параболическая аппроксимация")
    ax.axvline(-half_width, color=GRAY, ls="--", lw=0.9)
    ax.axvline(half_width, color=GRAY, ls="--", lw=0.9)
    ax.text(0, p0 * 1.04, "p₀", ha="center", color=RED)
    ax.annotate("a", xy=(half_width, 0), xytext=(half_width / 2, p0 * 0.25), arrowprops=dict(arrowstyle="->", color=GRAY), color=GRAY)
    ax.set_xlabel("Координата относительно центра пятна, мм")
    ax.set_ylabel("Контактное давление, МПа")
    ax.grid(True, ls="--", lw=0.5, color="#d4dae2")
    ax.legend(frameon=True)
    return save(fig, dst)


def figure_5_15(dst: Path) -> str:
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    setup_plain(ax, (-3.8, 3.8), (-3.75, 0.9))
    ring_sector(ax, (0, 0), 3.0, 2.5, 205, 335, "#f8fbff")
    add_ring_mesh(ax, (0, 0), 205, 335)
    angles = np.linspace(260, 280, 11)
    for a in angles:
        weight = 1.0 - ((a - 270) / 10) ** 2
        rad = math.radians(a)
        start = (3.25 * math.cos(rad), 3.25 * math.sin(rad))
        end = ((3.25 - 0.25 - 0.35 * weight) * math.cos(rad), (3.25 - 0.25 - 0.35 * weight) * math.sin(rad))
        arrow(ax, start, end, RED, lw=DIM_LW + 0.45 * weight)
    ax.text(0, -3.52, "p(s) = p₀(1 - (s/a)²), |s| ≤ a", ha="center", color=RED, fontsize=10)
    return save(fig, dst)


def composite_fields(number: str, field: str, labels: list[str], cases: list[str]) -> Callable[[Path], str]:
    def builder(dst: Path) -> str:
        sources = [RESULTS / case / field for case in cases]
        note = "composite: " + " + ".join(str(src.relative_to(ROOT)) for src in sources)
        return compose_side_by_side(sources, labels, dst, note)

    return builder


def figure_5_17(dst: Path) -> str:
    al = read_numeric_csv(RESULTS / "augmented_lagrangian_contact" / "contact_patch_profiles.csv")
    _, params = read_triplet_summary()
    half_width = params["surrogate_contact_half_width"]
    p0 = params["surrogate_max_pressure"]
    s = np.array([float(r["arc_coordinate"]) for r in al])
    p_eq = np.maximum(0.0, p0 * (1.0 - (s / half_width) ** 2))
    equivalent_force_kn = p_eq * np.array([float(r["facet_length"]) for r in al]) * 1.0e-3
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(
        [float(r["relative_angle_deg"]) for r in al],
        [float(r["integrated_normal_force_kn"]) for r in al],
        color=GREEN,
        lw=OBJECT_LW,
        label="контактная сила AL",
    )
    ax.plot(
        [float(r["relative_angle_deg"]) for r in al],
        equivalent_force_kn,
        color=RED,
        lw=OBJECT_LW,
        ls="--",
        label="эквивалентная нагрузка",
    )
    ax.set_xlabel("Относительный угол по внешней дуге, град")
    ax.set_ylabel("Нормальная сила на участке, кН")
    ax.grid(True, ls="--", lw=0.5, color="#d4dae2")
    ax.legend(frameon=True)
    return save(fig, dst)


def figure_5_20(dst: Path) -> str:
    al = read_numeric_csv(RESULTS / "augmented_lagrangian_contact" / "ring_contour_stress_profiles.csv")
    sur = read_numeric_csv(RESULTS / "no_contact_surrogate" / "ring_contour_stress_profiles.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4), sharex=True)
    for ax, component, title in zip(axes, ["sigma_rr_mpa", "sigma_tt_mpa"], [r"$\sigma_{rr}$", r"$\sigma_{\theta\theta}$"]):
        for rows, label, color, ls in [(al, "явный контакт AL", GREEN, "-"), (sur, "бесконтактная нагрузка", RED, "--")]:
            outer = [r for r in rows if str(r["contour"]) == "outer"]
            ax.plot([float(r["relative_angle_deg"]) for r in outer], [float(r[component]) for r in outer], color=color, ls=ls, lw=OBJECT_LW, label=label)
        ax.set_title(title)
        ax.set_xlabel("Угол, град")
        ax.grid(True, ls="--", lw=0.5, color="#d4dae2")
    axes[0].set_ylabel("Напряжение, МПа")
    axes[1].legend(frameon=True, fontsize=8)
    return save(fig, dst)


def figure_5_21(dst: Path) -> str:
    al = read_numeric_csv(RESULTS / "augmented_lagrangian_contact" / "ring_radial_section_profiles.csv")
    sur = read_numeric_csv(RESULTS / "no_contact_surrogate" / "ring_radial_section_profiles.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4), sharex=True)
    radius_key = next(k for k in al[0].keys() if str(k).startswith("radius_"))
    for ax, component, title in zip(axes, ["sigma_rr_mpa", "sigma_tt_mpa"], [r"$\sigma_{rr}$", r"$\sigma_{\theta\theta}$"]):
        for rows, label, color, ls in [(al, "явный контакт AL", GREEN, "-"), (sur, "бесконтактная нагрузка", RED, "--")]:
            ax.plot([float(r[radius_key]) for r in rows], [float(r[component]) for r in rows], color=color, ls=ls, lw=OBJECT_LW, label=label)
        ax.set_title(title)
        ax.set_xlabel("Радиус, мм")
        ax.grid(True, ls="--", lw=0.5, color="#d4dae2")
    axes[0].set_ylabel("Напряжение, МПа")
    axes[1].legend(frameon=True, fontsize=8)
    return save(fig, dst)


def build_tasks() -> list[FigureTask]:
    penalty = RESULTS / "penalty_contact"
    al = RESULTS / "augmented_lagrangian_contact"
    no_contact = RESULTS / "no_contact_surrogate"
    return [
        FigureTask("2.1", "Схема области тела с разбиением границы на участки Γᵤ и Γₜ", figure_2_1),
        FigureTask("2.2", "Аппроксимация поля перемещений внутри конечного элемента по узловым значениям", figure_2_2),
        FigureTask("2.3", "Четырехузловой изопараметрический конечный элемент в локальных координатах ξ, η", figure_2_3),
        FigureTask("2.4", "Структура матрицы деформаций четырехузлового элемента", figure_2_4),
        FigureTask("2.5", "Формирование элементной матрицы жесткости четырехузлового элемента", figure_2_5),
        FigureTask("2.6", "Схема квадратуры Гаусса 2×2 для четырехузлового элемента", figure_2_6),
        FigureTask("2.7", "Схема сборки локальных матриц в глобальную систему", figure_2_7),
        FigureTask("2.8", "Схема учета закреплений и внешних нагрузок в конечно-элементной модели", figure_2_8),
        FigureTask("2.9", "Точки Гаусса четырехузлового элемента и восстановление напряжений", figure_2_9),
        FigureTask("2.10", "Преобразование компонентов напряжений из декартовой системы координат в полярную", figure_2_10),
        FigureTask("2.11", "Типовая конечно-элементная сетка поперечного сечения массивной шины с локальным сгущением в зоне контакта", copy_image(penalty / "computational_mesh.png")),
        FigureTask("2.12", "Характерные зоны концентрации напряжений в поперечном сечении шины", figure_2_12),
        FigureTask("3.1", "Расчетная схема кольцевого сектора шины и жесткой опорной плоскости", figure_3_1),
        FigureTask("3.2", "Конечно-элементная сетка с локальным сгущением в зоне контакта", copy_image(penalty / "computational_mesh.png")),
        FigureTask("3.3", "Потенциальная и активная контактные области на внешнем контуре сектора", figure_3_3),
        FigureTask("3.4", "Переход от начальной конфигурации к текущей и смысл градиента деформации", figure_3_4),
        FigureTask("3.5", "Геометрический смысл нормального зазора", figure_3_5),
        FigureTask("3.6", "Возможные состояния точки контактной поверхности: отрыв, касание и активный контакт", figure_3_6),
        FigureTask("3.7", "Дискретизация внешнего контура сектора контактными фасетками", figure_3_7),
        FigureTask("3.8", "Точки Гаусса на контактной фасетке и проверка нормального зазора", figure_3_8),
        FigureTask("3.9", "Схема ньютоновского итерационного процесса для контактной задачи", figure_3_9),
        FigureTask("3.10", "Поэтапное нагружение и адаптивное деление неудачного шага", figure_3_10),
        FigureTask("3.11", "Место линейного решателя внутри итерации Ньютона", figure_3_11),
        FigureTask("4.1", "Схема сравнения штрафного метода и метода расширенного Лагранжа", figure_4_1),
        FigureTask("4.2", "Внешние итерации метода расширенного Лагранжа и обновление контактных множителей", figure_4_2),
        FigureTask("4.3", "Внешне-внутренний итерационный цикл метода расширенного Лагранжа", figure_4_3),
        FigureTask("4.4", "Сравнение максимальной пенетрации для штрафного метода и метода расширенного Лагранжа", figure_4_4),
        FigureTask("4.5", "Сравнение времени расчета для штрафного метода и метода расширенного Лагранжа", figure_4_5),
        FigureTask("4.6", "Распределение контактного давления по дуге для штрафного метода и метода расширенного Лагранжа", figure_4_6),
        FigureTask("4.7", "Активная контактная зона для штрафного метода и метода расширенного Лагранжа", figure_4_7),
        FigureTask("5.1", "Переход от явной контактной задачи к эквивалентной бесконтактной нагрузке", figure_5_1),
        FigureTask("5.2", "Поле перемещений в расчете явного контакта штрафным методом", copy_image(penalty / "displacement_magnitude.png")),
        FigureTask("5.3", "Модуль контактной силы в расчете штрафным методом", contact_force_boundary("penalty_contact")),
        FigureTask("5.4", "Распределение эквивалентных напряжений по Мизесу в расчете штрафным методом", copy_image(penalty / "von_mises_stress.png")),
        FigureTask("5.5", "Распределение определителя градиента деформации в расчете штрафным методом", copy_image(penalty / "jacobian_determinant.png")),
        FigureTask("5.6", "Поле пенетрации в расчете штрафным методом", penetration_boundary("penalty_contact")),
        FigureTask("5.7", "Профиль контактного давления по дуге в расчете штрафным методом", contact_pressure_profile("penalty_contact")),
        FigureTask("5.8", "Поле перемещений в расчете методом расширенного Лагранжа", copy_image(al / "displacement_magnitude.png")),
        FigureTask("5.9", "Модуль контактной силы в расчете методом расширенного Лагранжа", copy_image(al / "contact_force_magnitude.png")),
        FigureTask("5.10", "Распределение эквивалентных напряжений по Мизесу в расчете методом расширенного Лагранжа", copy_image(al / "von_mises_stress.png")),
        FigureTask("5.11", "Распределение определителя градиента деформации в расчете методом расширенного Лагранжа", copy_image(al / "jacobian_determinant.png")),
        FigureTask("5.12", "Активная контактная область в расчете методом расширенного Лагранжа", figure_5_12),
        FigureTask("5.13", "Профиль контактного давления по дуге в расчете методом расширенного Лагранжа", contact_pressure_profile("augmented_lagrangian_contact")),
        FigureTask("5.14", "Параболическая аппроксимация контактного давления на внешнем контуре шины", figure_5_14),
        FigureTask("5.15", "Приложение восстановленного параболического давления в бесконтактной задаче", figure_5_15),
        FigureTask("5.16", "Сравнение полей перемещений для явного контакта и бесконтактной аппроксимации", composite_fields("5.16", "displacement_magnitude.png", ["Явный контакт AL", "Бесконтактная нагрузка"], ["augmented_lagrangian_contact", "no_contact_surrogate"])),
        FigureTask("5.17", "Сравнение модулей контактных и эквивалентных нагрузочных сил", figure_5_17),
        FigureTask("5.18", "Сравнение эквивалентных напряжений по Мизесу", composite_fields("5.18", "von_mises_stress.png", ["Явный контакт AL", "Бесконтактная нагрузка"], ["augmented_lagrangian_contact", "no_contact_surrogate"])),
        FigureTask("5.19", "Сравнение плотности энергии деформации", composite_fields("5.19", "strain_energy_density.png", ["Явный контакт AL", "Бесконтактная нагрузка"], ["augmented_lagrangian_contact", "no_contact_surrogate"])),
        FigureTask("5.20", "Профили напряжений по внешнему контуру шины", figure_5_20),
        FigureTask("5.21", "Профили напряжений по радиальному сечению", figure_5_21),
    ]


def find_docx() -> Path | None:
    downloads = Path.home() / "Downloads"
    candidates = [p for p in downloads.glob("*FEM*нумерац*.docx") if "(1)" not in p.name]
    if not candidates:
        candidates = list(downloads.glob("*FEM*.docx"))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def extract_doc_context() -> list[tuple[str, list[str]]]:
    path = find_docx()
    if path is None:
        return []
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    paragraphs: list[str] = []
    for paragraph in root.findall(".//w:body/w:p", ns):
        text = "".join((node.text or "") for node in paragraph.findall(".//w:t", ns)).strip()
        if text:
            paragraphs.append(text)
    contexts: list[tuple[str, list[str]]] = []
    for index, text in enumerate(paragraphs):
        if text.startswith("Рисунок "):
            contexts.append((text, paragraphs[max(0, index - 3): min(len(paragraphs), index + 4)]))
    return contexts


def write_manifest(tasks: list[FigureTask]) -> None:
    with (OUT / "manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["number", "file", "vector_file", "caption", "source"])
        writer.writeheader()
        writer.writerows(MANIFEST)

    with (OUT / "README.md").open("w", encoding="utf-8-sig") as handle:
        handle.write("# Пронумерованные рисунки для диплома\n\n")
        handle.write("Документ Word не изменялся. Расчетные рисунки взяты из `results/main_scale_hyperelastic_reference_triplet_coarse_symmetric_anchor`, где соблюдено условие симметрии.\n\n")
        handle.write("Сгенерированные SVG/PNG-схемы приведены к более строгому инженерному стилю: тонкие основные и выносные линии, прямоугольные блоки, умеренные цветовые акценты и штриховка опорных поверхностей по логике ЕСКД.\n\n")
        handle.write("| Рисунок | PNG | SVG | Источник |\n")
        handle.write("| --- | --- | --- | --- |\n")
        for item in MANIFEST:
            vector = f"`{item['vector_file']}`" if item["vector_file"] else "—"
            handle.write(f"| {item['number']} | `{item['file']}` | {vector} | {item['source']} |\n")

    contexts = extract_doc_context()
    if contexts:
        with (OUT / "figure_context.md").open("w", encoding="utf-8-sig") as handle:
            handle.write("# Контекст подписей из docx\n\n")
            for caption, context in contexts:
                handle.write(f"## {caption}\n\n")
                for line in context:
                    handle.write(f"- {line}\n")
                handle.write("\n")


def main() -> None:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = build_tasks()
    for task in tasks:
        dst = OUT / file_name(task.number)
        svg_dst = dst.with_suffix(".svg")
        if dst.exists():
            dst.unlink()
        if svg_dst.exists():
            svg_dst.unlink()
        source = task.builder(dst)
        vector_file = svg_dst.name if svg_dst.exists() else ""
        MANIFEST.append(
            {
                "number": task.number,
                "file": dst.name,
                "vector_file": vector_file,
                "caption": task.caption,
                "source": source,
            }
        )
    write_manifest(tasks)
    print(f"Generated {len(tasks)} figures in {OUT}")


if __name__ == "__main__":
    main()
