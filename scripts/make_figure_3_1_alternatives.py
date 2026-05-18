#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures_numbered" / "alternatives"

INK = "#202020"
GRAY = "#5c5c5c"
LIGHT_GRAY = "#f4f5f6"
MESH_GRAY = "#bfc5cc"
BLUE = "#2f5f8f"
RED = "#8f2f2f"


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif", "Liberation Serif"],
            "font.size": 10,
            "axes.linewidth": 0.8,
            "mathtext.fontset": "dejavuserif",
            "savefig.dpi": 300,
            "svg.fonttype": "none",
        }
    )


def setup(ax: plt.Axes) -> None:
    ax.set_xlim(-4.4, 4.4)
    ax.set_ylim(-3.9, 0.95)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def polar(radius: float, angle_deg: float) -> tuple[float, float]:
    angle = math.radians(angle_deg)
    return radius * math.cos(angle), radius * math.sin(angle)


def arc_points(radius: float, start: float, end: float, n: int = 100) -> tuple[np.ndarray, np.ndarray]:
    angles = np.deg2rad(np.linspace(start, end, n))
    return radius * np.cos(angles), radius * np.sin(angles)


def add_sector(ax: plt.Axes) -> None:
    theta1, theta2 = 240, 300
    r_inner, r_outer = 2.5, 3.0

    sector = patches.Wedge(
        (0, 0),
        r_outer,
        theta1,
        theta2,
        width=r_outer - r_inner,
        facecolor=LIGHT_GRAY,
        edgecolor=INK,
        linewidth=1.15,
        joinstyle="miter",
    )
    ax.add_patch(sector)

    for radius in np.linspace(r_inner, r_outer, 5):
        x, y = arc_points(radius, theta1, theta2)
        ax.plot(x, y, color=MESH_GRAY, lw=0.45)

    for angle in np.linspace(theta1, theta2, 9):
        xi, yi = polar(r_inner, angle)
        xo, yo = polar(r_outer, angle)
        ax.plot([xi, xo], [yi, yo], color=MESH_GRAY, lw=0.45)

    for angle in (theta1, theta2):
        xo, yo = polar(3.32, angle)
        ax.plot([0, xo], [0, yo], color="#9da3aa", lw=0.75, ls=(0, (4, 4)))

    ax.scatter([0], [0], s=12, color=INK, zorder=4)
    ax.text(0.08, 0.08, "O", ha="left", va="bottom", fontsize=9, color=INK)


def add_dimension_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    lw: float = 0.9,
    scale: float = 9,
) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "<->",
            "color": color,
            "lw": lw,
            "shrinkA": 0,
            "shrinkB": 0,
            "mutation_scale": scale,
        },
    )


def add_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    lw: float = 0.9,
    scale: float = 9,
) -> None:
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "->",
            "color": color,
            "lw": lw,
            "shrinkA": 0,
            "shrinkB": 0,
            "mutation_scale": scale,
        },
    )


def add_callout(
    ax: plt.Axes,
    text: str,
    target: tuple[float, float],
    text_xy: tuple[float, float],
    color: str = INK,
    ha: str = "center",
) -> None:
    ax.annotate(
        text,
        xy=target,
        xytext=text_xy,
        ha=ha,
        va="center",
        fontsize=9.2,
        color=color,
        linespacing=1.15,
        arrowprops={
            "arrowstyle": "-",
            "color": color,
            "lw": 0.75,
            "shrinkA": 4,
            "shrinkB": 4,
            "connectionstyle": "arc3,rad=0.0",
        },
    )


def add_plane(ax: plt.Axes) -> None:
    plane_y = -3.35
    ax.plot([-3.85, 3.85], [plane_y, plane_y], color=INK, lw=1.05)
    for x0 in np.arange(-3.75, 3.9, 0.35):
        ax.plot([x0 - 0.18, x0 + 0.18], [plane_y - 0.16, plane_y], color="#8a8a8a", lw=0.55)

    ax.text(2.65, -3.68, "жесткая опорная плоскость", ha="center", va="center", fontsize=9.2, color=GRAY)


def add_annotations(ax: plt.Axes) -> None:
    # Opening angle.
    ax.add_patch(patches.Arc((0, 0), 1.35, 1.35, theta1=240, theta2=300, color=BLUE, lw=0.95))
    ax.text(0, -0.84, r"$\varphi = 60^\circ$", ha="center", va="center", fontsize=10.3, color=BLUE)

    # Radius dimension lines are kept inside the empty sector center, labels are outside.
    add_arrow(ax, (0, 0), polar(2.5, 248), BLUE, lw=0.8, scale=8)
    add_arrow(ax, (0, 0), polar(3.0, 292), BLUE, lw=0.8, scale=8)
    add_callout(ax, r"$R_i = 250$ мм", polar(2.5, 248), (-3.45, -1.12), BLUE, ha="left")
    add_callout(ax, r"$R_o = 300$ мм", polar(3.0, 292), (2.72, -1.12), BLUE, ha="left")

    # Initial gap.
    plane_y = -3.35
    outer_bottom = -3.0
    x_gap = 1.05
    add_dimension_arrow(ax, (x_gap, plane_y), (x_gap, outer_bottom), RED, lw=0.85, scale=8)
    ax.plot([0.18, x_gap + 0.14], [outer_bottom, outer_bottom], color=RED, lw=0.65, ls=(0, (3, 3)))
    ax.plot([0.18, x_gap + 0.14], [plane_y, plane_y], color=RED, lw=0.65, ls=(0, (3, 3)))
    ax.text(1.23, -3.18, r"$g_0 = 2$ мм", ha="left", va="center", fontsize=9.2, color=RED)

    # Prescribed displacement.
    for x in (-0.42, 0.0, 0.42):
        add_arrow(ax, (x, -2.31), (x, -2.78), RED, lw=0.85, scale=8)
    add_callout(ax, "заданное\nперемещение $u_y$", (0.42, -2.47), (2.42, -2.24), RED, ha="left")

    ax.text(-3.2, 0.34, "кольцевой сектор шины", ha="left", va="center", fontsize=9.4, color=GRAY)

    # A compact reference axis in the margin helps read the direction convention.
    ax.plot([-4.05, -3.45], [-3.18, -3.18], color=GRAY, lw=0.75)
    ax.plot([-4.05, -4.05], [-3.18, -2.58], color=GRAY, lw=0.75)
    add_arrow(ax, (-3.45, -3.18), (-3.3, -3.18), GRAY, lw=0.75, scale=7)
    add_arrow(ax, (-4.05, -2.58), (-4.05, -2.43), GRAY, lw=0.75, scale=7)
    ax.text(-3.25, -3.18, "x", ha="left", va="center", fontsize=8.5, color=GRAY)
    ax.text(-4.05, -2.35, "y", ha="center", va="bottom", fontsize=8.5, color=GRAY)


def build_vector_refined() -> tuple[Path, Path]:
    configure_style()
    OUT.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    setup(ax)
    add_sector(ax)
    add_plane(ax)
    add_annotations(ax)

    png = OUT / "fig_3_01_vector_refined.png"
    svg = OUT / "fig_3_01_vector_refined.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, svg


def main() -> None:
    png, svg = build_vector_refined()
    print(png)
    print(svg)


if __name__ == "__main__":
    main()
