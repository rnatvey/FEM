#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

try:
    import pyvista as pv
except ImportError as exc:  # pragma: no cover - import guard for runtime environment
    raise SystemExit(
        "PyVista is required for this postprocessor. Install pyvista and vtk in the Python environment."
    ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read FEM VTU + metrics.json outputs, run simple checks, and save plots."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        default=["results/ring_contact_study"],
        help="Case directories or roots that contain exported case directories.",
    )
    parser.add_argument(
        "--warp-factor",
        type=float,
        default=1.0,
        help="Visual warp factor applied to displacement plots.",
    )
    parser.add_argument(
        "--summary-name",
        default="summary_metrics.png",
        help="Filename for the summary metrics plot saved in the common root.",
    )
    return parser.parse_args()


def discover_case_directories(paths: list[str]) -> list[Path]:
    case_dirs: list[Path] = []
    for raw_path in paths:
        root = Path(raw_path)
        if not root.exists():
            continue

        if (root / "metrics.json").exists() and list(root.glob("*.vtu")):
            case_dirs.append(root)
            continue

        for metrics_path in root.rglob("metrics.json"):
            case_dir = metrics_path.parent
            if list(case_dir.glob("*.vtu")):
                case_dirs.append(case_dir)

    unique_case_dirs = sorted(set(case_dirs))
    return unique_case_dirs


def load_contact_facet_rows(case_dir: Path, metrics: dict[str, Any]) -> list[dict[str, Any]]:
    contact_file_name = metric_path(metrics, "files", "contact_facets")
    contact_path = case_dir / contact_file_name if contact_file_name else case_dir / "contact_facets.csv"
    if not contact_path.exists():
        return []

    float_fields = {
        "reference_midpoint_x",
        "reference_midpoint_y",
        "deformed_midpoint_x",
        "deformed_midpoint_y",
        "normal_x",
        "normal_y",
        "thickness",
        "facet_length",
        "active_length",
        "integrated_area",
        "active_area",
        "average_gap",
        "average_penetration",
        "maximum_penetration",
        "integrated_normal_force",
        "average_pressure",
    }
    int_fields = {"facet_id", "element_id", "surface_index", "active", "active_gauss_points"}

    rows: list[dict[str, Any]] = []
    with contact_path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        for raw_row in reader:
            row: dict[str, Any] = {}
            for key, value in raw_row.items():
                if value is None:
                    row[key] = value
                elif key in float_fields:
                    row[key] = float(value)
                elif key in int_fields:
                    row[key] = int(value)
                else:
                    row[key] = value
            row["active"] = bool(row.get("active", 0))
            rows.append(row)

    return rows


def load_case(case_dir: Path) -> dict[str, Any]:
    metrics_path = case_dir / "metrics.json"
    with metrics_path.open("r", encoding="utf-8") as stream:
        metrics = json.load(stream)

    vtu_files = sorted(case_dir.glob("*.vtu"))
    if not vtu_files:
        raise FileNotFoundError(f"No VTU file found in {case_dir}")

    mesh = pv.read(vtu_files[0])
    return {
        "case_dir": case_dir,
        "metrics": metrics,
        "mesh": mesh,
        "vtu_path": vtu_files[0],
        "contact_facet_rows": load_contact_facet_rows(case_dir, metrics),
    }


def metric_path(metrics: dict[str, Any], *keys: str, default: Any = None) -> Any:
    current: Any = metrics
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def wrap_angle(angle: float) -> float:
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


def wrap_angles(angles: np.ndarray) -> np.ndarray:
    return (angles + math.pi) % (2.0 * math.pi) - math.pi


def load_ring_metadata(case: dict[str, Any]) -> dict[str, Any] | None:
    metrics = case["metrics"]
    extra = metrics.get("extra", {})
    if extra.get("geometry_family") != "tire_ring":
        return None

    required_numeric_keys = [
        "ring_center_x",
        "ring_center_y",
        "ring_inner_radius",
        "ring_outer_radius",
        "ring_contact_center_angle_rad",
    ]
    if any(key not in extra for key in required_numeric_keys):
        return None

    return {
        "center": np.array([float(extra["ring_center_x"]), float(extra["ring_center_y"])]),
        "inner_radius": float(extra["ring_inner_radius"]),
        "outer_radius": float(extra["ring_outer_radius"]),
        "start_angle": float(extra.get("ring_start_angle_rad", math.nan)),
        "end_angle": float(extra.get("ring_end_angle_rad", math.nan)),
        "contact_center_angle": float(extra["ring_contact_center_angle_rad"]),
        "contact_half_angle": float(extra.get("ring_contact_half_angle_rad", math.nan)),
        "radial_step_hint": float(extra.get("mesh_min_radial_step", 0.0)),
        "angular_step_hint": float(extra.get("mesh_min_angular_step", 0.0)),
        "plane_normal": np.array(
            [
                float(extra.get("rigid_plane_normal_x", 0.0)),
                float(extra.get("rigid_plane_normal_y", 1.0)),
            ]
        ),
    }


def cluster_levels(values: np.ndarray, tolerance: float) -> list[float]:
    if values.size == 0:
        return []

    sorted_values = np.sort(values.astype(float))
    clusters: list[list[float]] = [[float(sorted_values[0])]]
    for value in sorted_values[1:]:
        if abs(float(value) - clusters[-1][-1]) <= tolerance:
            clusters[-1].append(float(value))
        else:
            clusters.append([float(value)])
    return [float(sum(cluster) / len(cluster)) for cluster in clusters]


def closest_numeric_level(levels: list[float], target: float) -> float:
    return min(levels, key=lambda value: abs(value - target))


def closest_angle_level(levels: list[float], target: float) -> float:
    return min(levels, key=lambda value: abs(wrap_angle(value - target)))


def point_array(mesh: pv.DataSet, name: str, default_components: int = 1) -> np.ndarray:
    if name not in mesh.point_data:
        if default_components == 1:
            return np.zeros(mesh.n_points)
        return np.zeros((mesh.n_points, default_components))

    array = np.asarray(mesh.point_data[name])
    if default_components == 1:
        return array.reshape(mesh.n_points)
    if array.ndim == 1:
        return array.reshape(mesh.n_points, default_components)
    return array


def write_rows_csv(path: Path, headers: list[str], rows: list[list[Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        writer.writerows(rows)


def build_ring_point_fields(case: dict[str, Any], ring_metadata: dict[str, Any]) -> dict[str, np.ndarray]:
    mesh = case["mesh"]
    points = np.asarray(mesh.points)[:, :2]
    relative_points = points - ring_metadata["center"]
    radii = np.linalg.norm(relative_points, axis=1)
    angles = np.arctan2(relative_points[:, 1], relative_points[:, 0])
    relative_angles = wrap_angles(angles - ring_metadata["contact_center_angle"])
    relative_angles_deg = np.rad2deg(relative_angles)

    displacement = point_array(mesh, "displacement", default_components=3)[:, :2]
    stress = point_array(mesh, "stress_2d", default_components=3)
    contact_force = point_array(mesh, "contact_force", default_components=3)[:, :2]
    signed_distance = point_array(mesh, "rigid_plane_signed_distance")
    penetration = point_array(mesh, "rigid_plane_penetration")

    cosine = np.cos(angles)
    sine = np.sin(angles)

    sigma_xx = stress[:, 0]
    sigma_yy = stress[:, 1]
    tau_xy = stress[:, 2]
    sigma_rr = cosine * cosine * sigma_xx + sine * sine * sigma_yy + 2.0 * sine * cosine * tau_xy
    sigma_tt = sine * sine * sigma_xx + cosine * cosine * sigma_yy - 2.0 * sine * cosine * tau_xy
    tau_rt = -sine * cosine * sigma_xx + sine * cosine * sigma_yy + (cosine * cosine - sine * sine) * tau_xy

    radial_displacement = cosine * displacement[:, 0] + sine * displacement[:, 1]
    tangential_displacement = -sine * displacement[:, 0] + cosine * displacement[:, 1]

    plane_normal = ring_metadata["plane_normal"]
    plane_normal_norm = np.linalg.norm(plane_normal)
    if plane_normal_norm > 0.0:
        plane_normal = plane_normal / plane_normal_norm
    contact_normal_force = contact_force @ plane_normal

    return {
        "points": points,
        "radii": radii,
        "angles": angles,
        "relative_angles": relative_angles,
        "relative_angles_deg": relative_angles_deg,
        "sigma_rr": sigma_rr,
        "sigma_tt": sigma_tt,
        "tau_rt": tau_rt,
        "radial_displacement": radial_displacement,
        "tangential_displacement": tangential_displacement,
        "contact_normal_force": contact_normal_force,
        "signed_distance": signed_distance,
        "penetration": penetration,
    }


def build_ring_contact_facet_fields(
    case: dict[str, Any],
    ring_metadata: dict[str, Any],
) -> list[dict[str, Any]]:
    rows = case.get("contact_facet_rows", [])
    if not rows:
        return []

    enriched_rows: list[dict[str, Any]] = []
    for row in rows:
        midpoint = np.array(
            [float(row["reference_midpoint_x"]), float(row["reference_midpoint_y"])],
            dtype=float,
        )
        relative_point = midpoint - ring_metadata["center"]
        angle = math.atan2(float(relative_point[1]), float(relative_point[0]))
        relative_angle = wrap_angle(angle - ring_metadata["contact_center_angle"])
        enriched_row = dict(row)
        enriched_row["reference_midpoint"] = midpoint
        enriched_row["arc_coordinate"] = ring_metadata["outer_radius"] * relative_angle
        enriched_row["relative_angle_rad"] = relative_angle
        enriched_row["relative_angle_deg"] = math.degrees(relative_angle)
        enriched_rows.append(enriched_row)

    enriched_rows.sort(key=lambda row: float(row["arc_coordinate"]))
    return enriched_rows


def summarize_contact_facet_rows(
    rows: list[dict[str, Any]],
    ring_metadata: dict[str, Any] | None = None,
) -> dict[str, float]:
    if not rows:
        return {}

    total_normal_force = float(sum(float(row["integrated_normal_force"]) for row in rows))
    active_contact_area = float(sum(float(row["active_area"]) for row in rows))
    contact_patch_length = float(sum(float(row["active_length"]) for row in rows))
    max_average_pressure = float(max(float(row["average_pressure"]) for row in rows))
    mean_active_pressure = (
        total_normal_force / active_contact_area if active_contact_area > 0.0 else 0.0
    )
    active_facet_count = float(sum(1.0 for row in rows if bool(row["active"])))

    summary = {
        "total_normal_force": total_normal_force,
        "active_contact_area": active_contact_area,
        "contact_patch_length": contact_patch_length,
        "max_average_pressure": max_average_pressure,
        "mean_active_pressure": mean_active_pressure,
        "active_facet_count": active_facet_count,
    }

    if total_normal_force > 0.0:
        center_of_pressure_x = sum(
            float(row["integrated_normal_force"]) * float(row["reference_midpoint_x"])
            for row in rows
        ) / total_normal_force
        center_of_pressure_y = sum(
            float(row["integrated_normal_force"]) * float(row["reference_midpoint_y"])
            for row in rows
        ) / total_normal_force
        summary["center_of_pressure_x"] = center_of_pressure_x
        summary["center_of_pressure_y"] = center_of_pressure_y

        if ring_metadata is not None:
            cp_relative = np.array(
                [center_of_pressure_x, center_of_pressure_y], dtype=float
            ) - ring_metadata["center"]
            cp_angle = math.atan2(float(cp_relative[1]), float(cp_relative[0]))
            cp_relative_angle = wrap_angle(cp_angle - ring_metadata["contact_center_angle"])
            summary["center_of_pressure_arc"] = ring_metadata["outer_radius"] * cp_relative_angle

    return summary


def save_ring_contour_profiles(case: dict[str, Any], ring_metadata: dict[str, Any]) -> None:
    case_dir = case["case_dir"]
    fields = build_ring_point_fields(case, ring_metadata)

    radial_tolerance = max(
        1.0e-8,
        0.35 * max(
            ring_metadata["radial_step_hint"],
            0.01 * (ring_metadata["outer_radius"] - ring_metadata["inner_radius"]),
        ),
    )
    radial_levels = cluster_levels(fields["radii"], radial_tolerance * 0.5)
    if len(radial_levels) < 3:
        return

    contour_targets = {
        "inner": ring_metadata["inner_radius"],
        "mid": 0.5 * (ring_metadata["inner_radius"] + ring_metadata["outer_radius"]),
        "outer": ring_metadata["outer_radius"],
    }
    contour_styles = {
        "inner": ("#4C78A8", "Inner contour"),
        "mid": ("#F58518", "Mid-surface"),
        "outer": ("#54A24B", "Outer contour"),
    }

    csv_rows: list[list[Any]] = []
    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    axis_specs = [
        ("sigma_rr", "Radial Stress", axes[0]),
        ("sigma_tt", "Circumferential Stress", axes[1]),
        ("tau_rt", "Shear Stress", axes[2]),
    ]

    for contour_key, target_radius in contour_targets.items():
        selected_radius = closest_numeric_level(radial_levels, target_radius)
        indices = np.where(np.abs(fields["radii"] - selected_radius) <= radial_tolerance)[0]
        if indices.size < 2:
            continue

        order = np.argsort(fields["relative_angles"][indices])
        ordered_indices = indices[order]
        x_values = fields["relative_angles_deg"][ordered_indices]

        for field_name, _, axis in axis_specs:
            axis.plot(
                x_values,
                fields[field_name][ordered_indices],
                color=contour_styles[contour_key][0],
                linewidth=2.0,
                label=contour_styles[contour_key][1],
            )

        for ordered_index in ordered_indices:
            csv_rows.append(
                [
                    contour_key,
                    fields["relative_angles_deg"][ordered_index],
                    fields["radii"][ordered_index],
                    fields["sigma_rr"][ordered_index],
                    fields["sigma_tt"][ordered_index],
                    fields["tau_rt"][ordered_index],
                ]
            )

    axes[0].set_title("Stress Profiles Along Tire Contours")
    axes[0].set_ylabel(r"$\sigma_{rr}$")
    axes[1].set_ylabel(r"$\sigma_{\theta\theta}$")
    axes[2].set_ylabel(r"$\tau_{r\theta}$")
    axes[2].set_xlabel("Angle relative to contact center [deg]")
    for axis in axes:
        axis.grid(True, alpha=0.3)
        axis.legend()

    fig.tight_layout()
    fig.savefig(case_dir / "ring_contour_stress_profiles.png", dpi=180)
    plt.close(fig)

    write_rows_csv(
        case_dir / "ring_contour_stress_profiles.csv",
        ["contour", "relative_angle_deg", "radius", "sigma_rr", "sigma_tt", "tau_rt"],
        csv_rows,
    )


def save_ring_radial_profiles(case: dict[str, Any], ring_metadata: dict[str, Any]) -> None:
    case_dir = case["case_dir"]
    fields = build_ring_point_fields(case, ring_metadata)

    angle_tolerance = max(
        1.0e-8,
        0.45 * max(ring_metadata["angular_step_hint"], math.radians(0.5)),
    )
    angle_levels = cluster_levels(fields["relative_angles"], angle_tolerance * 0.5)
    if not angle_levels:
        return

    selected_angle = closest_angle_level(angle_levels, 0.0)
    indices = np.where(np.abs(wrap_angles(fields["relative_angles"] - selected_angle)) <= angle_tolerance)[0]
    if indices.size < 2:
        return

    order = np.argsort(fields["radii"][indices])
    ordered_indices = indices[order]
    radii = fields["radii"][ordered_indices]

    csv_rows = [
        [
            radii[i],
            fields["relative_angles_deg"][ordered_indices[i]],
            fields["sigma_rr"][ordered_indices[i]],
            fields["sigma_tt"][ordered_indices[i]],
            fields["radial_displacement"][ordered_indices[i]],
        ]
        for i in range(len(ordered_indices))
    ]

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    axes[0].plot(radii, fields["sigma_rr"][ordered_indices], color="#4C78A8", linewidth=2.0)
    axes[1].plot(radii, fields["sigma_tt"][ordered_indices], color="#F58518", linewidth=2.0)
    axes[2].plot(radii, fields["radial_displacement"][ordered_indices], color="#54A24B", linewidth=2.0)

    axes[0].set_title("Radial Section Profiles Through Contact Symmetry Plane")
    axes[0].set_ylabel(r"$\sigma_{rr}$")
    axes[1].set_ylabel(r"$\sigma_{\theta\theta}$")
    axes[2].set_ylabel(r"$u_r$")
    axes[2].set_xlabel("Radius")

    for axis in axes:
        axis.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(case_dir / "ring_radial_section_profiles.png", dpi=180)
    plt.close(fig)

    write_rows_csv(
        case_dir / "ring_radial_section_profiles.csv",
        ["radius", "relative_angle_deg", "sigma_rr", "sigma_tt", "radial_displacement"],
        csv_rows,
    )


def save_contact_patch_profiles(case: dict[str, Any], ring_metadata: dict[str, Any]) -> None:
    case_dir = case["case_dir"]
    contact_facet_rows = build_ring_contact_facet_fields(case, ring_metadata)
    if contact_facet_rows:
        arc_coordinate = np.array([float(row["arc_coordinate"]) for row in contact_facet_rows], dtype=float)
        relative_angle_deg = np.array(
            [float(row["relative_angle_deg"]) for row in contact_facet_rows],
            dtype=float,
        )
        average_penetration = np.array(
            [float(row["average_penetration"]) for row in contact_facet_rows],
            dtype=float,
        )
        maximum_penetration = np.array(
            [float(row["maximum_penetration"]) for row in contact_facet_rows],
            dtype=float,
        )
        integrated_normal_force = np.array(
            [float(row["integrated_normal_force"]) for row in contact_facet_rows],
            dtype=float,
        )
        average_pressure = np.array(
            [float(row["average_pressure"]) for row in contact_facet_rows],
            dtype=float,
        )
        active_length = np.array(
            [float(row["active_length"]) for row in contact_facet_rows],
            dtype=float,
        )
        facet_length = np.array(
            [float(row["facet_length"]) for row in contact_facet_rows],
            dtype=float,
        )
        active_mask = np.array([bool(row["active"]) for row in contact_facet_rows], dtype=bool)

        csv_rows = [
            [
                arc_coordinate[i],
                relative_angle_deg[i],
                int(active_mask[i]),
                facet_length[i],
                active_length[i],
                average_penetration[i],
                maximum_penetration[i],
                integrated_normal_force[i],
                average_pressure[i],
            ]
            for i in range(len(contact_facet_rows))
        ]

        fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
        axes[0].plot(arc_coordinate, average_penetration, color="#E45756", linewidth=2.0)
        axes[0].plot(
            arc_coordinate,
            maximum_penetration,
            color="#B279A2",
            linewidth=1.5,
            linestyle="--",
        )
        axes[1].plot(arc_coordinate, integrated_normal_force, color="#4C78A8", linewidth=2.0)
        axes[2].plot(arc_coordinate, average_pressure, color="#72B7B2", linewidth=2.0)
        axes[2].fill_between(
            arc_coordinate,
            0.0,
            average_pressure,
            where=active_mask,
            alpha=0.25,
            color="#72B7B2",
        )

        axes[0].set_title("Facet-Level Contact Patch Profiles")
        axes[0].set_ylabel("Penetration")
        axes[1].set_ylabel("Integrated normal force")
        axes[2].set_ylabel("Average pressure")
        axes[2].set_xlabel("Arc coordinate relative to contact center")

        for axis in axes:
            axis.grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(case_dir / "contact_patch_profiles.png", dpi=180)
        plt.close(fig)

        write_rows_csv(
            case_dir / "contact_patch_profiles.csv",
            [
                "arc_coordinate",
                "relative_angle_deg",
                "active",
                "facet_length",
                "active_length",
                "average_penetration",
                "maximum_penetration",
                "integrated_normal_force",
                "average_pressure",
            ],
            csv_rows,
        )
        return

    fields = build_ring_point_fields(case, ring_metadata)

    radial_tolerance = max(
        1.0e-8,
        0.35 * max(
            ring_metadata["radial_step_hint"],
            0.01 * (ring_metadata["outer_radius"] - ring_metadata["inner_radius"]),
        ),
    )
    radial_levels = cluster_levels(fields["radii"], radial_tolerance * 0.5)
    if not radial_levels:
        return

    outer_radius = closest_numeric_level(radial_levels, ring_metadata["outer_radius"])
    indices = np.where(np.abs(fields["radii"] - outer_radius) <= radial_tolerance)[0]
    if indices.size < 2:
        return

    order = np.argsort(fields["relative_angles"][indices])
    ordered_indices = indices[order]
    ordered_points = fields["points"][ordered_indices]
    ordered_relative_angles = fields["relative_angles"][ordered_indices]
    arc_coordinate = ring_metadata["outer_radius"] * ordered_relative_angles

    segment_lengths = np.linalg.norm(np.diff(ordered_points, axis=0), axis=1)
    tributary_lengths = np.zeros(len(ordered_indices))
    if len(ordered_indices) > 1:
        tributary_lengths[0] = 0.5 * segment_lengths[0]
        tributary_lengths[-1] = 0.5 * segment_lengths[-1]
    if len(ordered_indices) > 2:
        tributary_lengths[1:-1] = 0.5 * (segment_lengths[:-1] + segment_lengths[1:])

    contact_force = fields["contact_normal_force"][ordered_indices]
    traction_estimate = np.divide(
        contact_force,
        np.maximum(tributary_lengths, 1.0e-12),
    )

    csv_rows = [
        [
            arc_coordinate[i],
            math.degrees(ordered_relative_angles[i]),
            fields["penetration"][ordered_indices[i]],
            fields["signed_distance"][ordered_indices[i]],
            contact_force[i],
            traction_estimate[i],
        ]
        for i in range(len(ordered_indices))
    ]

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    axes[0].plot(arc_coordinate, fields["penetration"][ordered_indices], color="#E45756", linewidth=2.0)
    axes[1].plot(arc_coordinate, contact_force, color="#4C78A8", linewidth=2.0)
    axes[2].plot(arc_coordinate, traction_estimate, color="#72B7B2", linewidth=2.0)

    axes[0].set_title("Contact Patch Profiles Along Outer Contour")
    axes[0].set_ylabel("Penetration")
    axes[1].set_ylabel("Normal contact force")
    axes[2].set_ylabel("Normal traction estimate")
    axes[2].set_xlabel("Arc coordinate relative to contact center")

    for axis in axes:
        axis.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(case_dir / "contact_patch_profiles.png", dpi=180)
    plt.close(fig)

    write_rows_csv(
        case_dir / "contact_patch_profiles.csv",
        [
            "arc_coordinate",
            "relative_angle_deg",
            "penetration",
            "signed_distance",
            "contact_normal_force",
            "normal_traction_estimate",
        ],
        csv_rows,
    )


def run_checks(case: dict[str, Any]) -> list[str]:
    mesh = case["mesh"]
    metrics = case["metrics"]
    contact_facet_rows = case.get("contact_facet_rows", [])
    checks: list[str] = []

    expected_nodes = metric_path(metrics, "counts", "nodes")
    expected_elements = metric_path(metrics, "counts", "elements")

    if expected_nodes is not None and int(expected_nodes) != mesh.n_points:
        raise ValueError(
            f"{case['case_dir']}: metrics.json nodes={expected_nodes}, VTU points={mesh.n_points}"
        )
    checks.append(f"nodes={mesh.n_points}")

    if expected_elements is not None and int(expected_elements) != mesh.n_cells:
        raise ValueError(
            f"{case['case_dir']}: metrics.json elements={expected_elements}, VTU cells={mesh.n_cells}"
        )
    checks.append(f"elements={mesh.n_cells}")

    if "displacement_magnitude" in mesh.point_data:
        max_displacement = float(mesh.point_data["displacement_magnitude"].max())
        metrics_max_displacement = metric_path(metrics, "extrema", "max_displacement_magnitude")
        if metrics_max_displacement is not None:
            delta = abs(max_displacement - float(metrics_max_displacement))
            if delta > 1.0e-8 * max(1.0, abs(max_displacement), abs(float(metrics_max_displacement))):
                raise ValueError(
                    f"{case['case_dir']}: displacement magnitude mismatch between VTU and metrics.json"
                )
        checks.append(f"max_displacement={max_displacement:.6e}")

    if "rigid_plane_penetration" in mesh.point_data:
        max_penetration = float(mesh.point_data["rigid_plane_penetration"].max())
        metrics_max_penetration = metric_path(metrics, "contact", "max_nodal_penetration")
        if metrics_max_penetration is not None:
            delta = abs(max_penetration - float(metrics_max_penetration))
            if delta > 1.0e-8 * max(1.0, abs(max_penetration), abs(float(metrics_max_penetration))):
                raise ValueError(
                    f"{case['case_dir']}: penetration mismatch between VTU and metrics.json"
                )
        checks.append(f"max_penetration={max_penetration:.6e}")

    if contact_facet_rows:
        contact_summary = summarize_contact_facet_rows(contact_facet_rows)
        total_normal_force = float(contact_summary["total_normal_force"])
        metrics_total_normal_force = metric_path(metrics, "contact", "total_normal_force")
        if metrics_total_normal_force is not None:
            delta = abs(total_normal_force - float(metrics_total_normal_force))
            if delta > 1.0e-8 * max(1.0, abs(total_normal_force), abs(float(metrics_total_normal_force))):
                raise ValueError(
                    f"{case['case_dir']}: total normal force mismatch between contact_facets.csv and metrics.json"
                )

        contact_patch_length = float(contact_summary["contact_patch_length"])
        metrics_patch_length = metric_path(metrics, "contact", "contact_patch_length")
        if metrics_patch_length is not None:
            delta = abs(contact_patch_length - float(metrics_patch_length))
            if delta > 1.0e-8 * max(1.0, abs(contact_patch_length), abs(float(metrics_patch_length))):
                raise ValueError(
                    f"{case['case_dir']}: contact patch length mismatch between contact_facets.csv and metrics.json"
                )

        checks.append(f"facet_total_force={total_normal_force:.6e}")
        checks.append(f"contact_patch_length={contact_patch_length:.6e}")

    return checks


def save_case_metric_overview(case: dict[str, Any]) -> None:
    case_dir = case["case_dir"]
    metrics = case["metrics"]

    timing_labels = ["assembly", "solve", "total"]
    timing_values = [
        float(metric_path(metrics, "timings", "assembly_time_seconds", default=0.0)),
        float(metric_path(metrics, "timings", "solve_time_seconds", default=0.0)),
        float(metric_path(metrics, "timings", "total_time_seconds", default=0.0)),
    ]

    summary_lines = [
        f"nodes: {metric_path(metrics, 'counts', 'nodes', default='n/a')}",
        f"elements: {metric_path(metrics, 'counts', 'elements', default='n/a')}",
        f"total_dofs: {metric_path(metrics, 'counts', 'total_dofs', default='n/a')}",
        f"free_dofs: {metric_path(metrics, 'counts', 'free_dofs', default='n/a')}",
        f"matrix_nnz: {metric_path(metrics, 'matrix', 'nnz', default='n/a')}",
        f"linear_iterations: {metric_path(metrics, 'iterations', 'linear_iterations', default='n/a')}",
        f"nonlinear_iterations: {metric_path(metrics, 'iterations', 'nonlinear_iterations', default='n/a')}",
        f"active_contact_facets: {metric_path(metrics, 'contact', 'active_contact_facets', default='n/a')}",
        f"max_penetration: {metric_path(metrics, 'contact', 'max_penetration', default='n/a')}",
        f"contact_patch_length: {metric_path(metrics, 'contact', 'contact_patch_length', default='n/a')}",
        f"max_avg_pressure: {metric_path(metrics, 'contact', 'max_facet_average_pressure', default='n/a')}",
        f"total_normal_force: {metric_path(metrics, 'contact', 'total_normal_force', default='n/a')}",
        f"max_displacement: {metric_path(metrics, 'extrema', 'max_displacement_magnitude', default='n/a')}",
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    bars = axes[0].bar(timing_labels, timing_values, color=["#4C78A8", "#F58518", "#54A24B"])
    axes[0].set_ylabel("Time [s]")
    axes[0].set_title("Timing Overview")
    axes[0].grid(True, axis="y", alpha=0.3)

    for bar, value in zip(bars, timing_values):
        axes[0].text(
            bar.get_x() + bar.get_width() / 2.0,
            value,
            f"{value:.3e}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    axes[1].axis("off")
    axes[1].set_title("Case Metrics")
    axes[1].text(
        0.0,
        1.0,
        "\n".join(summary_lines),
        va="top",
        ha="left",
        family="monospace",
        fontsize=10,
    )

    fig.tight_layout()
    fig.savefig(case_dir / "case_overview.png", dpi=180)
    plt.close(fig)


def save_case_plots(case: dict[str, Any], warp_factor: float) -> None:
    case_dir = case["case_dir"]
    mesh = case["mesh"]

    plot_specs = [
        ("displacement_magnitude", "displacement_magnitude.png", "Displacement Magnitude", "viridis"),
        ("sigma_yy", "sigma_yy.png", "Sigma YY", "coolwarm"),
        ("von_mises_stress", "von_mises_stress.png", "Von Mises Stress", "plasma"),
        ("reaction_force_magnitude", "reaction_force_magnitude.png", "Reaction Force Magnitude", "magma"),
        ("active_contact_facet", "active_contact_facet.png", "Active Contact Facets", "binary"),
        (
            "candidate_contact_facet",
            "candidate_contact_facet.png",
            "Candidate Contact Facets",
            "binary",
        ),
    ]
    if "rigid_plane_penetration" in mesh.point_data:
        plot_specs.append(
            ("rigid_plane_penetration", "penetration.png", "Rigid-Plane Penetration", "inferno")
        )
    if "contact_force_magnitude" in mesh.point_data:
        plot_specs.append(
            ("contact_force_magnitude", "contact_force_magnitude.png", "Contact Force Magnitude", "magma")
        )
    if "rigid_plane_signed_distance" in mesh.point_data:
        plot_specs.append(
            ("rigid_plane_signed_distance", "signed_distance.png", "Rigid-Plane Signed Distance", "coolwarm")
        )

    deformed_mesh = mesh
    if "displacement" in mesh.point_data:
        deformed_mesh = mesh.warp_by_vector("displacement", factor=warp_factor)

    available_scalars = set(deformed_mesh.array_names)

    for scalar_name, file_name, title, cmap in plot_specs:
        if scalar_name not in available_scalars:
            continue

        plotter = pv.Plotter(off_screen=True, window_size=(1400, 900))
        plotter.add_text(title, font_size=12)
        plotter.add_mesh(
            deformed_mesh,
            scalars=scalar_name,
            show_edges=True,
            cmap=cmap,
            scalar_bar_args={"title": title},
        )
        plotter.view_xy()
        plotter.camera.zoom(1.2)
        plotter.screenshot(case_dir / file_name)
        plotter.close()


def save_ring_specific_plots(case: dict[str, Any]) -> None:
    ring_metadata = load_ring_metadata(case)
    if ring_metadata is None:
        return

    save_ring_contour_profiles(case, ring_metadata)
    save_ring_radial_profiles(case, ring_metadata)
    save_contact_patch_profiles(case, ring_metadata)


def build_summary_plot(cases: list[dict[str, Any]], output_root: Path, output_name: str) -> bool:
    if not cases:
        return False

    rows: list[dict[str, Any]] = []
    for case in cases:
        metrics = case["metrics"]
        extra = metrics.get("extra", {})
        mesh_label = extra.get("mesh_label", case["case_dir"].name)
        penalty = float(extra.get("penalty_parameter", math.nan))
        max_penetration = float(metric_path(metrics, "contact", "max_penetration", default=math.nan))
        total_time = float(metric_path(metrics, "timings", "total_time_seconds", default=math.nan))
        rows.append(
            {
                "mesh_label": mesh_label,
                "penalty": penalty,
                "max_penetration": max_penetration,
                "total_time": total_time,
            }
        )

    valid_rows = [
        row
        for row in rows
        if math.isfinite(row["penalty"]) and row["penalty"] > 0.0
    ]
    if not valid_rows:
        return False

    grouped_labels = sorted({row["mesh_label"] for row in valid_rows})

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for mesh_label in grouped_labels:
        group = sorted(
            [row for row in valid_rows if row["mesh_label"] == mesh_label],
            key=lambda row: row["penalty"],
        )
        penalties = [row["penalty"] for row in group]
        penetrations = [row["max_penetration"] for row in group]
        total_times = [row["total_time"] for row in group]

        axes[0].plot(penalties, penetrations, marker="o", label=mesh_label)
        axes[1].plot(penalties, total_times, marker="o", label=mesh_label)

    axes[0].set_xscale("log")
    axes[0].set_xlabel("Penalty Parameter")
    axes[0].set_ylabel("Max Penetration")
    axes[0].set_title("Penalty vs Max Penetration")
    axes[0].grid(True, which="both", alpha=0.3)

    axes[1].set_xscale("log")
    axes[1].set_xlabel("Penalty Parameter")
    axes[1].set_ylabel("Total Solve Time [s]")
    axes[1].set_title("Penalty vs Total Time")
    axes[1].grid(True, which="both", alpha=0.3)

    for axis in axes:
        axis.legend()

    fig.tight_layout()
    fig.savefig(output_root / output_name, dpi=180)
    plt.close(fig)
    return True


def build_contact_summary_plot(cases: list[dict[str, Any]], output_root: Path) -> bool:
    if not cases:
        return False

    rows: list[dict[str, Any]] = []
    for case in cases:
        metrics = case["metrics"]
        extra = metrics.get("extra", {})
        mesh_label = extra.get("mesh_label", case["case_dir"].name)
        penalty = float(extra.get("penalty_parameter", math.nan))
        ring_metadata = load_ring_metadata(case)
        facet_summary = summarize_contact_facet_rows(
            case.get("contact_facet_rows", []),
            ring_metadata,
        )
        total_normal_force = float(
            facet_summary.get(
                "total_normal_force",
                metric_path(metrics, "contact", "total_normal_force", default=math.nan),
            )
        )
        contact_patch_length = float(
            facet_summary.get(
                "contact_patch_length",
                metric_path(metrics, "contact", "contact_patch_length", default=math.nan),
            )
        )
        max_average_pressure = float(
            facet_summary.get(
                "max_average_pressure",
                metric_path(metrics, "contact", "max_facet_average_pressure", default=math.nan),
            )
        )
        mean_active_pressure = float(
            facet_summary.get(
                "mean_active_pressure",
                metric_path(metrics, "contact", "mean_active_pressure", default=math.nan),
            )
        )

        row = {
            "case_name": case["case_dir"].name,
            "mesh_label": mesh_label,
            "penalty": penalty,
            "total_normal_force": total_normal_force,
            "contact_patch_length": contact_patch_length,
            "max_average_pressure": max_average_pressure,
            "mean_active_pressure": mean_active_pressure,
        }

        if "center_of_pressure_arc" in facet_summary:
            row["center_of_pressure_arc"] = float(facet_summary["center_of_pressure_arc"])
        else:
            cp_x = metric_path(metrics, "contact", "center_of_pressure_x", default=math.nan)
            cp_y = metric_path(metrics, "contact", "center_of_pressure_y", default=math.nan)
            if ring_metadata is not None and math.isfinite(float(cp_x)) and math.isfinite(float(cp_y)):
                cp_relative = np.array([float(cp_x), float(cp_y)], dtype=float) - ring_metadata["center"]
                cp_angle = math.atan2(float(cp_relative[1]), float(cp_relative[0]))
                cp_relative_angle = wrap_angle(cp_angle - ring_metadata["contact_center_angle"])
                row["center_of_pressure_arc"] = ring_metadata["outer_radius"] * cp_relative_angle
            else:
                row["center_of_pressure_arc"] = math.nan

        rows.append(row)

    valid_rows = [
        row
        for row in rows
        if math.isfinite(row["penalty"]) and row["penalty"] > 0.0
    ]
    if not valid_rows:
        return False

    grouped_labels = sorted({row["mesh_label"] for row in valid_rows})

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    axis_specs = [
        ("contact_patch_length", "Contact patch length", axes[0, 0]),
        ("max_average_pressure", "Max facet average pressure", axes[0, 1]),
        ("total_normal_force", "Total normal force", axes[1, 0]),
        ("mean_active_pressure", "Mean active pressure", axes[1, 1]),
    ]

    for mesh_label in grouped_labels:
        group = sorted(
            [row for row in valid_rows if row["mesh_label"] == mesh_label],
            key=lambda row: row["penalty"],
        )
        penalties = [row["penalty"] for row in group]
        for field_name, title, axis in axis_specs:
            values = [row[field_name] for row in group]
            axis.plot(penalties, values, marker="o", linewidth=2.0, label=mesh_label)
            axis.set_xscale("log")
            axis.set_xlabel("Penalty Parameter")
            axis.set_ylabel(title)
            axis.set_title(f"Penalty vs {title}")
            axis.grid(True, which="both", alpha=0.3)

    for axis in axes.flatten():
        axis.legend()

    fig.tight_layout()
    fig.savefig(output_root / "summary_contact_metrics.png", dpi=180)
    plt.close(fig)

    write_rows_csv(
        output_root / "summary_contact_metrics.csv",
        [
            "case_name",
            "mesh_label",
            "penalty",
            "contact_patch_length",
            "max_average_pressure",
            "mean_active_pressure",
            "total_normal_force",
            "center_of_pressure_arc",
        ],
        [
            [
                row["case_name"],
                row["mesh_label"],
                row["penalty"],
                row["contact_patch_length"],
                row["max_average_pressure"],
                row["mean_active_pressure"],
                row["total_normal_force"],
                row["center_of_pressure_arc"],
            ]
            for row in valid_rows
        ],
    )
    return True


def main() -> int:
    args = parse_args()
    pv.set_plot_theme("document")

    case_dirs = discover_case_directories(args.paths)
    if not case_dirs:
        raise SystemExit("No exported case directories were found.")

    loaded_cases = [load_case(case_dir) for case_dir in case_dirs]
    common_root = Path(args.paths[0])

    for case in loaded_cases:
        checks = run_checks(case)
        save_case_plots(case, args.warp_factor)
        save_case_metric_overview(case)
        save_ring_specific_plots(case)
        print(f"{case['case_dir']}: " + ", ".join(checks))

    if build_summary_plot(loaded_cases, common_root, args.summary_name):
        print(f"Saved summary plot to {common_root / args.summary_name}")
    if build_contact_summary_plot(loaded_cases, common_root):
        print(f"Saved contact summary plots to {common_root / 'summary_contact_metrics.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
