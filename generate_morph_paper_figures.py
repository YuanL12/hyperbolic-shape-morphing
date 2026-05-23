#!/usr/bin/env python3
"""Generate paper figures showing selected morph frames as SVG grids."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from generate_genus3_klein_embedding import geodesic_points


FRAME_INDICES = [0, 3, 6, 9, 12, 15, 18, 21, 24]

MORPHS = [
    (
        "genus2_triangulation",
        "boundary_morph_genus2_bolza_triangulation_mean_value_all1_to_half10_frames",
        "Genus 2 triangulation",
    ),
    (
        "genus2_cells",
        "boundary_morph_genus2_bolza_cells_mean_value_all1_to_half10_frames",
        "Genus 2 cells",
    ),
    (
        "genus3_triangulation",
        "boundary_morph_genus3_triangulation_mean_value_all1_to_half10_frames",
        "Genus 3 triangulation",
    ),
    (
        "genus3_cells",
        "boundary_morph_genus3_cells_mean_value_all1_to_half10_frames",
        "Genus 3 cells",
    ),
]


def load_frame(frame_dir: Path, frame_index: int) -> Dict[str, object]:
    path = frame_dir / f"frame_{frame_index:03d}.json"
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def project(point: Sequence[float], x0: float, y0: float, scale: float) -> Tuple[float, float]:
    return x0 + point[0] * scale, y0 - point[1] * scale


def path_for_geodesic(
    vertices: Sequence[Sequence[float]],
    i: int,
    j: int,
    x0: float,
    y0: float,
    scale: float,
    samples: int = 14,
) -> str:
    points = geodesic_points(vertices[i], vertices[j], samples=samples)
    projected = [project(point, x0, y0, scale) for point in points]
    first = projected[0]
    rest = " ".join(f"L {x:.2f} {y:.2f}" for x, y in projected[1:])
    return f"M {first[0]:.2f} {first[1]:.2f} {rest}"


def face_fill(face_size: int) -> str:
    if face_size == 3:
        return "#e7eef8"
    if face_size == 4:
        return "#f4d6a3"
    if face_size == 5:
        return "#d9e9c7"
    return "#eadcf4"


def draw_frame(
    data: Dict[str, object],
    frame_index: int,
    x: float,
    y: float,
    cell: float,
) -> List[str]:
    vertices = data["vertices"]
    edges = data["edges"]
    faces = data.get("faces", [])
    fixed = data.get("fixed", [])
    margin = 18.0
    label_height = 20.0
    radius = 0.5 * (cell - 2.0 * margin - label_height)
    x0 = x + 0.5 * cell
    y0 = y + margin + radius
    scale = radius
    lines = [
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{cell:.2f}" height="{cell:.2f}" fill="#ffffff"/>',
        f'<circle cx="{x0:.2f}" cy="{y0:.2f}" r="{radius:.2f}" fill="#fbfcfd" stroke="#111827" stroke-width="1.1"/>',
    ]

    for face in faces:
        if len(face) > 5:
            continue
        polygon = " ".join(
            f"{px:.2f},{py:.2f}"
            for px, py in (project(vertices[idx], x0, y0, scale) for idx in face)
        )
        lines.append(
            f'<polygon points="{polygon}" fill="{face_fill(len(face))}" fill-opacity="0.36" stroke="none"/>'
        )

    for i, j in edges:
        path = path_for_geodesic(vertices, int(i), int(j), x0, y0, scale)
        lines.append(
            f'<path d="{path}" fill="none" stroke="#15803d" stroke-width="1.35" stroke-linecap="round" stroke-opacity="0.95"/>'
        )

    for k, i in enumerate(fixed):
        j = fixed[(k + 1) % len(fixed)]
        path = path_for_geodesic(vertices, int(i), int(j), x0, y0, scale, samples=24)
        lines.append(
            f'<path d="{path}" fill="none" stroke="#000000" stroke-width="1.85" stroke-dasharray="6 5" stroke-linecap="round"/>'
        )

    fixed_set = set(int(idx) for idx in fixed)
    for idx, point in enumerate(vertices):
        px, py = project(point, x0, y0, scale)
        if idx in fixed_set:
            lines.append(
                f'<circle cx="{px:.2f}" cy="{py:.2f}" r="1.95" fill="#dc2626" stroke="#7f1d1d" stroke-width="0.45"/>'
            )
        else:
            lines.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="0.9" fill="#111827"/>')

    t = frame_index / 24.0
    lines.append(
        f'<text x="{x0:.2f}" y="{y + cell - 8.0:.2f}" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#111827">'
        f'frame {frame_index}, t={t:.2f}</text>'
    )
    return lines


def write_grid_svg(frame_dir: Path, output_svg: Path, title: str) -> None:
    cell = 250.0
    title_height = 34.0
    width = 3.0 * cell
    height = title_height + 3.0 * cell
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{width:.0f}" height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2.0:.2f}" y="22" text-anchor="middle" '
        f'font-family="Helvetica,Arial,sans-serif" font-size="16" font-weight="700" fill="#111827">{title}</text>',
    ]

    for grid_index, frame_index in enumerate(FRAME_INDICES):
        row = grid_index // 3
        col = grid_index % 3
        frame = load_frame(frame_dir, frame_index)
        lines.extend(draw_frame(frame, frame_index, col * cell, title_height + row * cell, cell))

    lines.append("</svg>")
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    output_svg.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="output/paper_figures")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    for name, frame_dir, title in MORPHS:
        svg_path = output_dir / f"{name}_morph_grid.svg"
        write_grid_svg(Path(frame_dir), svg_path, title)
        print(f"wrote {svg_path}")


if __name__ == "__main__":
    main()
