#!/usr/bin/env python3
"""Generate paper figures showing selected morph frames as PDF grids."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from generate_genus3_klein_embedding import geodesic_points


FRAME_INDICES = [0, 3, 6, 9, 12, 15, 18, 21, 24]
FRAME_ROOT = Path("output/frames")

MORPHS_MEAN_VALUE_FIXED_CORNER = [
    (
        "genus2_triangulation_mean_value_fixed_corner",
        "boundary_morph_genus2_bolza_triangulation_mean_value_fixed_corner_frames",
        "Genus 2 triangulation, mean-value, fixed corner",
    ),
    (
        "genus2_cells_mean_value_fixed_corner",
        "boundary_morph_genus2_bolza_cells_mean_value_fixed_corner_frames",
        "Genus 2 cells, mean-value, fixed corner",
    ),
    (
        "genus3_triangulation_mean_value_fixed_corner",
        "boundary_morph_genus3_triangulation_mean_value_fixed_corner_frames",
        "Genus 3 triangulation, mean-value, fixed corner",
    ),
    (
        "genus3_cells_mean_value_fixed_corner",
        "boundary_morph_genus3_cells_mean_value_fixed_corner_frames",
        "Genus 3 cells, mean-value, fixed corner",
    ),
]

MORPHS_MEAN_VALUE_RELAXED_CORNER = [
    (
        "genus2_triangulation_mean_value_relaxed_corner",
        "boundary_morph_genus2_bolza_triangulation_mean_value_relaxed_corner_frames",
        "Genus 2 triangulation, mean-value, relaxed corner",
    ),
    (
        "genus2_cells_mean_value_relaxed_corner",
        "boundary_morph_genus2_bolza_cells_mean_value_relaxed_corner_frames",
        "Genus 2 cells, mean-value, relaxed corner",
    ),
    (
        "genus3_triangulation_mean_value_relaxed_corner",
        "boundary_morph_genus3_triangulation_mean_value_relaxed_corner_frames",
        "Genus 3 triangulation, mean-value, relaxed corner",
    ),
    (
        "genus3_cells_mean_value_relaxed_corner",
        "boundary_morph_genus3_cells_mean_value_relaxed_corner_frames",
        "Genus 3 cells, mean-value, relaxed corner",
    ),
]


def load_frame(frame_dir: Path, frame_index: int) -> Dict[str, object]:
    path = frame_dir / f"frame_{frame_index:03d}.json"
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def project(point: Sequence[float], x0: float, y0: float, scale: float) -> Tuple[float, float]:
    return x0 + point[0] * scale, y0 - point[1] * scale


def face_fill(face_size: int) -> str:
    if face_size == 3:
        return "#e7eef8"
    if face_size == 4:
        return "#f4d6a3"
    if face_size == 5:
        return "#d9e9c7"
    return "#eadcf4"


def hex_to_rgb(color: str) -> Tuple[float, float, float]:
    value = color.lstrip("#")
    return tuple(int(value[i : i + 2], 16) / 255.0 for i in (0, 2, 4))  # type: ignore[return-value]


def pdf_y(y: float, height: float) -> float:
    return height - y


def pdf_graphics_state_name(opacity: float) -> str:
    return f"GS{int(round(opacity * 100)):03d}"


def pdf_polyline(
    points: Sequence[Tuple[float, float]],
    height: float,
    stroke: str,
    width: float,
    opacity: float = 1.0,
) -> str:
    if not points:
        return ""
    r, g, b = hex_to_rgb(stroke)
    lines = [
        "q",
        f"/{pdf_graphics_state_name(opacity)} gs",
        "1 J 1 j",
        f"{width:.3f} w",
        f"{r:.4f} {g:.4f} {b:.4f} RG",
        f"{points[0][0]:.3f} {pdf_y(points[0][1], height):.3f} m",
    ]
    lines.extend(f"{x:.3f} {pdf_y(y, height):.3f} l" for x, y in points[1:])
    lines.extend(["S", "Q"])
    return "\n".join(lines) + "\n"


def pdf_polygon(
    points: Sequence[Tuple[float, float]],
    height: float,
    fill: str,
    opacity: float,
) -> str:
    if not points:
        return ""
    r, g, b = hex_to_rgb(fill)
    lines = [
        "q",
        f"/{pdf_graphics_state_name(opacity)} gs",
        f"{r:.4f} {g:.4f} {b:.4f} rg",
        f"{points[0][0]:.3f} {pdf_y(points[0][1], height):.3f} m",
    ]
    lines.extend(f"{x:.3f} {pdf_y(y, height):.3f} l" for x, y in points[1:])
    lines.extend(["h f", "Q"])
    return "\n".join(lines) + "\n"


def pdf_circle_path(cx: float, cy: float, radius: float, height: float) -> str:
    c = 0.5522847498307936 * radius
    y = pdf_y(cy, height)
    return "\n".join(
        [
            f"{cx + radius:.3f} {y:.3f} m",
            f"{cx + radius:.3f} {y + c:.3f} {cx + c:.3f} {y + radius:.3f} {cx:.3f} {y + radius:.3f} c",
            f"{cx - c:.3f} {y + radius:.3f} {cx - radius:.3f} {y + c:.3f} {cx - radius:.3f} {y:.3f} c",
            f"{cx - radius:.3f} {y - c:.3f} {cx - c:.3f} {y - radius:.3f} {cx:.3f} {y - radius:.3f} c",
            f"{cx + c:.3f} {y - radius:.3f} {cx + radius:.3f} {y - c:.3f} {cx + radius:.3f} {y:.3f} c",
            "h",
        ]
    )


def pdf_circle(
    cx: float,
    cy: float,
    radius: float,
    height: float,
    fill: str,
    stroke: str | None = None,
    stroke_width: float = 1.0,
) -> str:
    fr, fg, fb = hex_to_rgb(fill)
    lines = ["q", f"{fr:.4f} {fg:.4f} {fb:.4f} rg"]
    if stroke is not None:
        sr, sg, sb = hex_to_rgb(stroke)
        lines.extend([f"{sr:.4f} {sg:.4f} {sb:.4f} RG", f"{stroke_width:.3f} w"])
    lines.append(pdf_circle_path(cx, cy, radius, height))
    lines.append("B" if stroke is not None else "f")
    lines.append("Q")
    return "\n".join(lines) + "\n"


def pdf_text(
    text: str,
    x: float,
    y: float,
    height: float,
    size: float = 10.0,
    color: str = "#111827",
) -> str:
    escaped = text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
    r, g, b = hex_to_rgb(color)
    # Approximate center alignment with Helvetica's average digit/letter width.
    x_adjusted = x - 0.25 * size * len(text)
    return (
        "q\n"
        f"{r:.4f} {g:.4f} {b:.4f} rg\n"
        "BT\n"
        f"/F1 {size:.1f} Tf\n"
        f"{x_adjusted:.3f} {pdf_y(y, height):.3f} Td\n"
        f"({escaped}) Tj\n"
        "ET\n"
        "Q\n"
    )


def reference_domain(data: Dict[str, object]) -> Tuple[List[int], List[Sequence[float]]]:
    domain = data.get("reference_fundamental_domain")
    if not isinstance(domain, dict):
        raise ValueError("Frame is missing reference_fundamental_domain.")
    corner_indices = domain.get("corner_indices")
    domain_vertices = domain.get("vertices")
    if not isinstance(corner_indices, list) or not isinstance(domain_vertices, list) or not domain_vertices:
        raise ValueError(
            "reference_fundamental_domain must contain corner_indices and vertices."
        )
    return [int(idx) for idx in corner_indices], domain_vertices


def draw_frame_pdf(
    data: Dict[str, object],
    frame_index: int,
    x: float,
    y: float,
    cell: float,
    page_height: float,
) -> str:
    vertices = data["vertices"]
    edges = data["edges"]
    faces = data.get("faces", [])
    reference_indices, reference_vertices = reference_domain(data)
    boundary_count = 0
    if len(reference_indices) >= 2:
        boundary_count = len(reference_indices) * (reference_indices[1] - reference_indices[0])

    margin = 18.0
    label_height = 20.0
    radius = 0.5 * (cell - 2.0 * margin - label_height)
    x0 = x + 0.5 * cell
    y0 = y + margin + radius
    scale = radius
    chunks = [
        pdf_polygon(
            [(x, y), (x + cell, y), (x + cell, y + cell), (x, y + cell)],
            page_height,
            "#ffffff",
            1.0,
        ),
        pdf_circle(x0, y0, radius, page_height, "#fbfcfd", "#374151", 0.75),
    ]

    for face in faces:
        if len(face) > 5:
            continue
        points = [project(vertices[idx], x0, y0, scale) for idx in face]
        chunks.append(pdf_polygon(points, page_height, face_fill(len(face)), 0.36))

    for i, j in edges:
        points = [
            project(point, x0, y0, scale)
            for point in geodesic_points(vertices[int(i)], vertices[int(j)], samples=14)
        ]
        chunks.append(pdf_polyline(points, page_height, "#15803d", 1.25, 0.85))

    for k, point in enumerate(reference_vertices):
        next_point = reference_vertices[(k + 1) % len(reference_vertices)]
        points = [
            project(p, x0, y0, scale)
            for p in geodesic_points(point, next_point, samples=24)
        ]
        chunks.append(pdf_polyline(points, page_height, "#4b5563", 0.95, 0.62))

    for idx, point in enumerate(vertices):
        px, py = project(point, x0, y0, scale)
        fill = "#dc2626" if idx < boundary_count else "#111827"
        chunks.append(pdf_circle(px, py, 0.9, page_height, fill))

    t = frame_index / 24.0
    chunks.append(
        pdf_text(
            f"frame {frame_index}, t={t:.2f}",
            x0,
            y + cell - 8.0,
            page_height,
            size=10.0,
        )
    )
    return "".join(chunks)


def write_pdf(objects: Sequence[bytes], output_pdf: Path) -> None:
    header = b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n"
    offsets = []
    body = bytearray()
    for object_id, data in enumerate(objects, start=1):
        offsets.append(len(header) + len(body))
        body.extend(f"{object_id} 0 obj\n".encode("ascii"))
        body.extend(data)
        body.extend(b"\nendobj\n")
    xref_offset = len(header) + len(body)
    xref = bytearray()
    xref.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    xref.extend(b"0000000000 65535 f \n")
    for offset in offsets:
        xref.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
    trailer = (
        f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
        f"startxref\n{xref_offset}\n%%EOF\n"
    ).encode("ascii")
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output_pdf.write_bytes(header + body + xref + trailer)


def write_grid_pdf(frame_dir: Path, output_pdf: Path, title: str) -> None:
    cell = 250.0
    width = 3.0 * cell
    height = 3.0 * cell
    content = [
        pdf_polygon(
            [(0.0, 0.0), (width, 0.0), (width, height), (0.0, height)],
            height,
            "#ffffff",
            1.0,
        )
    ]
    for grid_index, frame_index in enumerate(FRAME_INDICES):
        row = grid_index // 3
        col = grid_index % 3
        frame = load_frame(frame_dir, frame_index)
        content.append(draw_frame_pdf(frame, frame_index, col * cell, row * cell, cell, height))

    stream = "".join(content).encode("ascii")
    alpha_objects = [
        (36, b"<< /Type /ExtGState /CA 0.36 /ca 0.36 >>"),
        (62, b"<< /Type /ExtGState /CA 0.62 /ca 0.62 >>"),
        (85, b"<< /Type /ExtGState /CA 0.85 /ca 0.85 >>"),
        (100, b"<< /Type /ExtGState /CA 1 /ca 1 >>"),
    ]
    ext_resources = " ".join(
        f"/{pdf_graphics_state_name(value / 100.0)} {object_id} 0 R"
        for object_id, (value, _) in enumerate(alpha_objects, start=6)
    )
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width:.0f} {height:.0f}] "
            f"/Resources << /Font << /F1 5 0 R >> /ExtGState << {ext_resources} >> >> "
            "/Contents 4 0 R >>"
        ).encode("ascii"),
        b"<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n" + stream + b"\nendstream",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    ]
    objects.extend(data for _, data in alpha_objects)
    write_pdf(objects, output_pdf)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="output/paper_figures")
    parser.add_argument(
        "--variant",
        choices=(
            "mean-value-fixed_corner",
            "mean-value-relaxed_corner",
        ),
        default="mean-value-fixed_corner",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    if args.variant == "mean-value-relaxed_corner":
        morphs = MORPHS_MEAN_VALUE_RELAXED_CORNER
    else:
        morphs = MORPHS_MEAN_VALUE_FIXED_CORNER
    for name, frame_dir, title in morphs:
        pdf_path = output_dir / f"{name}_morph_grid.pdf"
        write_grid_pdf(FRAME_ROOT / frame_dir, pdf_path, title)
        print(f"wrote {pdf_path}")


if __name__ == "__main__":
    main()
