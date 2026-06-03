#!/usr/bin/env python3
"""Generate black/white Poincare checkerboard art from a tiled disk triangulation."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from embedding_utils import poincare_geodesic_points
from poincare_harmonic_map import MobiusIsometry, to_complex


Point = Tuple[float, float]


def svg_escape(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def project(point: Sequence[float], size: int, radius: float) -> Point:
    return (
        0.5 * size + radius * float(point[0]),
        0.5 * size - radius * float(point[1]),
    )


def path_polygon(points: Sequence[Point]) -> str:
    if not points:
        return ""
    commands = [f"M {points[0][0]:.3f} {points[0][1]:.3f}"]
    commands.extend(f"L {x:.3f} {y:.3f}" for x, y in points[1:])
    commands.append("Z")
    return " ".join(commands)


def apply_mobius(transform: MobiusIsometry, point: Sequence[float]) -> Point:
    image = transform.apply(to_complex(point))
    return (image.real, image.imag)


def compose_mobius(
    left: MobiusIsometry,
    right: MobiusIsometry,
) -> MobiusIsometry:
    """Return left after right as SU(1,1) disk-isometry parameters."""
    return MobiusIsometry(
        left.a * right.a + left.b * right.b.conjugate(),
        left.a * right.b + left.b * right.a.conjugate(),
    )


def geodesic_path(
    start: Sequence[float],
    end: Sequence[float],
    *,
    size: int,
    radius: float,
    samples: int,
) -> str:
    return path_polygon(
        [
            project(point, size, radius)
            for point in poincare_geodesic_points(start, end, samples=samples)
        ]
    )[:-2]


def reference_vertices(data: Dict[str, object]) -> List[Sequence[float]]:
    reference = data.get("reference_fundamental_domain")
    if isinstance(reference, dict):
        vertices = reference.get("vertices")
        if isinstance(vertices, list) and vertices:
            return vertices

    fixed = data.get("fixed", [])
    vertices = data.get("vertices", [])
    if not isinstance(fixed, list) or not isinstance(vertices, list):
        raise ValueError("Input needs fixed vertices or reference_fundamental_domain.")
    return [vertices[int(index)] for index in fixed]


def transformed_domain_path(
    domain_vertices: Sequence[Sequence[float]],
    transform: MobiusIsometry,
    *,
    size: int,
    radius: float,
) -> str:
    projected: List[Point] = []
    transformed = [apply_mobius(transform, vertex) for vertex in domain_vertices]
    for idx, start in enumerate(transformed):
        end = transformed[(idx + 1) % len(transformed)]
        segment = [
            project(point, size, radius)
            for point in poincare_geodesic_points(start, end, samples=48)
        ]
        if projected:
            projected.extend(segment[1:])
        else:
            projected.extend(segment)
    return path_polygon(projected)


def _transform_key(transform: MobiusIsometry) -> Tuple[int, int, int, int]:
    scale = 10**9
    return (
        round(transform.a.real * scale),
        round(transform.a.imag * scale),
        round(transform.b.real * scale),
        round(transform.b.imag * scale),
    )


def _side_pairing_generators(data: Dict[str, object]) -> List[MobiusIsometry]:
    metadata = data.get("metadata", {})
    pairings = metadata.get("pairings") if isinstance(metadata, dict) else None
    constraints = data.get("constraints", [])
    if not isinstance(pairings, list) or not isinstance(constraints, list):
        return []

    side_subdivisions = int(metadata.get("side_subdivisions", 4))
    side_steps = side_subdivisions + 1
    by_slave_side: Dict[int, Dict[str, object]] = {}
    for constraint in constraints:
        if not isinstance(constraint, dict):
            continue
        slave = int(constraint["slave"])
        by_slave_side.setdefault(slave // side_steps, constraint)

    generators: List[MobiusIsometry] = []
    for pairing in pairings:
        if not isinstance(pairing, dict):
            continue
        slave_side = int(pairing["slave_side"])
        constraint = by_slave_side.get(slave_side)
        if constraint is None:
            continue
        gamma = MobiusIsometry(to_complex(constraint["a"]), to_complex(constraint["b"]))
        generators.extend([gamma, gamma.inverse()])
    return generators


def tile_transforms(data: Dict[str, object], tile_depth: int) -> List[MobiusIsometry]:
    identity = MobiusIsometry(1.0 + 0.0j, 0.0j)
    if tile_depth <= 0:
        return [identity]

    generators = _side_pairing_generators(data)
    transforms = [identity]
    seen = {_transform_key(identity)}
    frontier = [identity]
    for _ in range(tile_depth):
        next_frontier: List[MobiusIsometry] = []
        for current in frontier:
            for generator in generators:
                candidate = compose_mobius(generator, current)
                key = _transform_key(candidate)
                if key in seen:
                    continue
                seen.add(key)
                transforms.append(candidate)
                next_frontier.append(candidate)
        frontier = next_frontier
    return transforms


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="input/example_genus2_bolza_symmetric_center_fan.json",
        help="Input triangulation JSON.",
    )
    parser.add_argument(
        "--color-reference",
        help=(
            "Optional embedding whose face centers define the black/white "
            "pattern. Use this for morph frames so colors move with face IDs."
        ),
    )
    parser.add_argument(
        "--output",
        default="output/poincare_checkerboard_genus2_depth2.svg",
        help="Output SVG path.",
    )
    parser.add_argument("--size", type=int, default=1400, help="SVG width and height.")
    parser.add_argument("--tile-depth", type=int, default=2, help="Side-pairing expansion depth.")
    parser.add_argument(
        "--curve-samples",
        type=int,
        default=18,
        help="Samples per Poincare geodesic triangle edge.",
    )
    parser.add_argument(
        "--show-domain",
        action="store_true",
        help="Draw red fundamental-domain guide curves.",
    )
    parser.add_argument(
        "--pattern",
        choices=("spiral8", "ring_alternating"),
        default="spiral8",
        help="Black/white pattern used to color the triangulation faces.",
    )
    return parser.parse_args()


def load_data(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def geodesic_polygon_path(
    points: Sequence[Sequence[float]],
    *,
    size: int,
    radius: float,
    samples: int,
) -> str:
    projected: List[Point] = []
    for idx, start in enumerate(points):
        end = points[(idx + 1) % len(points)]
        segment = [
            project(point, size, radius)
            for point in poincare_geodesic_points(start, end, samples=samples)
        ]
        if projected:
            projected.extend(segment[1:])
        else:
            projected.extend(segment)
    return path_polygon(projected)


def face_center(vertices: Sequence[Sequence[float]], face: Sequence[int]) -> Point:
    return (
        sum(float(vertices[int(idx)][0]) for idx in face) / len(face),
        sum(float(vertices[int(idx)][1]) for idx in face) / len(face),
    )


def angle_sector(point: Sequence[float], sectors: int) -> int:
    theta = math.atan2(float(point[1]), float(point[0]))
    return int(math.floor((theta + math.pi) / (2.0 * math.pi) * sectors)) % sectors


def vertex_ring_and_sector(vertex: int, cycle_size: int) -> Tuple[int, int]:
    if vertex == 0:
        return 0, 0
    shifted = vertex - 1
    return shifted // cycle_size + 1, shifted % cycle_size


def ring_alternating_color(
    face: Sequence[int],
    *,
    cycle_size: int,
) -> str:
    noncenter = [int(vertex) for vertex in face if int(vertex) != 0]
    if not noncenter:
        return "#f2f2f2"
    ring_sector = [vertex_ring_and_sector(vertex, cycle_size) for vertex in noncenter]
    band = min(ring for ring, _ in ring_sector)
    sector = min(sector for _, sector in ring_sector)
    parity = (band + sector) % 2
    return "#050505" if parity == 0 else "#f2f2f2"


def spiral8_color(
    vertices: Sequence[Sequence[float]],
    face: Sequence[int],
) -> str:
    center = face_center(vertices, face)
    radius = math.hypot(center[0], center[1])
    theta = math.atan2(center[1], center[0])
    # Eightfold logarithmic-spiral bands.  The radial term rotates the black
    # wedge as it moves outward, producing lightning/arrow-like arms.
    phase = 8.0 * theta + 11.5 * radius + 0.22 * math.sin(24.0 * radius)
    wave = math.sin(phase)
    ring_gate = 0.18 * math.sin(6.0 * math.pi * radius + 0.6)
    center_gate = -0.28 if radius < 0.16 else 0.0
    return "#050505" if wave > ring_gate + center_gate else "#f2f2f2"


def face_color(
    vertices: Sequence[Sequence[float]],
    face: Sequence[int],
    *,
    pattern: str,
    cycle_size: int,
) -> str:
    if pattern == "ring_alternating":
        return ring_alternating_color(face, cycle_size=cycle_size)
    return spiral8_color(vertices, face)


def write_svg(
    data: Dict[str, object],
    output: Path,
    *,
    size: int,
    tile_depth: int,
    curve_samples: int,
    show_domain: bool,
    pattern: str,
) -> None:
    vertices = data.get("vertices")
    faces = data.get("faces")
    if not isinstance(vertices, list) or not isinstance(faces, list):
        raise ValueError("Input needs vertices and faces lists.")

    margin = 64.0
    disk_radius = 0.5 * size - margin
    output_key = "".join(ch if ch.isalnum() else "_" for ch in output.stem)
    disk_clip = f"diskClip_{output_key}"
    transforms = tile_transforms(data, tile_depth)
    domain_vertices = reference_vertices(data)
    tile_clip_ids = [f"tileClip_{output_key}_{idx}" for idx, _ in enumerate(transforms)]
    title = f"{data.get('name', 'poincare_checkerboard')} depth={tile_depth}"
    color_vertices = data.get("color_reference_vertices", vertices)
    if not isinstance(color_vertices, list):
        raise ValueError("color_reference_vertices must be a list when provided.")
    metadata = data.get("metadata", {})
    cycle_size = 40
    if isinstance(metadata, dict):
        cycle_size = int(metadata.get("boundary_cycle_size", cycle_size))

    lines: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" viewBox="0 0 {size} {size}">',
        "<defs>",
        f'<clipPath id="{disk_clip}">',
        f'<circle cx="{size / 2:.3f}" cy="{size / 2:.3f}" r="{disk_radius:.3f}"/>',
        "</clipPath>",
    ]
    for clip_id, transform in zip(tile_clip_ids, transforms):
        tile_path = transformed_domain_path(domain_vertices, transform, size=size, radius=disk_radius)
        lines.extend([f'<clipPath id="{clip_id}">', f'<path d="{tile_path}"/>', "</clipPath>"])
    lines.extend(
        [
            "</defs>",
            '<rect width="100%" height="100%" fill="#ffffff"/>',
            f'<circle cx="{size / 2:.3f}" cy="{size / 2:.3f}" r="{disk_radius:.3f}" fill="#d7d7d7"/>',
            f'<g clip-path="url(#{disk_clip})">',
        ]
    )

    sorted_faces = sorted(
        enumerate(faces),
        key=lambda item: math.hypot(*face_center(vertices, item[1])),
        reverse=True,
    )
    for transform_index, transform in enumerate(transforms):
        transformed_vertices = [apply_mobius(transform, vertex) for vertex in vertices]
        lines.append(f'<g clip-path="url(#{tile_clip_ids[transform_index]})">')
        for face_index, face in sorted_faces:
            center = face_center(transformed_vertices, face)
            if center[0] * center[0] + center[1] * center[1] > 0.995:
                continue
            face_points = [transformed_vertices[int(idx)] for idx in face]
            d = geodesic_polygon_path(
                face_points,
                size=size,
                radius=disk_radius,
                samples=curve_samples,
            )
            fill = face_color(
                color_vertices,
                face,
                pattern=pattern,
                cycle_size=cycle_size,
            )
            stroke = "#111111" if fill == "#f2f2f2" else "#eeeeee"
            lines.append(
                f'<path d="{d}" fill="{fill}" stroke="{stroke}" '
                'stroke-opacity="0.38" stroke-width="0.65" stroke-linejoin="round"/>'
            )
        lines.append("</g>")

    if show_domain:
        for transform in transforms:
            transformed = [apply_mobius(transform, vertex) for vertex in domain_vertices]
            for idx, start in enumerate(transformed):
                end = transformed[(idx + 1) % len(transformed)]
                d = geodesic_path(start, end, size=size, radius=disk_radius, samples=80)
                lines.append(
                    f'<path d="{d}" fill="none" stroke="#ef2424" stroke-width="3.2" '
                    'stroke-opacity="0.75" stroke-linecap="round"/>'
                )

    lines.extend(
        [
            "</g>",
            f'<circle cx="{size / 2:.3f}" cy="{size / 2:.3f}" r="{disk_radius:.3f}" '
            'fill="none" stroke="#111827" stroke-width="4"/>',
            f'<text x="{size / 2:.3f}" y="{size - 24:.3f}" text-anchor="middle" '
            'font-family="Avenir, Helvetica, Arial, sans-serif" font-size="24" fill="#334155">',
            svg_escape(title),
            "</text>",
            "</svg>",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    data = load_data(Path(args.input))
    if args.color_reference:
        reference = load_data(Path(args.color_reference))
        data["color_reference_vertices"] = reference["vertices"]
    write_svg(
        data,
        Path(args.output),
        size=args.size,
        tile_depth=args.tile_depth,
        curve_samples=args.curve_samples,
        show_domain=args.show_domain,
        pattern=args.pattern,
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
