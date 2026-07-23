#!/usr/bin/env python3
"""Generate a less-regular genus-3 Klein-quartic-style disk embedding."""

from __future__ import annotations

import argparse
import cmath
import json
import math
import random
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from hyper_morph import DirectedEdgeWeightCalculator, HarmonicMapSolver
from hyper_morph.embedding import validate_disk_embedding
from hyper_morph.solver import MobiusIsometry, from_complex
from hyper_morph.weights import compose_mobius


KLEIN_SIDE_LABELS = [1, 7, 3, 2, 5, 4, 7, 6, 2, 1, 4, 3, 6, 5]


def regular_polygon_interior_angle(vertex_radius: float, sides: int) -> float:
    vertices = [
        vertex_radius * complex(
            math.cos(2.0 * math.pi * k / sides),
            math.sin(2.0 * math.pi * k / sides),
        )
        for k in range(sides)
    ]
    prev_z = vertices[-1]
    z = vertices[0]
    next_z = vertices[1]

    def tangent_to_neighbor(neighbor: complex) -> complex:
        det = z.real * neighbor.imag - z.imag * neighbor.real
        rhs_z = 0.5 * (abs(z) ** 2 + 1.0)
        rhs_n = 0.5 * (abs(neighbor) ** 2 + 1.0)
        center = complex(
            (rhs_z * neighbor.imag - rhs_n * z.imag) / det,
            (z.real * rhs_n - neighbor.real * rhs_z) / det,
        )
        radius_vec = z - center
        tangent = 1j * radius_vec
        if ((z + 1e-6 * tangent) - neighbor).real ** 2 + (
            (z + 1e-6 * tangent) - neighbor
        ).imag ** 2 > ((z - 1e-6 * tangent) - neighbor).real ** 2 + (
            (z - 1e-6 * tangent) - neighbor
        ).imag ** 2:
            tangent = -tangent
        return tangent / abs(tangent)

    t_prev = tangent_to_neighbor(prev_z)
    t_next = tangent_to_neighbor(next_z)
    angle = abs(cmath.phase(t_next / t_prev))
    return min(angle, 2.0 * math.pi - angle)


def radius_for_regular_klein_14gon() -> float:
    target_angle = 2.0 * math.pi / 7.0
    lo, hi = 0.1, 0.98
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if regular_polygon_interior_angle(mid, 14) > target_angle:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def disk_translation_to_zero(z: complex) -> MobiusIsometry:
    scale = 1.0 / math.sqrt(1.0 - abs(z) ** 2)
    return MobiusIsometry(scale + 0.0j, -scale * z)


def disk_translation_from_zero(z: complex) -> MobiusIsometry:
    scale = 1.0 / math.sqrt(1.0 - abs(z) ** 2)
    return MobiusIsometry(scale + 0.0j, scale * z)


def disk_rotation(theta: float) -> MobiusIsometry:
    return MobiusIsometry(cmath.exp(0.5j * theta), 0.0j)


def side_pairing_isometry(
    master_start: complex,
    master_end: complex,
    slave_start: complex,
    slave_end: complex,
) -> MobiusIsometry:
    """Map a master side to a slave side with reversed orientation."""

    to_zero_master = disk_translation_to_zero(master_start)
    to_zero_slave = disk_translation_to_zero(slave_end)
    master_dir = to_zero_master.apply(master_end)
    slave_dir = to_zero_slave.apply(slave_start)
    theta = cmath.phase(slave_dir / master_dir)
    return compose_mobius(
        disk_translation_from_zero(slave_end),
        compose_mobius(disk_rotation(theta), to_zero_master),
    )


def build_edges(faces: Sequence[Sequence[int]]) -> List[List[int]]:
    edges = set()
    for i, j, k in faces:
        for u, v in ((i, j), (j, k), (k, i)):
            edges.add((min(u, v), max(u, v)))
    return [[i, j] for i, j in sorted(edges)]


def orient2d(a: complex, b: complex, c: complex) -> float:
    return (b.real - a.real) * (c.imag - a.imag) - (b.imag - a.imag) * (c.real - a.real)


def circumcircle_contains(a: complex, b: complex, c: complex, p: complex) -> bool:
    ax = a.real - p.real
    ay = a.imag - p.imag
    bx = b.real - p.real
    by = b.imag - p.imag
    cx = c.real - p.real
    cy = c.imag - p.imag
    det = (
        (ax * ax + ay * ay) * (bx * cy - by * cx)
        - (bx * bx + by * by) * (ax * cy - ay * cx)
        + (cx * cx + cy * cy) * (ax * by - ay * bx)
    )
    if orient2d(a, b, c) < 0.0:
        det = -det
    return det > 1e-12


def bowyer_watson(points: Sequence[complex]) -> List[List[int]]:
    min_x = min(p.real for p in points)
    max_x = max(p.real for p in points)
    min_y = min(p.imag for p in points)
    max_y = max(p.imag for p in points)
    span = max(max_x - min_x, max_y - min_y)
    center = complex(0.5 * (min_x + max_x), 0.5 * (min_y + max_y))
    super_points = [
        center + complex(-12.0 * span, -4.0 * span),
        center + complex(0.0, 13.0 * span),
        center + complex(12.0 * span, -4.0 * span),
    ]
    work = list(points) + super_points
    n = len(points)
    triangles: List[Tuple[int, int, int]] = [(n, n + 1, n + 2)]

    for idx in range(n):
        p = work[idx]
        bad = [
            tri for tri in triangles
            if circumcircle_contains(work[tri[0]], work[tri[1]], work[tri[2]], p)
        ]
        edge_count: Dict[Tuple[int, int], int] = {}
        for tri in bad:
            for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                edge = (min(u, v), max(u, v))
                edge_count[edge] = edge_count.get(edge, 0) + 1
        bad_set = set(bad)
        triangles = [tri for tri in triangles if tri not in bad_set]
        for u, v in edge_count:
            if edge_count[(u, v)] != 1:
                continue
            tri = (u, v, idx)
            if orient2d(work[tri[0]], work[tri[1]], work[tri[2]]) < 0.0:
                tri = (v, u, idx)
            triangles.append(tri)

    return [
        [i, j, k] for i, j, k in triangles
        if i < n and j < n and k < n and abs(orient2d(work[i], work[j], work[k])) > 1e-12
    ]


def point_in_convex_polygon(point: complex, polygon: Sequence[complex]) -> bool:
    return all(
        orient2d(polygon[k], polygon[(k + 1) % len(polygon)], point) >= -1e-12
        for k in range(len(polygon))
    )


def make_boundary_vertices(sides: int, side_subdivisions: int, radius: float) -> List[complex]:
    side_steps = side_subdivisions + 1
    corners = [
        radius * complex(math.cos(2.0 * math.pi * k / sides), math.sin(2.0 * math.pi * k / sides))
        for k in range(sides)
    ]
    vertices = []
    for side in range(sides):
        start = corners[side]
        end = corners[(side + 1) % sides]
        for offset in range(side_steps):
            t = offset / side_steps
            vertices.append((1.0 - t) * start + t * end)
    return vertices


def make_delaunay_vertices_and_faces(
    *,
    sides: int,
    side_subdivisions: int,
    interior_count: int,
    boundary_radius: float,
    seed: int,
) -> Tuple[List[complex], List[List[int]], int]:
    boundary = make_boundary_vertices(sides, side_subdivisions, boundary_radius)
    rng = random.Random(seed)
    min_x = min(p.real for p in boundary)
    max_x = max(p.real for p in boundary)
    min_y = min(p.imag for p in boundary)
    max_y = max(p.imag for p in boundary)
    points = list(boundary)
    attempts = 0
    while len(points) < len(boundary) + interior_count:
        attempts += 1
        if attempts > interior_count * 1000:
            raise RuntimeError("Could not generate enough interior points.")
        p = complex(rng.uniform(min_x, max_x), rng.uniform(min_y, max_y))
        if not point_in_convex_polygon(p, boundary[:: side_subdivisions + 1]):
            continue
        if min(abs(p - q) for q in points) < 0.035:
            continue
        points.append(p)

    faces = bowyer_watson(points)
    boundary_count = len(boundary)

    # The boundary points are inserted first, in cyclic order. Keep only
    # triangles inside the boundary polygon; the Delaunay hull is the 14-gon.
    clean_faces = []
    for face in faces:
        centroid = sum((points[i] for i in face), 0.0j) / 3.0
        if point_in_convex_polygon(centroid, boundary[:: side_subdivisions + 1]):
            clean_faces.append(face)
    return points, clean_faces, boundary_count


def make_ring_vertices(
    *,
    sides: int,
    side_subdivisions: int,
    rings: int,
    boundary_radius: float,
    irregularity: float,
) -> List[complex]:
    cycle_size = sides * (side_subdivisions + 1)
    vertices = [0.0j]
    for ring in range(1, rings + 1):
        rho = boundary_radius * ring / rings
        for k in range(cycle_size):
            theta = 2.0 * math.pi * k / cycle_size
            if ring != rings:
                taper = (ring / rings) * (1.0 - 0.08 * ring / rings)
                rho_k = rho * (
                    1.0
                    + irregularity
                    * taper
                    * (
                        0.35 * math.sin(3.0 * theta + 0.8 * ring)
                        + 0.22 * math.cos(5.0 * theta - 0.4 * ring)
                    )
                )
                theta_k = theta + irregularity * taper * (
                    0.20 * math.sin(4.0 * theta - 0.3 * ring)
                    - 0.14 * math.cos(7.0 * theta + 0.2 * ring)
                )
            else:
                rho_k = rho
                theta_k = theta
            vertices.append(rho_k * complex(math.cos(theta_k), math.sin(theta_k)))
    return vertices


def perturb_roots(
    data: Dict[str, object],
    *,
    irregularity: float,
) -> List[List[float]]:
    calculator = DirectedEdgeWeightCalculator(data)
    fixed = set(data["fixed"])
    seed = [complex(x, y) for x, y in data["vertices"]]

    for idx, z in enumerate(seed):
        if idx in fixed or idx not in calculator.root_index:
            continue
        x = z.real
        y = z.imag
        taper = max(0.0, 1.0 - abs(z) ** 2)
        delta = irregularity * taper * complex(
            0.58 * math.sin(4.7 * y + 0.5) + 0.23 * math.cos(6.1 * x - 0.8),
            0.46 * math.cos(5.3 * x + 0.4) - 0.26 * math.sin(4.2 * y - 1.1),
        )
        swirl = irregularity * 0.28 * taper * complex(-y, x) * math.sin(3.0 * x + 5.0 * y)
        candidate = z + delta + swirl
        if abs(candidate) >= 0.97:
            candidate *= 0.97 / abs(candidate)
        seed[idx] = candidate

    return [from_complex(z) for z in calculator.resolve_positions(seed)]


def outside_regular_domain(vertices: Sequence[Sequence[float]], fixed: Sequence[int]) -> List[int]:
    corners = [vertices[i] for i in fixed]

    def circle(a: Sequence[float], b: Sequence[float]) -> Tuple[float, float, float]:
        ax, ay = a
        bx, by = b
        det = ax * by - ay * bx
        rhs_a = 0.5 * (ax * ax + ay * ay + 1.0)
        rhs_b = 0.5 * (bx * bx + by * by + 1.0)
        cx = (rhs_a * by - rhs_b * ay) / det
        cy = (ax * rhs_b - bx * rhs_a) / det
        radius2 = (ax - cx) ** 2 + (ay - cy) ** 2
        return cx, cy, radius2

    side_circles = [circle(corners[k], corners[(k + 1) % len(corners)]) for k in range(len(corners))]
    outside = []
    for idx, (x, y) in enumerate(vertices):
        for cx, cy, radius2 in side_circles:
            origin_sign = cx * cx + cy * cy - radius2
            point_sign = (x - cx) ** 2 + (y - cy) ** 2 - radius2
            if origin_sign * point_sign < -1e-12:
                outside.append(idx)
                break
    return outside


def generate(args: argparse.Namespace) -> Dict[str, object]:
    sides = 14
    side_steps = args.side_subdivisions + 1
    cycle_size = sides * side_steps
    boundary_radius = radius_for_regular_klein_14gon()
    vertices, faces, boundary_count = make_delaunay_vertices_and_faces(
        sides=sides,
        side_subdivisions=args.side_subdivisions,
        interior_count=args.interior_points,
        boundary_radius=boundary_radius,
        seed=args.seed,
    )

    boundary_start = 0
    fixed = [side * side_steps for side in range(sides)]
    label_to_sides: Dict[int, List[int]] = {}
    for side, label in enumerate(KLEIN_SIDE_LABELS):
        label_to_sides.setdefault(label, []).append(side)

    constraints = []
    for label in sorted(label_to_sides):
        master_side, slave_side = label_to_sides[label]
        master_start = vertices[boundary_start + master_side * side_steps]
        master_end = vertices[boundary_start + ((master_side + 1) * side_steps) % cycle_size]
        slave_start = vertices[boundary_start + slave_side * side_steps]
        slave_end = vertices[boundary_start + ((slave_side + 1) * side_steps) % cycle_size]
        gamma = side_pairing_isometry(master_start, master_end, slave_start, slave_end)
        for offset in range(1, side_steps):
            slave = boundary_start + slave_side * side_steps + offset
            master = boundary_start + master_side * side_steps + (side_steps - offset)
            constraints.append(
                {
                    "slave": slave,
                    "master": master,
                    "a": from_complex(gamma.a),
                    "b": from_complex(gamma.b),
                }
            )

    data: Dict[str, object] = {
        "name": "genus3_klein_quartic_14gon_irregular",
        "description": (
            "Less-regular genus-3 Klein-quartic-style cut-open 14-gon with "
            "random interior points and local Bowyer-Watson Delaunay triangulation. "
            "Equal-numbered sides follow the standard Klein quartic figure; "
            "paired side interior vertices are constrained by disk Mobius isometries."
        ),
        "vertices": [from_complex(z) for z in vertices],
        "faces": faces,
        "edges": build_edges(faces),
        "fixed": fixed,
        "constraints": constraints,
        "metadata": {
            "surface": "Klein quartic",
            "genus": 3,
            "fundamental_domain": "regular 14-gon",
            "side_labels": KLEIN_SIDE_LABELS,
            "boundary_radius": boundary_radius,
            "side_subdivisions": args.side_subdivisions,
            "interior_points": args.interior_points,
            "random_seed": args.seed,
            "boundary_cycle_size": cycle_size,
            "boundary_vertex_count": boundary_count,
            "pairing_note": "Equal-numbered sides are associated with reversed boundary orientation.",
        },
    }

    solve_data = dict(data)
    solve_data["iterations"] = args.solve_iterations
    solve_data["step_size"] = args.solve_step_size
    solve_data["tolerance"] = args.solve_tolerance
    solution = HarmonicMapSolver(solve_data).solve()
    data["vertices"] = solution["vertices"]
    data["edges"] = solution["edges"]
    data["metadata"].update(
        {
            "regular_seed_solution_energy": solution["energy"],
            "regular_seed_solution_mean_free_gradient_norm": solution[
                "mean_free_gradient_norm"
            ],
        }
    )

    data["vertices"] = perturb_roots(data, irregularity=args.irregularity)
    report = validate_disk_embedding(data)
    outside = outside_regular_domain(data["vertices"], fixed)
    degree = [0] * len(data["vertices"])
    for i, j in data["edges"]:
        degree[i] += 1
        degree[j] += 1
    data["metadata"].update(
        {
            "irregular_embedding": True,
            "irregularity": args.irregularity,
            "max_vertex_degree": max(degree),
            "max_degree_vertices": [
                idx for idx, value in enumerate(degree) if value == max(degree)
            ][:20],
            "fundamental_domain_outside_vertex_count": len(outside),
            "fundamental_domain_outside_vertex_samples": outside[:30],
            "embedding_validation": report,
        }
    )
    if not report["is_embedding"]:
        raise ValueError(f"Generated embedding is invalid: {report}")
    return data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="examples/input/example_genus3_klein_quartic_embedded_irregular.json",
    )
    parser.add_argument("--side-subdivisions", type=int, default=4)
    parser.add_argument("--interior-points", type=int, default=180)
    parser.add_argument("--seed", type=int, default=7321)
    parser.add_argument("--irregularity", type=float, default=0.045)
    parser.add_argument("--solve-iterations", type=int, default=5000)
    parser.add_argument("--solve-step-size", type=float, default=0.006)
    parser.add_argument("--solve-tolerance", type=float, default=1e-9)
    parser.add_argument(
        "--preview-svg",
        default="examples/output/genus3_klein_initial_embedding.svg",
        help="SVG preview of the generated initial embedding.",
    )
    return parser.parse_args()


def geodesic_points(a: Sequence[float], b: Sequence[float], samples: int = 32) -> List[Tuple[float, float]]:
    det = a[0] * b[1] - a[1] * b[0]
    if abs(det) < 1e-10:
        return [(float(a[0]), float(a[1])), (float(b[0]), float(b[1]))]
    rhs_a = 0.5 * (a[0] * a[0] + a[1] * a[1] + 1.0)
    rhs_b = 0.5 * (b[0] * b[0] + b[1] * b[1] + 1.0)
    center = (
        (rhs_a * b[1] - rhs_b * a[1]) / det,
        (a[0] * rhs_b - b[0] * rhs_a) / det,
    )
    radius = math.hypot(a[0] - center[0], a[1] - center[1])
    start = math.atan2(a[1] - center[1], a[0] - center[0])
    end = math.atan2(b[1] - center[1], b[0] - center[0])
    delta = end - start
    while delta <= -math.pi:
        delta += 2.0 * math.pi
    while delta > math.pi:
        delta -= 2.0 * math.pi
    mid = (
        center[0] + radius * math.cos(start + 0.5 * delta),
        center[1] + radius * math.sin(start + 0.5 * delta),
    )
    if mid[0] * mid[0] + mid[1] * mid[1] > 1.0 + 1e-7:
        delta += -2.0 * math.pi if delta > 0.0 else 2.0 * math.pi
    return [
        (
            center[0] + radius * math.cos(start + delta * k / samples),
            center[1] + radius * math.sin(start + delta * k / samples),
        )
        for k in range(samples + 1)
    ]


def write_preview_svg(data: Dict[str, object], output: Path) -> None:
    vertices = data["vertices"]
    edges = data["edges"]
    fixed = data["fixed"]
    size = 1000
    margin = 64
    scale = (size - 2 * margin) / 2.0

    def project(point: Sequence[float]) -> Tuple[float, float]:
        return (
            size / 2.0 + point[0] * scale,
            size / 2.0 - point[1] * scale,
        )

    def path_for(points: Sequence[Tuple[float, float]]) -> str:
        projected = [project(p) for p in points]
        first = projected[0]
        rest = " ".join(f"L {x:.3f} {y:.3f}" for x, y in projected[1:])
        return f"M {first[0]:.3f} {first[1]:.3f} {rest}"

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="1000" viewBox="0 0 1000 1000">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<circle cx="{size/2}" cy="{size/2}" r="{scale}" fill="#fbfcfd" stroke="#111827" stroke-width="2"/>',
    ]
    for i, j in edges:
        d = path_for(geodesic_points(vertices[i], vertices[j], samples=20))
        lines.append(f'<path d="{d}" fill="none" stroke="#188038" stroke-width="1.15" stroke-linecap="round"/>')
    for k, i in enumerate(fixed):
        j = fixed[(k + 1) % len(fixed)]
        d = path_for(geodesic_points(vertices[i], vertices[j], samples=36))
        lines.append(f'<path d="{d}" fill="none" stroke="#000000" stroke-width="2.2" stroke-dasharray="8 7" stroke-linecap="round"/>')
    fixed_set = set(fixed)
    for idx, point in enumerate(vertices):
        x, y = project(point)
        if idx in fixed_set:
            lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="4.2" fill="#dc2626" stroke="#7f1d1d" stroke-width="1"/>')
        else:
            lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="1.9" fill="#111827"/>')
    lines.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    data = generate(args)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    preview = Path(args.preview_svg)
    write_preview_svg(data, preview)
    report = data["metadata"]["embedding_validation"]
    print(
        f"wrote {output}: vertices={len(data['vertices'])} "
        f"faces={len(data['faces'])} edges={len(data['edges'])} "
        f"constraints={len(data['constraints'])}"
    )
    print(f"preview: {preview}")
    print(
        "embedding check: "
        f"crossings={report['euclidean_chord_crossing_count']} "
        f"geodesic_crossings={report['poincare_geodesic_crossing_count']} "
        f"outside_regular_domain={data['metadata']['fundamental_domain_outside_vertex_count']} "
        f"max_degree={data['metadata']['max_vertex_degree']}"
    )


if __name__ == "__main__":
    main()
