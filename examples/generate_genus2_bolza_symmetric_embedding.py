#!/usr/bin/env python3
"""Generate a symmetric genus-2 Bolza triangulation with a central fan vertex."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from hyper_morph.embedding import poincare_geodesic_points, validate_disk_embedding
from hyper_morph.solver import from_complex

from generate_genus2_bolza_embedding import BOLZA_PAIRINGS
from generate_genus3_klein_embedding import build_edges, side_pairing_isometry


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="examples/input/example_genus2_bolza_symmetric_center_fan.json",
        help="Output JSON path.",
    )
    parser.add_argument("--rings", type=int, default=3, help="Number of rings outside the center.")
    parser.add_argument(
        "--side-subdivisions",
        type=int,
        default=4,
        help="Interior subdivision vertices per octagon side.",
    )
    parser.add_argument(
        "--boundary-radius",
        type=float,
        default=0.8408964152537146,
        help="Euclidean radius of the regular Bolza octagon vertices.",
    )
    return parser.parse_args()


def ring_vertices(
    *,
    sides: int,
    side_steps: int,
    rings: int,
    boundary_radius: float,
) -> List[complex]:
    cycle_size = sides * side_steps
    corners = [
        boundary_radius * complex(
            math.cos(2.0 * math.pi * idx / sides),
            math.sin(2.0 * math.pi * idx / sides),
        )
        for idx in range(sides)
    ]
    boundary: List[complex] = []
    for side in range(sides):
        samples = poincare_geodesic_points(
            (corners[side].real, corners[side].imag),
            (corners[(side + 1) % sides].real, corners[(side + 1) % sides].imag),
            samples=side_steps,
        )
        boundary.extend(complex(x, y) for x, y in samples[:-1])
    if len(boundary) != cycle_size:
        raise RuntimeError("Unexpected boundary cycle size.")

    vertices = [0.0j]
    for ring in range(1, rings + 1):
        t = ring / rings
        for point in boundary:
            vertices.append(t * point)
    return vertices


def ring_index(ring: int, idx: int, cycle_size: int) -> int:
    if ring == 0:
        return 0
    return 1 + (ring - 1) * cycle_size + idx % cycle_size


def ring_faces(*, rings: int, cycle_size: int) -> List[List[int]]:
    faces: List[List[int]] = []
    for idx in range(cycle_size):
        faces.append([0, ring_index(1, idx, cycle_size), ring_index(1, idx + 1, cycle_size)])
    for ring in range(1, rings):
        for idx in range(cycle_size):
            inner_a = ring_index(ring, idx, cycle_size)
            inner_b = ring_index(ring, idx + 1, cycle_size)
            outer_a = ring_index(ring + 1, idx, cycle_size)
            outer_b = ring_index(ring + 1, idx + 1, cycle_size)
            faces.append([inner_a, outer_a, outer_b])
            faces.append([inner_a, outer_b, inner_b])
    return faces


def reference_domain(boundary_radius: float) -> Dict[str, object]:
    vertices = [
        [
            boundary_radius * math.cos(2.0 * math.pi * idx / 8.0),
            boundary_radius * math.sin(2.0 * math.pi * idx / 8.0),
        ]
        for idx in range(8)
    ]
    return {
        "name": "bolza_regular_octagon",
        "surface": "Bolza surface",
        "polygon": "regular octagon",
        "model": "poincare_disk",
        "corner_indices": list(range(8)),
        "vertices": vertices,
        "description": "Drawing-only reference fundamental polygon.",
    }


def pairing_constraints(
    vertices: Sequence[complex],
    *,
    boundary_start: int,
    side_steps: int,
    boundary_count: int,
) -> List[Dict[str, object]]:
    constraints: List[Dict[str, object]] = []
    for pairing in BOLZA_PAIRINGS:
        master_side = int(pairing["master_side"])
        slave_side = int(pairing["slave_side"])
        master_start = vertices[boundary_start + master_side * side_steps]
        master_end = vertices[boundary_start + ((master_side + 1) * side_steps) % boundary_count]
        slave_start = vertices[boundary_start + slave_side * side_steps]
        slave_end = vertices[boundary_start + ((slave_side + 1) * side_steps) % boundary_count]
        gamma = side_pairing_isometry(master_start, master_end, slave_start, slave_end)
        for offset in range(1, side_steps):
            constraints.append(
                {
                    "slave": boundary_start + slave_side * side_steps + offset,
                    "master": boundary_start + master_side * side_steps + (side_steps - offset),
                    "a": from_complex(gamma.a),
                    "b": from_complex(gamma.b),
                }
            )
    return constraints


def generate(args: argparse.Namespace) -> Dict[str, object]:
    sides = 8
    side_steps = args.side_subdivisions + 1
    cycle_size = sides * side_steps
    vertices = ring_vertices(
        sides=sides,
        side_steps=side_steps,
        rings=args.rings,
        boundary_radius=args.boundary_radius,
    )
    faces = ring_faces(rings=args.rings, cycle_size=cycle_size)
    boundary_start = ring_index(args.rings, 0, cycle_size)
    fixed = [boundary_start + side * side_steps for side in range(sides)]
    data: Dict[str, object] = {
        "name": "genus2_bolza_symmetric_center_fan",
        "description": (
            "Symmetric genus-2 Bolza octagon triangulation with concentric rings "
            "and a high-degree central fan vertex."
        ),
        "vertices": [from_complex(vertex) for vertex in vertices],
        "faces": faces,
        "edges": build_edges(faces),
        "fixed": fixed,
        "constraints": pairing_constraints(
            vertices,
            boundary_start=boundary_start,
            side_steps=side_steps,
            boundary_count=cycle_size,
        ),
        "reference_fundamental_domain": reference_domain(args.boundary_radius),
        "metadata": {
            "surface": "Bolza surface",
            "genus": 2,
            "fundamental_domain": "regular 8-gon",
            "graph_construction": "concentric ring triangulation with center fan",
            "boundary_radius": args.boundary_radius,
            "side_subdivisions": args.side_subdivisions,
            "rings": args.rings,
            "boundary_cycle_size": cycle_size,
            "boundary_vertex_count": cycle_size,
            "central_vertex": 0,
            "central_vertex_degree": cycle_size,
            "pairings": BOLZA_PAIRINGS,
        },
    }
    data["metadata"]["embedding_validation"] = validate_disk_embedding(data)
    return data


def main() -> None:
    args = parse_args()
    data = generate(args)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    report = data["metadata"]["embedding_validation"]
    print(
        f"wrote {output}: vertices={len(data['vertices'])} "
        f"faces={len(data['faces'])} edges={len(data['edges'])} "
        f"center_degree={data['metadata']['central_vertex_degree']} "
        f"geodesic_crossings={report['poincare_geodesic_crossing_count']}"
    )


if __name__ == "__main__":
    main()
