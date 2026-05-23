#!/usr/bin/env python3
"""Generate an irregular genus-2 Bolza octagon Delaunay embedding."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List

from directed_edge_weights import DirectedEdgeWeightCalculator
from embedding_utils import validate_disk_embedding
from generate_genus3_klein_embedding import (
    build_edges,
    make_delaunay_vertices_and_faces,
    outside_regular_domain,
    perturb_roots,
    side_pairing_isometry,
    write_preview_svg,
)
from poincare_harmonic_map import HarmonicMapSolver, from_complex


BOLZA_PAIRINGS = [
    {"master_side": 0, "slave_side": 2},
    {"master_side": 1, "slave_side": 3},
    {"master_side": 4, "slave_side": 6},
    {"master_side": 5, "slave_side": 7},
]


def generate(args: argparse.Namespace) -> Dict[str, object]:
    sides = 8
    side_steps = args.side_subdivisions + 1
    boundary_count_expected = sides * side_steps
    vertices, faces, boundary_count = make_delaunay_vertices_and_faces(
        sides=sides,
        side_subdivisions=args.side_subdivisions,
        interior_count=args.interior_points,
        boundary_radius=args.boundary_radius,
        seed=args.seed,
    )
    if boundary_count != boundary_count_expected:
        raise RuntimeError("Unexpected boundary count.")

    fixed = [side * side_steps for side in range(sides)]
    constraints = []
    for pairing in BOLZA_PAIRINGS:
        master_side = int(pairing["master_side"])
        slave_side = int(pairing["slave_side"])
        master_start = vertices[master_side * side_steps]
        master_end = vertices[((master_side + 1) * side_steps) % boundary_count]
        slave_start = vertices[slave_side * side_steps]
        slave_end = vertices[((slave_side + 1) * side_steps) % boundary_count]
        gamma = side_pairing_isometry(master_start, master_end, slave_start, slave_end)
        for offset in range(1, side_steps):
            constraints.append(
                {
                    "slave": slave_side * side_steps + offset,
                    "master": master_side * side_steps + (side_steps - offset),
                    "a": from_complex(gamma.a),
                    "b": from_complex(gamma.b),
                }
            )

    data: Dict[str, object] = {
        "name": "genus2_bolza_octagon_delaunay_irregular",
        "description": (
            "Irregular genus-2 Bolza cut-open octagon generated from random "
            "interior points and a local Bowyer-Watson Delaunay triangulation. "
            "Boundary side vertices are constrained by disk Mobius side pairings."
        ),
        "vertices": [from_complex(z) for z in vertices],
        "faces": faces,
        "edges": build_edges(faces),
        "fixed": fixed,
        "constraints": constraints,
        "metadata": {
            "surface": "Bolza surface",
            "genus": 2,
            "fundamental_domain": "regular 8-gon",
            "boundary_radius": args.boundary_radius,
            "side_subdivisions": args.side_subdivisions,
            "interior_points": args.interior_points,
            "random_seed": args.seed,
            "boundary_cycle_size": boundary_count,
            "boundary_vertex_count": boundary_count,
            "pairings": BOLZA_PAIRINGS,
            "pairing_note": "Each master side maps to the slave side with reversed boundary orientation.",
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

    # Sanity check: every root has enough local neighbors for mean-value stars.
    DirectedEdgeWeightCalculator(data).compute()
    return data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="input/example_genus2_bolza_delaunay_irregular.json",
    )
    parser.add_argument("--side-subdivisions", type=int, default=4)
    parser.add_argument("--interior-points", type=int, default=150)
    parser.add_argument("--boundary-radius", type=float, default=0.8408964152537146)
    parser.add_argument("--seed", type=int, default=4219)
    parser.add_argument("--irregularity", type=float, default=0.04)
    parser.add_argument("--solve-iterations", type=int, default=5000)
    parser.add_argument("--solve-step-size", type=float, default=0.006)
    parser.add_argument("--solve-tolerance", type=float, default=1e-9)
    parser.add_argument(
        "--preview-svg",
        default="output/genus2_bolza_delaunay_initial_embedding.svg",
    )
    return parser.parse_args()


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
