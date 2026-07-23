#!/usr/bin/env python3
"""Copy endpoint inputs while relaxing selected fixed-corner orbits.

The copied files keep a drawing-only reference_fundamental_domain, but remove
the selected corner orbit(s) from fixed so the harmonic-map solver can move
them.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

from hyper_morph.morph import (
    reference_fundamental_domain_from,
    with_corner_attachment_constraints,
)


DEFAULT_ENDPOINTS = [
    "examples/input/example_genus2_bolza_delaunay_all_1_solution.json",
    "examples/input/example_genus2_bolza_delaunay_half_10_solution.json",
    "examples/input/example_genus2_bolza_convex_quad_pent_graph_all_1_solution.json",
    "examples/input/example_genus2_bolza_convex_quad_pent_graph_half_10_solution.json",
    "examples/input/example_genus3_klein_quartic_triangulation_all_1_solution.json",
    "examples/input/example_genus3_klein_quartic_triangulation_half_10_solution.json",
    "examples/input/example_genus3_klein_quartic_convex_quad_pent_graph_all_1_solution.json",
    "examples/input/example_genus3_klein_quartic_convex_quad_pent_graph_half_10_solution.json",
]


def load_json(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"{path} does not contain a JSON object.")
    return data


def relaxed_copy(
    data: Dict[str, object],
    corner_positions: List[int],
) -> Dict[str, object]:
    fixed = data.get("fixed", [])
    if not isinstance(fixed, list) or not fixed:
        raise ValueError("Input has no fixed corner list.")

    reference_corner_indices: List[int] = [int(idx) for idx in fixed]
    for corner_position in corner_positions:
        if corner_position < 0 or corner_position >= len(reference_corner_indices):
            raise ValueError(
                f"--corner-positions values must be in [0, {len(reference_corner_indices) - 1}] "
                "for this input."
            )
    relaxed_vertices = [reference_corner_indices[position] for position in corner_positions]

    with_corners = with_corner_attachment_constraints(data)
    base_constraints = list(data.get("constraints", []))
    corner_constraints = list(with_corners.get("constraints", []))[len(base_constraints) :]
    component_vertices = set(relaxed_vertices)
    changed = True
    while changed:
        changed = False
        for constraint in corner_constraints:
            slave = int(constraint["slave"])
            master = int(constraint["master"])
            if slave in component_vertices or master in component_vertices:
                before = len(component_vertices)
                component_vertices.update((slave, master))
                changed = changed or len(component_vertices) != before

    component_constraints = [
        constraint
        for constraint in corner_constraints
        if int(constraint["slave"]) in component_vertices
        and int(constraint["master"]) in component_vertices
    ]
    solver_fixed = [idx for idx in reference_corner_indices if idx not in component_vertices]

    out = dict(data)
    out["fixed"] = solver_fixed
    out["constraints"] = base_constraints + component_constraints
    out["reference_fundamental_domain"] = reference_fundamental_domain_from(data)

    metadata = dict(out.get("metadata", {})) if isinstance(out.get("metadata"), dict) else {}
    metadata.update(
        {
            "relaxed_corner_vertex": relaxed_vertices[0],
            "relaxed_corner_position": corner_positions[0],
            "relaxed_corner_vertices": relaxed_vertices,
            "relaxed_corner_positions": corner_positions,
            "relaxed_corner_component_vertices": sorted(component_vertices),
            "relaxed_corner_constraint_count": len(component_constraints),
            "relaxed_corner_note": (
                "This is a copied endpoint with the selected corner orbit "
                "removed from fixed and tied by corner Mobius constraints; "
                "reference_fundamental_domain preserves the original "
                "drawing-only fundamental polygon."
            ),
        }
    )
    out["metadata"] = metadata
    return out


def output_path_for(path: Path, suffix: str) -> Path:
    if path.name.endswith("_solution.json"):
        return path.with_name(path.name.replace("_solution.json", f"{suffix}_solution.json"))
    return path.with_name(f"{path.stem}{suffix}{path.suffix}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs",
        nargs="*",
        default=DEFAULT_ENDPOINTS,
        help="Endpoint JSON files to copy. Defaults to the eight README endpoints.",
    )
    parser.add_argument(
        "--corner-position",
        type=int,
        default=None,
        help="Deprecated alias for a single value in --corner-positions.",
    )
    parser.add_argument(
        "--corner-positions",
        type=int,
        nargs="+",
        default=[0, 1],
        help="Positions in the fixed-corner list whose orbits should be relaxed.",
    )
    parser.add_argument(
        "--suffix",
        default="_relaxed_corner_orbits",
        help="Suffix inserted before _solution.json or .json.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    corner_positions = [args.corner_position] if args.corner_position is not None else args.corner_positions
    for input_name in args.inputs:
        input_path = Path(input_name)
        data = load_json(input_path)
        out = relaxed_copy(data, corner_positions)
        output_path = output_path_for(input_path, args.suffix)
        with output_path.open("w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=2)
            fh.write("\n")
        relaxed = out["metadata"]["relaxed_corner_component_vertices"]
        print(f"wrote {output_path} with relaxed corner component {relaxed}")


if __name__ == "__main__":
    main()
