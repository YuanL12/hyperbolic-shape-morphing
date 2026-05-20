#!/usr/bin/env python3
"""Morph the harder genus-2 boundary-vertex example from 1/100 weights to 1."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

from poincare_harmonic_map import HarmonicMapSolver


DEFAULT_FRAME_COUNT = 25


def load_json(path: str) -> Dict[str, object]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError("Top-level JSON must be an object.")
    return data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=DEFAULT_FRAME_COUNT)
    # parser.add_argument(
    #     "--output-dir",
    #     default="boundary_morph_1_100_frames",
    #     help="Directory for per-frame JSON outputs.",
    # )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame_count = args.frames
    if frame_count < 2:
        raise ValueError("--frames must be at least 2")

    base = load_json("example_genus2_boundary_vertices.json")
    start = load_json("example_genus2_boundary_vertices_weighted_1_20.json")
    start_sol = load_json("example_genus2_boundary_vertices_weighted_1_20_solution.json")
    output_dir = "boundary_morph_1_20_to_1_frames"

    start_weights = [float(w) for w in start["edge_weights"]]
    if "edges" not in start:
        raise ValueError("Expected explicit edge list in weighted input.")

    out_dir = Path(output_dir)
    out_dir.mkdir(exist_ok=True)

    prev_vertices = [list(v) for v in start_sol["vertices"]]

    for frame_id in range(frame_count):
        # compute the edge weights at time t
        t = frame_id / (frame_count - 1)
        current_weights = [(1.0 - t) * w + t * 1.0 for w in start_weights]

        data = dict(base)
        data["vertices"] = prev_vertices
        data["edges"] = start["edges"]
        data["edge_weights"] = current_weights
        data["constraints"] = base.get("constraints", [])
        data["fixed"] = base.get("fixed", [])
        data["iterations"] = 2200
        data["step_size"] = 0.01
        data["tolerance"] = 1e-9

        solver = HarmonicMapSolver(data)
        result = solver.solve()
        result["morph_t"] = t
        result["edge_weights"] = current_weights
        result["faces"] = base["faces"]
        result["constraints"] = base.get("constraints", [])
        result["fixed"] = base.get("fixed", [])
        result["metadata"] = {
            "morph_direction": "weighted_1_100_to_constant",
            "start_max_weight": max(start_weights),
            "current_max_weight": max(current_weights),
        }

        frame_json = out_dir / f"frame_{frame_id:03d}.json"
        with open(frame_json, "w", encoding="utf-8") as fh:
            json.dump(result, fh, indent=2)
            fh.write("\n")

        prev_vertices = result["vertices"]
        print(
            f"frame {frame_id:03d}/{frame_count - 1:03d}: "
            f"t={t:.3f} energy={float(result['energy']):.9f} "
            f"maxw={max(current_weights):.3f}"
        )

    print(f"Wrote {out_dir}")


if __name__ == "__main__":
    main()
