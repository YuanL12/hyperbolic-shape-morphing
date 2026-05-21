#!/usr/bin/env python3
"""Morph genus-2 boundary vertices from mean-value directed weights to ones."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from directed_edge_weights import DirectedEdgeWeightCalculator
from embedding_utils import validate_disk_embedding
from poincare_harmonic_map import HarmonicMapSolver

DirectedWeight = Tuple[float, float]


def load_json(path: str) -> Dict[str, object]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError("Top-level JSON must be an object.")
    return data


def interpolate_directed_weights(
    start_weights: Sequence[DirectedWeight],
    target_weights: Sequence[DirectedWeight],
    t: float,
) -> List[DirectedWeight]:
    if len(start_weights) != len(target_weights):
        raise ValueError("start_weights and target_weights must have the same length.")
    return [
        (
            (1.0 - t) * start_i + t * target_i,
            (1.0 - t) * start_j + t * target_j,
        )
        for (start_i, start_j), (target_i, target_j) in zip(
            start_weights,
            target_weights,
        )
    ]


def average_edge_weights(directed_weights: Sequence[DirectedWeight]) -> List[float]:
    return [0.5 * (weight_i + weight_j) for weight_i, weight_j in directed_weights]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=25)
    parser.add_argument(
        "--base",
        default="input/example_genus2_boundary_vertices_embedded_1_20.json",
        help="Base unweighted mesh JSON used for vertices, faces/edges, constraints, and fixed vertices.",
    )
    parser.add_argument(
        "--output-dir",
        default="morph_frames",
        help="Directory for per-frame JSON outputs.",
    )
    parser.add_argument(
        "--normalization",
        choices=("unnormalized", "normalized"),
        default="unnormalized",
        help="Mean-value directed weights used for w0.",
    )
    parser.add_argument(
        "--directed-edge-weight-operator",
        choices=("mean_value",),
        default="mean_value",
        help="Operator used by poincare_harmonic_map.py for directed weights.",
    )
    parser.add_argument("--iterations", type=int, default=2200)
    parser.add_argument("--step-size", type=float, default=0.01)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument(
        "--skip-embedding-check",
        action="store_true",
        help="Skip validation that base vertices/faces form an embedding.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame_count = args.frames
    if frame_count < 2:
        raise ValueError("--frames must be at least 2")

    base = load_json(args.base)
    if not args.skip_embedding_check:
        embedding_report = validate_disk_embedding(base)
        if not embedding_report["is_embedding"]:
            raise ValueError(
                "Base input is not a valid disk embedding: "
                f"{embedding_report}"
            )
        print(
            "embedding check: "
            f"faces={embedding_report['face_count']} "
            f"edges={embedding_report['edge_count']} "
            "crossings=0"
        )

    initial_vertices = [list(v) for v in base["vertices"]]

    start_data = dict(base)
    start_data["vertices"] = initial_vertices
    start_data["constraints"] = base.get("constraints", [])
    start_data["fixed"] = base.get("fixed", [])

    # compute the directed edge weights from the initial vertices
    calculator = DirectedEdgeWeightCalculator(
        start_data,
        normalization=args.normalization,
    )
    start_directed_weights = calculator.compute()
    target_directed_weights = [(1.0, 1.0)] * len(start_directed_weights)

    # create the output directory
    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)

    prev_vertices = initial_vertices
    edges = [list(edge) for edge in calculator.edges]

    for frame_id in range(frame_count):
        t = frame_id / (frame_count - 1)
        current_directed_weights = interpolate_directed_weights(
            start_directed_weights,
            target_directed_weights,
            t,
        )
        current_edge_weights = average_edge_weights(current_directed_weights)

        data = dict(base)
        data["vertices"] = prev_vertices
        data["edges"] = edges
        data["edge_weights"] = current_edge_weights
        data["directed_edge_weights"] = [
            [weight_i, weight_j] for weight_i, weight_j in current_directed_weights
        ]
        data["directed_edge_weight_operator"] = args.directed_edge_weight_operator
        data["constraints"] = base.get("constraints", [])
        data["fixed"] = base.get("fixed", [])
        data["iterations"] = args.iterations
        data["step_size"] = args.step_size
        data["tolerance"] = args.tolerance

        solver = HarmonicMapSolver(data)
        result = solver.solve()
        result["morph_t"] = t
        result["edge_weights"] = current_edge_weights
        result["directed_edge_weights"] = data["directed_edge_weights"]
        result["directed_edge_weight_operator"] = args.directed_edge_weight_operator
        if "faces" in base:
            result["faces"] = base["faces"]
        result["constraints"] = base.get("constraints", [])
        result["fixed"] = base.get("fixed", [])
        result["metadata"] = {
            "morph_direction": "mean_value_directed_to_constant_directed_1",
            "directed_edge_weight_operator": args.directed_edge_weight_operator,
            "initial_base_input": args.base,
            "start_weight_normalization": args.normalization,
            "start_min_directed_weight": min(
                min(weight_i, weight_j) for weight_i, weight_j in start_directed_weights
            ),
            "start_max_directed_weight": max(
                max(weight_i, weight_j) for weight_i, weight_j in start_directed_weights
            ),
            "current_min_directed_weight": min(
                min(weight_i, weight_j)
                for weight_i, weight_j in current_directed_weights
            ),
            "current_max_directed_weight": max(
                max(weight_i, weight_j)
                for weight_i, weight_j in current_directed_weights
            ),
        }

        frame_json = out_dir / f"frame_{frame_id:03d}.json"
        with open(frame_json, "w", encoding="utf-8") as fh:
            json.dump(result, fh, indent=2)
            fh.write("\n")

        prev_vertices = result["vertices"]
        print(
            f"frame {frame_id:03d}/{frame_count - 1:03d}: "
            f"t={t:.3f} energy={float(result['energy']):.9f} "
            f"grad={float(result['mean_free_gradient_norm'] or 0.0):.3e} "
            f"max_directed_w={result['metadata']['current_max_directed_weight']:.3g}"
        )

    print(f"Wrote {out_dir}")


if __name__ == "__main__":
    main()
