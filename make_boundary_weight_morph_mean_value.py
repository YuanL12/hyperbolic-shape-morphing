#!/usr/bin/env python3
"""Morph boundary vertices between directed edge-weight fields."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from directed_edge_weights import DirectedEdgeWeightCalculator
from embedding_utils import validate_disk_embedding
from poincare_harmonic_map import MobiusIsometry, build_edges, from_complex, to_complex
from poincare_harmonic_map import HarmonicMapSolver

DirectedWeight = Tuple[float, float]
KLEIN_SIDE_LABELS = [1, 7, 3, 2, 5, 4, 7, 6, 2, 1, 4, 3, 6, 5]


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
        "--low-valence-star-policy",
        choices=("error", "unit"),
        default="error",
        help=(
            "Policy for mean-value stars with fewer than three attached "
            "neighbors. 'unit' assigns weight 1 to the available outgoing "
            "edge directions."
        ),
    )
    parser.add_argument(
        "--directed-edge-weight-operator",
        choices=("mean_value", "energy"),
        default="mean_value",
        help="Operator used by poincare_harmonic_map.py for directed weights.",
    )
    parser.add_argument(
        "--start-directed-weights",
        choices=("mean_value", "ones"),
        default="mean_value",
        help="Directed weights used at frame 0.",
    )
    parser.add_argument(
        "--target-directed-weights",
        choices=("ones", "half_10", "half_oriented_10", "mean_value_from_target"),
        default="ones",
        help=(
            "Directed weights used at the last frame. 'half_10' sets half the "
            "edges to (10, 10), half to (1, 1). 'half_oriented_10' sets half "
            "to (10, 1), half to (1, 10). 'mean_value_from_target' computes "
            "mean-value directed weights from --target-embedding."
        ),
    )
    parser.add_argument(
        "--target-embedding",
        help=(
            "Embedding JSON used when --target-directed-weights is "
            "'mean_value_from_target'. Topology is taken from --base."
        ),
    )
    parser.add_argument("--iterations", type=int, default=2200)
    parser.add_argument("--step-size", type=float, default=0.01)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument(
        "--line-search-objective",
        choices=("energy", "gradient_norm", "none"),
        default=None,
        help=(
            "Line-search objective passed to the solver. Defaults to the "
            "solver's operator-dependent choice."
        ),
    )
    parser.add_argument(
        "--convergence-criterion",
        choices=("gradient_norm", "relative_step"),
        default=None,
        help="Convergence criterion passed to the solver.",
    )
    parser.add_argument(
        "--skip-embedding-check",
        action="store_true",
        help="Skip validation that base vertices/faces form an embedding.",
    )
    parser.add_argument(
        "--skip-face-orientation-check",
        action="store_true",
        help="For non-triangle/coarsened maps, validate edge crossings only.",
    )
    return parser.parse_args()


def make_target_directed_weights(
    mode: str,
    edge_count: int,
) -> List[DirectedWeight]:
    if mode == "mean_value_from_target":
        raise ValueError(
            "mean_value_from_target requires --target-embedding and is handled "
            "in main()."
        )
    if mode == "ones":
        return [(1.0, 1.0)] * edge_count

    half_edge_count = edge_count // 2
    if mode == "half_10":
        return [
            (10.0, 10.0) if edge_idx < half_edge_count else (1.0, 1.0)
            for edge_idx in range(edge_count)
        ]
    if mode == "half_oriented_10":
        return [
            (10.0, 1.0) if edge_idx < half_edge_count else (1.0, 10.0)
            for edge_idx in range(edge_count)
        ]
    raise ValueError(f"Unknown target directed-weight mode: {mode!r}")


def with_corner_attachment_constraints(data: Dict[str, object]) -> Dict[str, object]:
    metadata = data.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}
    side_labels = metadata.get("side_labels")
    pairings = metadata.get("pairings")
    side_subdivisions = metadata.get("side_subdivisions")
    if not isinstance(side_labels, list) and not isinstance(pairings, list):
        fixed = data.get("fixed", [])
        if isinstance(fixed, list) and len(fixed) == 14:
            side_labels = KLEIN_SIDE_LABELS
        else:
            return data
    if side_subdivisions is None:
        fixed = data.get("fixed", [])
        if isinstance(fixed, list) and len(fixed) >= 2:
            side_subdivisions = int(fixed[1]) - int(fixed[0]) - 1
        else:
            return data

    side_steps = int(side_subdivisions) + 1
    side_count = (
        len(side_labels)
        if isinstance(side_labels, list)
        else len(data.get("fixed", []))
    )
    if side_count == 0:
        return data
    boundary_count = side_count * side_steps
    constraints = list(data.get("constraints", []))

    by_slave_side: Dict[int, Dict[str, object]] = {}
    for constraint in constraints:
        slave = int(constraint["slave"])
        if slave < boundary_count:
            by_slave_side.setdefault(slave // side_steps, constraint)

    transform_graph: Dict[int, List[Tuple[int, MobiusIsometry]]] = {}
    pairing_sides: List[Tuple[int, int]] = []
    if isinstance(pairings, list):
        for pairing in pairings:
            if not isinstance(pairing, dict):
                continue
            pairing_sides.append(
                (int(pairing["master_side"]), int(pairing["slave_side"]))
            )
    elif isinstance(side_labels, list):
        for label in sorted(set(int(value) for value in side_labels)):
            paired_sides = [
                side for side, side_label in enumerate(side_labels)
                if int(side_label) == label
            ]
            if len(paired_sides) != 2:
                continue
            pairing_sides.append((paired_sides[0], paired_sides[1]))

    for master_side, slave_side in pairing_sides:
        side_constraint = by_slave_side.get(slave_side)
        if side_constraint is None:
            continue
        gamma = MobiusIsometry(
            to_complex(side_constraint["a"]),
            to_complex(side_constraint["b"]),
        )
        inverse = gamma.inverse()
        master_start = master_side * side_steps
        master_end = ((master_side + 1) * side_steps) % boundary_count
        slave_start = slave_side * side_steps
        slave_end = ((slave_side + 1) * side_steps) % boundary_count

        for source, target, transform in (
            (master_start, slave_end, gamma),
            (slave_end, master_start, inverse),
            (master_end, slave_start, gamma),
            (slave_start, master_end, inverse),
        ):
            transform_graph.setdefault(source, []).append((target, transform))

    existing_slaves = {int(constraint["slave"]) for constraint in constraints}
    visited = set()
    corner_constraints = []
    for root in sorted(transform_graph):
        if root in visited:
            continue
        visited.add(root)
        queue = [root]
        while queue:
            current = queue.pop(0)
            for neighbor, transform in transform_graph.get(current, []):
                if neighbor in visited:
                    continue
                visited.add(neighbor)
                queue.append(neighbor)
                if neighbor not in existing_slaves:
                    existing_slaves.add(neighbor)
                    corner_constraints.append(
                        {
                            "slave": neighbor,
                            "master": current,
                            "a": from_complex(transform.a),
                            "b": from_complex(transform.b),
                        }
                    )

    if not corner_constraints:
        return data

    out = dict(data)
    out["constraints"] = constraints + corner_constraints
    return out


def main() -> None:
    args = parse_args()
    frame_count = args.frames
    if frame_count < 2:
        raise ValueError("--frames must be at least 2")

    base = load_json(args.base)
    if not args.skip_embedding_check:
        check_data = dict(base)
        if args.skip_face_orientation_check:
            check_data.pop("faces", None)
        embedding_report = validate_disk_embedding(check_data)
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
    edges = [
        list(edge) for edge in build_edges(
            faces=base.get("faces"),
            edges=base.get("edges"),
            n_vertices=len(initial_vertices),
        )
    ]

    start_data = dict(base)
    start_data["vertices"] = initial_vertices
    start_data["edges"] = edges
    start_data["constraints"] = base.get("constraints", [])
    start_data["fixed"] = base.get("fixed", [])
    weight_start_data = with_corner_attachment_constraints(start_data)

    if args.start_directed_weights == "mean_value":
        calculator = DirectedEdgeWeightCalculator(
            weight_start_data,
            normalization=args.normalization,
            low_valence_policy=args.low_valence_star_policy,
        )
        start_directed_weights = calculator.compute()
    else:
        start_directed_weights = [(1.0, 1.0)] * len(edges)

    if args.target_directed_weights == "mean_value_from_target":
        if not args.target_embedding:
            raise ValueError(
                "--target-embedding is required when --target-directed-weights "
                "is mean_value_from_target."
            )
        target_embedding = load_json(args.target_embedding)
        if len(target_embedding.get("vertices", [])) != len(initial_vertices):
            raise ValueError(
                "--target-embedding must have the same number of vertices as --base."
            )
        target_data = dict(base)
        target_data["vertices"] = [list(v) for v in target_embedding["vertices"]]
        target_data["edges"] = edges
        target_data["constraints"] = base.get("constraints", [])
        target_data["fixed"] = base.get("fixed", [])
        weight_target_data = with_corner_attachment_constraints(target_data)
        target_calculator = DirectedEdgeWeightCalculator(
            weight_target_data,
            normalization=args.normalization,
            low_valence_policy=args.low_valence_star_policy,
        )
        target_directed_weights = target_calculator.compute()
    else:
        target_directed_weights = make_target_directed_weights(
            args.target_directed_weights,
            len(start_directed_weights),
        )

    # create the output directory
    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)

    prev_vertices = initial_vertices

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
        if args.line_search_objective is not None:
            data["line_search_objective"] = args.line_search_objective
        if args.convergence_criterion is not None:
            data["convergence_criterion"] = args.convergence_criterion

        solver = HarmonicMapSolver(data)
        result = solver.solve()
        result["morph_t"] = t
        result["edge_weights"] = current_edge_weights
        result["directed_edge_weights"] = data["directed_edge_weights"]
        result["directed_edge_weight_operator"] = args.directed_edge_weight_operator
        result["line_search_objective"] = data.get("line_search_objective")
        result["convergence_criterion"] = data.get("convergence_criterion")
        if "faces" in base:
            result["faces"] = base["faces"]
        result["constraints"] = base.get("constraints", [])
        result["fixed"] = base.get("fixed", [])
        result["metadata"] = {
            "morph_direction": (
                f"{args.start_directed_weights}_directed_to_"
                f"{args.target_directed_weights}"
            ),
            "directed_edge_weight_operator": args.directed_edge_weight_operator,
            "initial_base_input": args.base,
            "start_weight_normalization": args.normalization,
            "low_valence_star_policy": args.low_valence_star_policy,
            "attach_corner_stars": True,
            "start_directed_weights": args.start_directed_weights,
            "target_directed_weights": args.target_directed_weights,
            "target_embedding": args.target_embedding,
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
