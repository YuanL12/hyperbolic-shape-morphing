#!/usr/bin/env python3
"""Morph boundary vertices between directed edge-weight fields."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .embedding import validate_disk_embedding
from .solver import HarmonicMapSolver, MobiusIsometry, build_edges, from_complex, to_complex
from .weights import DirectedEdgeWeightCalculator

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
        "--source-embedding",
        default="examples/input/example_genus2_bolza_delaunay_all_1_solution.json",
        help=(
            "Source embedding JSON used for frame 0 vertices, graph topology, "
            "constraints, and fixed vertices."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="examples/output/frames/morph_frames",
        help="Directory for per-frame JSON outputs.",
    )
    parser.add_argument(
        "--unnormalized",
        action="store_true",
        help="Use unnormalized directed mean-value weights.",
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
        "--start-directed-weights",
        choices=("mean_value", "ones"),
        default="mean_value",
        help="Directed weights used at frame 0.",
    )
    parser.add_argument(
        "--target-directed-weights",
        choices=(
            "ones",
            "half_10",
            "half_oriented_10",
            "mean_value",
        ),
        default="ones",
        help=(
            "Directed weights used at the last frame. 'half_10' sets half the "
            "edges to (10, 10), half to (1, 1). 'half_oriented_10' sets half "
            "to (10, 1), half to (1, 10). 'mean_value' computes mean-value "
            "directed weights from --target-embedding."
        ),
    )
    parser.add_argument(
        "--target-embedding",
        help=(
            "Embedding JSON used when --target-directed-weights is "
            "'mean_value'. Topology is taken from --source-embedding."
        ),
    )
    parser.add_argument("--iterations", type=int, default=2200)
    parser.add_argument("--step-size", type=float, default=0.01)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument(
        "--line-search-objective",
        choices=("energy", "gradient_norm", "none"),
        default="none",
        help=(
            "Line-search objective passed to the solver. Defaults to 'none' "
            "because the directed mean-value residual is not the gradient of "
            "the averaged diagnostic energy."
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
        help="Skip validation that source vertices/faces form an embedding.",
    )
    parser.add_argument(
        "--skip-face-orientation-check",
        action="store_true",
        help="For non-triangle/coarsened maps, validate edge crossings only.",
    )
    parser.add_argument(
        "--reference-fundamental-domain",
        action="store_true",
        help=(
            "Copy the source reference_fundamental_domain into each frame. "
            "This keeps the reference fundamental polygon fixed in drawings "
            "when solver corner constraints have been relaxed."
        ),
    )
    args = parser.parse_args()
    args.normalization = "unnormalized" if args.unnormalized else "normalized"
    return args


def make_target_directed_weights(
    mode: str,
    edge_count: int,
) -> List[DirectedWeight]:
    if mode == "mean_value":
        raise ValueError(
            "mean_value requires --target-embedding and is handled "
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


def infer_reference_domain_name(corner_count: int) -> Tuple[str, str, str]:
    if corner_count == 8:
        return ("bolza_regular_octagon", "Bolza surface", "regular octagon")
    if corner_count == 14:
        return ("klein_quartic_regular_14_gon", "Klein quartic surface", "regular 14-gon")
    return ("reference_fundamental_polygon", "unknown surface", f"{corner_count}-gon")


def reference_fundamental_domain_from(data: Dict[str, object]) -> Dict[str, object]:
    existing = data.get("reference_fundamental_domain")
    if isinstance(existing, dict):
        return existing

    corner_indices = data.get("fixed", [])
    if not isinstance(corner_indices, list):
        corner_indices = []
    corner_indices = [int(idx) for idx in corner_indices]

    domain_vertices = [list(data["vertices"][idx]) for idx in corner_indices]

    name, surface, polygon = infer_reference_domain_name(len(corner_indices))
    return {
        "name": name,
        "surface": surface,
        "polygon": polygon,
        "model": "poincare_disk",
        "corner_indices": corner_indices,
        "vertices": domain_vertices,
        "description": (
            "Drawing-only reference fundamental polygon. These coordinates do "
            "not constrain the moving graph embedding."
        ),
    }


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
    boundary_vertex_count = metadata.get("boundary_vertex_count")
    if isinstance(side_labels, list):
        side_count = len(side_labels)
    elif boundary_vertex_count is not None:
        side_count = int(boundary_vertex_count) // side_steps
    else:
        reference_domain = data.get("reference_fundamental_domain")
        reference_indices = (
            reference_domain.get("corner_indices", [])
            if isinstance(reference_domain, dict)
            else []
        )
        side_count = len(
            reference_indices
            if isinstance(reference_indices, list) and reference_indices
            else data.get("fixed", [])
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

    source_embedding = load_json(args.source_embedding)
    if not args.skip_embedding_check:
        check_data = dict(source_embedding)
        if args.skip_face_orientation_check:
            check_data.pop("faces", None)
        embedding_report = validate_disk_embedding(check_data)
        if not embedding_report["is_embedding"]:
            raise ValueError(
                "Source embedding is not a valid disk embedding: "
                f"{embedding_report}"
            )
        print(
            "embedding check: "
            f"faces={embedding_report['face_count']} "
            f"edges={embedding_report['edge_count']} "
            "crossings=0"
        )

    initial_vertices = [list(v) for v in source_embedding["vertices"]]
    source_metadata = source_embedding.get("metadata", {})
    if not isinstance(source_metadata, dict):
        source_metadata = {}
    reference_fundamental_domain = reference_fundamental_domain_from(source_embedding)
    edges = [
        list(edge) for edge in build_edges(
            faces=source_embedding.get("faces"),
            edges=source_embedding.get("edges"),
            n_vertices=len(initial_vertices),
        )
    ]

    start_data = dict(source_embedding)
    start_data["vertices"] = initial_vertices
    start_data["edges"] = edges
    start_data["constraints"] = source_embedding.get("constraints", [])
    start_data["fixed"] = source_embedding.get("fixed", [])
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

    target_directed_weight_mode = args.target_directed_weights
    target_vertices: Optional[List[List[float]]] = None

    if target_directed_weight_mode == "mean_value":
        if not args.target_embedding:
            raise ValueError(
                "--target-embedding is required when --target-directed-weights "
                "is mean_value."
            )
        target_embedding = load_json(args.target_embedding)
        if len(target_embedding.get("vertices", [])) != len(initial_vertices):
            raise ValueError(
                "--target-embedding must have the same number of vertices as "
                "--source-embedding."
            )
        target_vertices = [list(v) for v in target_embedding["vertices"]]
        target_data = dict(source_embedding)
        target_data["vertices"] = target_vertices
        target_data["edges"] = edges
        target_data["constraints"] = source_embedding.get("constraints", [])
        target_data["fixed"] = source_embedding.get("fixed", [])
        weight_target_data = with_corner_attachment_constraints(target_data)
        target_calculator = DirectedEdgeWeightCalculator(
            weight_target_data,
            normalization=args.normalization,
            low_valence_policy=args.low_valence_star_policy,
        )
        target_directed_weights = target_calculator.compute()
    else:
        target_directed_weights = make_target_directed_weights(
            target_directed_weight_mode,
            len(start_directed_weights),
        )

    # create the output directory
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prev_vertices = initial_vertices

    for frame_id in range(frame_count):
        t = frame_id / (frame_count - 1)
        current_directed_weights = interpolate_directed_weights(
            start_directed_weights,
            target_directed_weights,
            t,
        )
        current_edge_weights = average_edge_weights(current_directed_weights)

        is_source_endpoint = (
            frame_id == 0 and args.start_directed_weights == "mean_value"
        )
        is_target_endpoint = (
            frame_id == frame_count - 1
            and target_directed_weight_mode == "mean_value"
            and target_vertices is not None
        )
        if is_source_endpoint:
            seed_vertices = initial_vertices
            endpoint_passthrough = "source"
        elif is_target_endpoint:
            seed_vertices = target_vertices
            endpoint_passthrough = "target"
        else:
            seed_vertices = prev_vertices
            endpoint_passthrough = None

        data = dict(source_embedding)
        data["vertices"] = seed_vertices
        data["edges"] = edges
        data["edge_weights"] = current_edge_weights
        data["directed_edge_weights"] = [
            [weight_i, weight_j] for weight_i, weight_j in current_directed_weights
        ]
        data["constraints"] = source_embedding.get("constraints", [])
        data["fixed"] = source_embedding.get("fixed", [])
        data["iterations"] = 0 if endpoint_passthrough else args.iterations
        data["step_size"] = args.step_size
        data["tolerance"] = args.tolerance
        data["edge_force_model"] = "hyperbolic_mean_value"
        if args.line_search_objective is not None:
            data["line_search_objective"] = args.line_search_objective
        if args.convergence_criterion is not None:
            data["convergence_criterion"] = args.convergence_criterion

        solver = HarmonicMapSolver(data)
        result = solver.solve()
        result["morph_t"] = t
        result["edge_weights"] = current_edge_weights
        result["directed_edge_weights"] = data["directed_edge_weights"]
        if "faces" in source_embedding:
            result["faces"] = source_embedding["faces"]
        result["constraints"] = source_embedding.get("constraints", [])
        result["fixed"] = source_embedding.get("fixed", [])
        if args.reference_fundamental_domain:
            result["reference_fundamental_domain"] = reference_fundamental_domain
        result["metadata"] = {
            "morph_direction": (
                f"{args.start_directed_weights}_directed_to_"
                f"{target_directed_weight_mode}"
            ),
            "source_embedding": args.source_embedding,
            "start_weight_normalization": args.normalization,
            "low_valence_star_policy": args.low_valence_star_policy,
            "attach_corner_stars": True,
            "reference_fundamental_domain": args.reference_fundamental_domain,
            "relaxed_corner_vertex": source_metadata.get("relaxed_corner_vertex"),
            "relaxed_corner_position": source_metadata.get("relaxed_corner_position"),
            "relaxed_corner_vertices": source_metadata.get("relaxed_corner_vertices"),
            "relaxed_corner_positions": source_metadata.get("relaxed_corner_positions"),
            "relaxed_corner_component_vertices": source_metadata.get(
                "relaxed_corner_component_vertices"
            ),
            "relaxed_corner_constraint_count": source_metadata.get(
                "relaxed_corner_constraint_count"
            ),
            "start_directed_weights": args.start_directed_weights,
            "target_directed_weights": target_directed_weight_mode,
            "target_embedding": args.target_embedding,
            "endpoint_passthrough": endpoint_passthrough,
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
            f"stop={result['stop_reason']} "
            f"max_directed_w={result['metadata']['current_max_directed_weight']:.3g}"
        )

    print(f"Wrote {out_dir}")


if __name__ == "__main__":
    main()
