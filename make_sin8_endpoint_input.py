#!/usr/bin/env python3
"""Create a smooth angular sin(8 theta) weighted endpoint input."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Sequence


def load_json(path: str) -> Dict[str, object]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError("Top-level JSON must be an object.")
    return data


def midpoint_angle(a: Sequence[float], b: Sequence[float]) -> float:
    x = 0.5 * (float(a[0]) + float(b[0]))
    y = 0.5 * (float(a[1]) + float(b[1]))
    return math.atan2(y, x)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", help="Source graph or embedding JSON.")
    parser.add_argument("--output", required=True, help="Output weighted JSON.")
    parser.add_argument("--min-weight", type=float, default=1.0)
    parser.add_argument("--max-weight", type=float, default=50.0)
    parser.add_argument("--frequency", type=float, default=8.0)
    parser.add_argument(
        "--phase",
        type=float,
        default=0.0,
        help="Angular phase shift in radians inside sin(frequency * theta + phase).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.min_weight <= 0.0 or args.max_weight <= 0.0:
        raise ValueError("Weights must be positive.")
    if args.max_weight < args.min_weight:
        raise ValueError("--max-weight must be at least --min-weight.")

    data = load_json(args.source)
    vertices = data.get("vertices")
    edges = data.get("edges")
    if not isinstance(vertices, list) or not isinstance(edges, list):
        raise ValueError("Input JSON must contain 'vertices' and 'edges' lists.")

    weights = []
    for edge in edges:
        i, j = int(edge[0]), int(edge[1])
        theta = midpoint_angle(vertices[i], vertices[j])
        unit = 0.5 + 0.5 * math.sin(args.frequency * theta + args.phase)
        weights.append(args.min_weight + (args.max_weight - args.min_weight) * unit)

    data["edge_weights"] = weights
    data["directed_edge_weights"] = [[weight, weight] for weight in weights]
    data.setdefault("metadata", {})
    if isinstance(data["metadata"], dict):
        data["metadata"]["endpoint_weight_pattern"] = (
            f"smooth sin({args.frequency:g} theta + {args.phase:g}) "
            f"weights from {args.min_weight:g} to {args.max_weight:g}"
        )
        data["metadata"]["endpoint_weight_min"] = min(weights)
        data["metadata"]["endpoint_weight_max"] = max(weights)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    print(
        f"Wrote {output}: min_weight={min(weights):.6g} "
        f"max_weight={max(weights):.6g}"
    )


if __name__ == "__main__":
    main()
