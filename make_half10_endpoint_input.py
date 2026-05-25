#!/usr/bin/env python3
"""Create a half-10 weighted endpoint input from an embedding JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict


def load_json(path: str) -> Dict[str, object]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError("Top-level JSON must be an object.")
    return data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", help="Source graph or embedding JSON.")
    parser.add_argument("--output", required=True, help="Output weighted JSON.")
    parser.add_argument(
        "--heavy-weight",
        type=float,
        default=10.0,
        help="Weight assigned to the first half of edges.",
    )
    parser.add_argument(
        "--light-weight",
        type=float,
        default=1.0,
        help="Weight assigned to the second half of edges.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = load_json(args.source)
    edges = data.get("edges")
    if not isinstance(edges, list):
        raise ValueError("Input JSON must contain an 'edges' list.")

    edge_count = len(edges)
    half_edge_count = edge_count // 2
    weights = [
        args.heavy_weight if edge_idx < half_edge_count else args.light_weight
        for edge_idx in range(edge_count)
    ]
    data["edge_weights"] = weights
    data["directed_edge_weights"] = [[weight, weight] for weight in weights]
    data.setdefault("metadata", {})
    if isinstance(data["metadata"], dict):
        data["metadata"]["endpoint_weight_pattern"] = (
            f"first half {args.heavy_weight:g}, second half {args.light_weight:g}"
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
