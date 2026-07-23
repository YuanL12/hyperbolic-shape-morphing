#!/usr/bin/env python3
"""Generate a minimal Poincaré-disk mesh for the solver example."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("examples/input/simple_disk.json"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = {
        "vertices": [
            [0.18, 0.12],
            [-0.45, -0.45],
            [0.45, -0.45],
            [0.45, 0.45],
            [-0.45, 0.45],
        ],
        "faces": [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]],
        "fixed": [1, 2, 3, 4],
        "iterations": 300,
        "step_size": 0.05,
        "tolerance": 1e-10,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
