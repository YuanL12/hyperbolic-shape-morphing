#!/usr/bin/env python3
"""Extract and preview an essentially 3-connected graph on the genus-3 surface."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from generate_genus3_klein_embedding import (
    KLEIN_SIDE_LABELS,
    write_preview_svg,
)


def quotient_parent(vertex_count: int, side_subdivisions: int) -> List[int]:
    side_steps = side_subdivisions + 1
    boundary_count = len(KLEIN_SIDE_LABELS) * side_steps
    parent = list(range(vertex_count))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        root_a = find(a)
        root_b = find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    for label in sorted(set(KLEIN_SIDE_LABELS)):
        master_side, slave_side = [
            side for side, side_label in enumerate(KLEIN_SIDE_LABELS)
            if side_label == label
        ]
        for offset in range(1, side_steps):
            union(
                master_side * side_steps + offset,
                slave_side * side_steps + (side_steps - offset),
            )
        union(
            master_side * side_steps,
            ((slave_side + 1) * side_steps) % boundary_count,
        )
        union(
            ((master_side + 1) * side_steps) % boundary_count,
            slave_side * side_steps,
        )

    for vertex in range(vertex_count):
        find(vertex)
    return parent


def quotient_graph(
    vertex_count: int,
    edges: Sequence[Sequence[int]],
    side_subdivisions: int,
) -> Tuple[int, List[Tuple[int, int]], int]:
    parent = quotient_parent(vertex_count, side_subdivisions)
    quotient_index: Dict[int, int] = {}

    def quotient_vertex(vertex: int) -> int:
        root = parent[vertex]
        if root not in quotient_index:
            quotient_index[root] = len(quotient_index)
        return quotient_index[root]

    quotient_edges = set()
    loop_count = 0
    for i, j in edges:
        qi = quotient_vertex(int(i))
        qj = quotient_vertex(int(j))
        if qi == qj:
            loop_count += 1
            continue
        quotient_edges.add((min(qi, qj), max(qi, qj)))

    return len(quotient_index), sorted(quotient_edges), loop_count


def connectivity_report(vertex_count: int, edges: Sequence[Tuple[int, int]]) -> Dict[str, object]:
    adjacency: Dict[int, set[int]] = defaultdict(set)
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)

    vertices = set(range(vertex_count))

    def connected(removed: Sequence[int] = ()) -> bool:
        removed_set = set(removed)
        starts = [vertex for vertex in vertices if vertex not in removed_set]
        if not starts:
            return True
        seen = {starts[0]}
        queue = deque([starts[0]])
        while queue:
            vertex = queue.popleft()
            for neighbor in adjacency[vertex]:
                if neighbor not in removed_set and neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        return len(seen) == len(starts)

    cut_vertices = [
        vertex for vertex in sorted(vertices)
        if not connected([vertex])
    ]
    cut_pairs = []
    for first, second in itertools.combinations(sorted(vertices), 2):
        if not connected([first, second]):
            cut_pairs.append([first, second])
            if len(cut_pairs) >= 20:
                break

    degrees = [len(adjacency[vertex]) for vertex in sorted(vertices)]
    return {
        "quotient_vertex_count": vertex_count,
        "quotient_edge_count": len(edges),
        "min_quotient_degree": min(degrees) if degrees else 0,
        "max_quotient_degree": max(degrees) if degrees else 0,
        "quotient_degree_histogram": sorted(Counter(degrees).items()),
        "is_connected": connected(),
        "cut_vertex_count": len(cut_vertices),
        "cut_vertex_samples": cut_vertices[:20],
        "has_two_vertex_cut": bool(cut_pairs),
        "two_vertex_cut_samples": cut_pairs,
        "passes_finite_3_connected_quotient_check": (
            connected()
            and not cut_vertices
            and not cut_pairs
            and min(degrees, default=0) >= 3
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        default="examples/input/example_genus3_klein_quartic_embedded_irregular.json",
    )
    parser.add_argument(
        "--output",
        default="examples/input/example_genus3_klein_quartic_essential_graph.json",
    )
    parser.add_argument(
        "--preview-svg",
        default="examples/output/genus3_klein_essential_graph.svg",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with open(args.base, "r", encoding="utf-8") as fh:
        base = json.load(fh)

    side_subdivisions = int(base["metadata"]["side_subdivisions"])
    quotient_vertices, quotient_edges, loop_count = quotient_graph(
        len(base["vertices"]),
        base["edges"],
        side_subdivisions,
    )
    report = connectivity_report(quotient_vertices, quotient_edges)

    graph = {
        "name": "genus3_klein_quartic_essentially_3_connected_graph",
        "description": (
            "Graph-only 1-skeleton extracted from the irregular genus-3 "
            "Klein-quartic-style embedding. Faces are intentionally omitted."
        ),
        "vertices": base["vertices"],
        "edges": base["edges"],
        "fixed": base.get("fixed", []),
        "constraints": base.get("constraints", []),
        "metadata": {
            "source": args.base,
            "surface": "Klein quartic",
            "genus": 3,
            "fundamental_domain": "regular 14-gon",
            "side_labels": KLEIN_SIDE_LABELS,
            "graph_construction": "irregular embedded triangulation 1-skeleton",
            "faces_omitted": True,
            "quotient_loop_edges_ignored": loop_count,
            **report,
        },
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        json.dump(graph, fh, indent=2)
        fh.write("\n")

    preview = Path(args.preview_svg)
    write_preview_svg(graph, preview)

    print(
        f"wrote {output}: vertices={len(graph['vertices'])} "
        f"edges={len(graph['edges'])}"
    )
    print(f"preview: {preview}")
    print(
        "quotient check: "
        f"vertices={report['quotient_vertex_count']} "
        f"edges={report['quotient_edge_count']} "
        f"min_degree={report['min_quotient_degree']} "
        f"max_degree={report['max_quotient_degree']} "
        f"cut_vertices={report['cut_vertex_count']} "
        f"has_two_vertex_cut={report['has_two_vertex_cut']} "
        f"passes={report['passes_finite_3_connected_quotient_check']}"
    )


if __name__ == "__main__":
    main()
