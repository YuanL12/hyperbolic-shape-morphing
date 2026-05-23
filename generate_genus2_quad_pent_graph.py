#!/usr/bin/env python3
"""Remove Delaunay edges from the genus-2 octagon to make quad/pent cells."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from embedding_utils import validate_disk_embedding
from generate_genus3_klein_embedding import write_preview_svg
from generate_genus3_quad_pent_graph import (
    Edge,
    coarsen_by_removing_edges,
    edge_key,
)


def quotient_parent(
    vertex_count: int,
    pairings: Sequence[Dict[str, int]],
    side_subdivisions: int,
    side_count: int,
) -> List[int]:
    side_steps = side_subdivisions + 1
    boundary_count = side_count * side_steps
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

    for pairing in pairings:
        master_side = int(pairing["master_side"])
        slave_side = int(pairing["slave_side"])
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
    pairings: Sequence[Dict[str, int]],
    side_subdivisions: int,
    side_count: int,
) -> Tuple[int, List[Tuple[int, int]], int]:
    parent = quotient_parent(vertex_count, pairings, side_subdivisions, side_count)
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

    cut_vertices = [vertex for vertex in sorted(vertices) if not connected([vertex])]
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


def geodesic_tangent(
    vertices: Sequence[Sequence[float]],
    a: int,
    b: int,
) -> Tuple[float, float]:
    ax, ay = vertices[a]
    bx, by = vertices[b]
    det = ax * by - ay * bx
    if abs(det) < 1e-12:
        vx = bx - ax
        vy = by - ay
    else:
        rhs_a = 0.5 * (ax * ax + ay * ay + 1.0)
        rhs_b = 0.5 * (bx * bx + by * by + 1.0)
        cx = (rhs_a * by - rhs_b * ay) / det
        cy = (ax * rhs_b - bx * rhs_a) / det
        vx = -(ay - cy)
        vy = ax - cx
        eps = 1e-6
        plus = (ax + eps * vx - bx) ** 2 + (ay + eps * vy - by) ** 2
        minus = (ax - eps * vx - bx) ** 2 + (ay - eps * vy - by) ** 2
        if minus < plus:
            vx = -vx
            vy = -vy
    norm = math.hypot(vx, vy)
    return vx / norm, vy / norm


def signed_angle(u: Tuple[float, float], v: Tuple[float, float]) -> float:
    return math.atan2(u[0] * v[1] - u[1] * v[0], u[0] * v[0] + u[1] * v[1])


def face_area(face: Sequence[int], vertices: Sequence[Sequence[float]]) -> float:
    return 0.5 * sum(
        vertices[face[k]][0] * vertices[face[(k + 1) % len(face)]][1]
        - vertices[face[k]][1] * vertices[face[(k + 1) % len(face)]][0]
        for k in range(len(face))
    )


def nonconvex_hyperbolic_faces(data: Dict[str, object]) -> List[Tuple[int, int]]:
    vertices = data["vertices"]
    bad = []
    for face_id, face in enumerate(data.get("faces", [])):
        if len(face) <= 3:
            continue
        orientation = 1 if face_area(face, vertices) > 0.0 else -1
        ok = True
        for k, vertex in enumerate(face):
            prev_vertex = face[k - 1]
            next_vertex = face[(k + 1) % len(face)]
            to_prev = geodesic_tangent(vertices, vertex, prev_vertex)
            incoming = (-to_prev[0], -to_prev[1])
            outgoing = geodesic_tangent(vertices, vertex, next_vertex)
            turn = signed_angle(incoming, outgoing)
            if orientation > 0 and turn <= 1e-8:
                ok = False
                break
            if orientation < 0 and turn >= -1e-8:
                ok = False
                break
        if not ok:
            bad.append((face_id, len(face)))
    return bad


def find_coarsening(
    base: Dict[str, object],
    *,
    quad_target: int,
    pentagon_target: int,
    seed_start: int,
    seed_limit: int,
) -> Tuple[List[List[int]], List[List[int]], List[Edge], Dict[str, object], int]:
    metadata = base["metadata"]
    pairings = metadata["pairings"]
    side_subdivisions = int(metadata["side_subdivisions"])
    side_count = len(base["fixed"])
    for seed in range(seed_start, seed_limit):
        faces, edges, removed_edges = coarsen_by_removing_edges(
            base,
            quad_target=quad_target,
            pentagon_target=pentagon_target,
            seed=seed,
        )
        trial = dict(base)
        trial["faces"] = faces
        trial["edges"] = edges
        if nonconvex_hyperbolic_faces(trial):
            continue
        quotient_vertices, quotient_edges, loop_count = quotient_graph(
            len(base["vertices"]),
            edges,
            pairings,
            side_subdivisions,
            side_count,
        )
        report = connectivity_report(quotient_vertices, quotient_edges)
        if report["passes_finite_3_connected_quotient_check"]:
            report["quotient_loop_edges_ignored"] = loop_count
            return faces, edges, removed_edges, report, seed
    raise RuntimeError("Could not find a convex 3-connected coarsening.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        default="input/example_genus2_bolza_delaunay_all_1_solution.json",
    )
    parser.add_argument(
        "--output",
        default="input/example_genus2_bolza_convex_quad_pent_graph.json",
    )
    parser.add_argument(
        "--preview-svg",
        default="output/genus2_bolza_convex_quad_pent_graph.svg",
    )
    parser.add_argument("--quads", type=int, default=55)
    parser.add_argument("--pentagons", type=int, default=14)
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--seed-limit", type=int, default=400)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with open(args.base, "r", encoding="utf-8") as fh:
        base = json.load(fh)

    faces, edges, removed_edges, report, seed = find_coarsening(
        base,
        quad_target=args.quads,
        pentagon_target=args.pentagons,
        seed_start=args.seed_start,
        seed_limit=args.seed_limit,
    )
    graph = {
        "name": "genus2_bolza_octagon_convex_quad_pent_3_connected_graph",
        "description": (
            "Planar embedded graph generated by removing selected Delaunay "
            "edges from the irregular genus-2 Bolza octagon embedding. It "
            "contains quadrilateral and pentagonal cells and passes the finite "
            "quotient 3-connectedness check."
        ),
        "vertices": base["vertices"],
        "edges": edges,
        "faces": faces,
        "fixed": base.get("fixed", []),
        "constraints": base.get("constraints", []),
        "metadata": {
            **base.get("metadata", {}),
            "source": args.base,
            "graph_construction": "edge removal from irregular embedded Delaunay triangulation",
            "removed_construction_edges": len(removed_edges),
            "removed_edge_samples": removed_edges[:30],
            "face_size_histogram": sorted(Counter(len(face) for face in faces).items()),
            "quad_target": args.quads,
            "pentagon_target": args.pentagons,
            "coarsening_random_seed": seed,
            "nonconvex_hyperbolic_face_count": 0,
            **report,
        },
    }

    edge_only = dict(graph)
    edge_only.pop("faces", None)
    embedding_report = validate_disk_embedding(edge_only)
    graph["metadata"]["edge_embedding_validation"] = embedding_report
    if not embedding_report["is_embedding"]:
        raise RuntimeError(f"Graph edges are not embedded: {embedding_report}")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fh:
        json.dump(graph, fh, indent=2)
        fh.write("\n")

    preview = Path(args.preview_svg)
    write_preview_svg(graph, preview)
    print(
        f"wrote {output}: vertices={len(graph['vertices'])} "
        f"edges={len(graph['edges'])} faces={len(graph['faces'])}"
    )
    print(f"preview: {preview}")
    print(f"face sizes: {graph['metadata']['face_size_histogram']}")
    print(
        "quotient check: "
        f"vertices={report['quotient_vertex_count']} "
        f"edges={report['quotient_edge_count']} "
        f"min_degree={report['min_quotient_degree']} "
        f"passes={report['passes_finite_3_connected_quotient_check']}"
    )


if __name__ == "__main__":
    main()
