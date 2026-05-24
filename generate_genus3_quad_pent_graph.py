#!/usr/bin/env python3
"""Remove construction edges to create a 3-connected graph with quad/pent faces."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from embedding_utils import validate_disk_embedding
from generate_genus3_essential_graph import connectivity_report, quotient_graph
from generate_genus3_klein_embedding import KLEIN_SIDE_LABELS, geodesic_points


Edge = Tuple[int, int]


def edge_key(i: int, j: int) -> Edge:
    return (min(i, j), max(i, j))


def face_edges(face: Sequence[int]) -> List[Edge]:
    return [
        edge_key(int(face[k]), int(face[(k + 1) % len(face)]))
        for k in range(len(face))
    ]


def face_boundary(cluster: Sequence[int], triangles: Sequence[Sequence[int]]) -> Optional[List[int]]:
    edge_counts: Dict[Edge, int] = defaultdict(int)
    for face_id in cluster:
        for edge in face_edges(triangles[face_id]):
            edge_counts[edge] += 1

    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    adjacency: Dict[int, List[int]] = defaultdict(list)
    for i, j in boundary_edges:
        adjacency[i].append(j)
        adjacency[j].append(i)

    if not adjacency or any(len(neighbors) != 2 for neighbors in adjacency.values()):
        return None

    start = min(adjacency)
    cycle = [start]
    prev = None
    current = start
    while True:
        candidates = [vertex for vertex in adjacency[current] if vertex != prev]
        if not candidates:
            return None
        nxt = candidates[0]
        if nxt == start:
            break
        if nxt in cycle:
            return None
        cycle.append(nxt)
        prev, current = current, nxt

    if len(cycle) != len(boundary_edges):
        return None
    return cycle


def build_dual_adjacency(triangles: Sequence[Sequence[int]]) -> Dict[int, List[int]]:
    edge_to_faces: Dict[Edge, List[int]] = defaultdict(list)
    for face_id, face in enumerate(triangles):
        for edge in face_edges(face):
            edge_to_faces[edge].append(face_id)

    dual: Dict[int, List[int]] = defaultdict(list)
    for incident_faces in edge_to_faces.values():
        if len(incident_faces) != 2:
            continue
        first, second = incident_faces
        dual[first].append(second)
        dual[second].append(first)
    return dual


def coarsen_by_removing_edges(
    base: Dict[str, object],
    *,
    quad_target: int,
    pentagon_target: int,
    seed: int,
) -> Tuple[List[List[int]], List[List[int]], List[Edge]]:
    triangles = base["faces"]
    dual = build_dual_adjacency(triangles)
    rng = random.Random(seed)
    unassigned = set(range(len(triangles)))
    clusters: List[List[int]] = []

    for _ in range(pentagon_target):
        options = []
        for center in list(unassigned):
            neighbors = [neighbor for neighbor in dual[center] if neighbor in unassigned]
            for first, second in itertools.combinations(neighbors, 2):
                cluster = [center, first, second]
                boundary = face_boundary(cluster, triangles)
                if boundary is not None and len(boundary) == 5:
                    options.append(cluster)
        if not options:
            break
        cluster = rng.choice(options)
        clusters.append(cluster)
        unassigned.difference_update(cluster)

    quad_count = 0
    while quad_count < quad_target:
        options = []
        faces = list(unassigned)
        rng.shuffle(faces)
        for face_id in faces:
            for neighbor in dual[face_id]:
                if neighbor not in unassigned:
                    continue
                cluster = [face_id, neighbor]
                boundary = face_boundary(cluster, triangles)
                if boundary is not None and len(boundary) == 4:
                    options.append(cluster)
            if len(options) > 50:
                break
        if not options:
            break
        cluster = rng.choice(options)
        clusters.append(cluster)
        unassigned.difference_update(cluster)
        quad_count += 1

    clusters.extend([[face_id] for face_id in sorted(unassigned)])

    coarse_faces = []
    removed_edges = set()
    for cluster in clusters:
        boundary = face_boundary(cluster, triangles)
        if boundary is None:
            raise RuntimeError(f"Invalid merged face cluster: {cluster}")
        coarse_faces.append(boundary)

        edge_counts: Dict[Edge, int] = defaultdict(int)
        for face_id in cluster:
            for edge in face_edges(triangles[face_id]):
                edge_counts[edge] += 1
        removed_edges.update(edge for edge, count in edge_counts.items() if count == 2)

    base_edges = {edge_key(int(i), int(j)) for i, j in base["edges"]}
    coarse_edges = sorted(base_edges - removed_edges)
    return coarse_faces, [[i, j] for i, j in coarse_edges], sorted(removed_edges)


def find_3_connected_coarsening(
    base: Dict[str, object],
    *,
    quad_target: int,
    pentagon_target: int,
    seed_start: int,
    seed_limit: int,
) -> Tuple[List[List[int]], List[List[int]], List[Edge], Dict[str, object], int]:
    side_subdivisions = int(base["metadata"]["side_subdivisions"])
    for seed in range(seed_start, seed_limit):
        faces, edges, removed_edges = coarsen_by_removing_edges(
            base,
            quad_target=quad_target,
            pentagon_target=pentagon_target,
            seed=seed,
        )
        quotient_vertices, quotient_edges, loop_count = quotient_graph(
            len(base["vertices"]),
            edges,
            side_subdivisions,
        )
        report = connectivity_report(quotient_vertices, quotient_edges)
        if report["passes_finite_3_connected_quotient_check"]:
            report["quotient_loop_edges_ignored"] = loop_count
            return faces, edges, removed_edges, report, seed
    raise RuntimeError(
        "Could not find a 3-connected quotient coarsening in the requested seed range."
    )


def write_preview_svg(data: Dict[str, object], output: Path) -> None:
    vertices = data["vertices"]
    edges = data["edges"]
    faces = data["faces"]
    fixed = data.get("fixed", [])
    size = 1000
    margin = 64
    scale = (size - 2 * margin) / 2.0

    def project(point: Sequence[float]) -> Tuple[float, float]:
        return (
            size / 2.0 + point[0] * scale,
            size / 2.0 - point[1] * scale,
        )

    def path_for(points: Sequence[Sequence[float]]) -> str:
        projected = [project(point) for point in points]
        first = projected[0]
        rest = " ".join(f"L {x:.3f} {y:.3f}" for x, y in projected[1:])
        return f"M {first[0]:.3f} {first[1]:.3f} {rest}"

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="1000" viewBox="0 0 1000 1000">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<circle cx="{size / 2}" cy="{size / 2}" r="{scale}" fill="#fbfcfd" stroke="#111827" stroke-width="2"/>',
    ]

    for face in faces:
        if len(face) == 3:
            fill = "#e8eef8"
        elif len(face) == 4:
            fill = "#f3d8a8"
        elif len(face) == 5:
            fill = "#d7e8c7"
        else:
            fill = "#eadcf4"
        polygon = " ".join(
            f"{x:.3f},{y:.3f}" for x, y in (project(vertices[idx]) for idx in face)
        )
        lines.append(
            f'<polygon points="{polygon}" fill="{fill}" fill-opacity="0.52" stroke="none"/>'
        )

    for i, j in edges:
        d = path_for(geodesic_points(vertices[i], vertices[j], samples=20))
        lines.append(
            f'<path d="{d}" fill="none" stroke="#166534" stroke-width="1.65" stroke-linecap="round"/>'
        )

    for k, i in enumerate(fixed):
        j = fixed[(k + 1) % len(fixed)]
        d = path_for(geodesic_points(vertices[i], vertices[j], samples=36))
        lines.append(
            f'<path d="{d}" fill="none" stroke="#000000" stroke-width="2.2" stroke-dasharray="8 7" stroke-linecap="round"/>'
        )

    fixed_set = set(fixed)
    for idx, point in enumerate(vertices):
        x, y = project(point)
        if idx in fixed_set:
            lines.append(
                f'<circle cx="{x:.3f}" cy="{y:.3f}" r="4.3" fill="#dc2626" stroke="#7f1d1d" stroke-width="1"/>'
            )
        else:
            lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="2.1" fill="#111827"/>')

    lines.append("</svg>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        default="input/example_genus3_klein_quartic_embedded_irregular.json",
    )
    parser.add_argument(
        "--output",
        default="input/example_genus3_klein_quartic_convex_quad_pent_graph.json",
    )
    parser.add_argument(
        "--preview-svg",
        default="output/genus3_klein_quad_pent_graph.svg",
    )
    parser.add_argument("--quads", type=int, default=80)
    parser.add_argument("--pentagons", type=int, default=20)
    parser.add_argument("--seed-start", type=int, default=38)
    parser.add_argument("--seed-limit", type=int, default=39)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with open(args.base, "r", encoding="utf-8") as fh:
        base = json.load(fh)

    faces, edges, removed_edges, report, seed = find_3_connected_coarsening(
        base,
        quad_target=args.quads,
        pentagon_target=args.pentagons,
        seed_start=args.seed_start,
        seed_limit=args.seed_limit,
    )
    face_histogram = sorted(Counter(len(face) for face in faces).items())

    graph = {
        "name": "genus3_klein_quartic_quad_pent_3_connected_graph",
        "description": (
            "Planar embedded graph generated by removing selected construction "
            "edges from the irregular genus-3 Klein-quartic embedding. The graph "
            "includes quadrilateral and pentagonal faces and passes the finite "
            "quotient 3-connectedness check."
        ),
        "vertices": base["vertices"],
        "edges": edges,
        "faces": faces,
        "fixed": base.get("fixed", []),
        "constraints": base.get("constraints", []),
        "metadata": {
            "source": args.base,
            "surface": "Klein quartic",
            "genus": 3,
            "fundamental_domain": "regular 14-gon",
            "side_labels": KLEIN_SIDE_LABELS,
            "graph_construction": "edge removal from irregular embedded triangulation",
            "removed_construction_edges": len(removed_edges),
            "removed_edge_samples": removed_edges[:30],
            "face_size_histogram": face_histogram,
            "quad_target": args.quads,
            "pentagon_target": args.pentagons,
            "random_seed": seed,
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
    print(f"face sizes: {face_histogram}")
    print(
        "quotient check: "
        f"vertices={report['quotient_vertex_count']} "
        f"edges={report['quotient_edge_count']} "
        f"min_degree={report['min_quotient_degree']} "
        f"cut_vertices={report['cut_vertex_count']} "
        f"has_two_vertex_cut={report['has_two_vertex_cut']} "
        f"passes={report['passes_finite_3_connected_quotient_check']}"
    )


if __name__ == "__main__":
    main()
