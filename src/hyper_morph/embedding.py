"""Utilities for validating disk mesh embeddings."""

from __future__ import annotations

import math
from typing import Dict, List, Mapping, Sequence, Tuple


Point = Sequence[float]
Edge = Tuple[int, int]
Face = Tuple[int, int, int]


def edges_from_faces(faces: Sequence[Sequence[int]]) -> List[Edge]:
    edges = set()
    for face in faces:
        if len(face) != 3:
            raise ValueError(f"Face must have three vertices, got {face!r}")
        i, j, k = (int(face[0]), int(face[1]), int(face[2]))
        for u, v in ((i, j), (j, k), (k, i)):
            if u != v:
                edges.add((min(u, v), max(u, v)))
    return sorted(edges)


def signed_area2(a: Point, b: Point, c: Point) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def point_norm2(point: Point) -> float:
    return point[0] * point[0] + point[1] * point[1]


def _segments_cross(a: Point, b: Point, c: Point, d: Point, eps: float) -> bool:
    def sign(value: float) -> int:
        return int(value > eps) - int(value < -eps)

    o1 = signed_area2(a, b, c)
    o2 = signed_area2(a, b, d)
    o3 = signed_area2(c, d, a)
    o4 = signed_area2(c, d, b)
    return sign(o1) * sign(o2) < 0 and sign(o3) * sign(o4) < 0


def _normalize_angle(delta: float) -> float:
    while delta <= -math.pi:
        delta += 2.0 * math.pi
    while delta > math.pi:
        delta -= 2.0 * math.pi
    return delta


def poincare_geodesic_points(
    a: Point,
    b: Point,
    samples: int = 36,
    eps: float = 1e-10,
) -> List[Tuple[float, float]]:
    det = a[0] * b[1] - a[1] * b[0]
    if abs(det) < eps:
        return [(float(a[0]), float(a[1])), (float(b[0]), float(b[1]))]

    rhs_a = 0.5 * (point_norm2(a) + 1.0)
    rhs_b = 0.5 * (point_norm2(b) + 1.0)
    center = (
        (rhs_a * b[1] - rhs_b * a[1]) / det,
        (a[0] * rhs_b - b[0] * rhs_a) / det,
    )
    radius = math.hypot(a[0] - center[0], a[1] - center[1])
    if not math.isfinite(radius) or radius <= eps:
        return [(float(a[0]), float(a[1])), (float(b[0]), float(b[1]))]

    start = math.atan2(a[1] - center[1], a[0] - center[0])
    end = math.atan2(b[1] - center[1], b[0] - center[0])
    delta = _normalize_angle(end - start)

    mid = (
        center[0] + radius * math.cos(start + 0.5 * delta),
        center[1] + radius * math.sin(start + 0.5 * delta),
    )
    if point_norm2(mid) > 1.0 + 1e-7:
        delta += -2.0 * math.pi if delta > 0.0 else 2.0 * math.pi

    return [
        (
            center[0] + radius * math.cos(start + delta * k / samples),
            center[1] + radius * math.sin(start + delta * k / samples),
        )
        for k in range(samples + 1)
    ]


def _bbox(points: Sequence[Point]) -> Tuple[float, float, float, float]:
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    return min(xs), min(ys), max(xs), max(ys)


def _bbox_overlaps(
    a: Tuple[float, float, float, float],
    b: Tuple[float, float, float, float],
) -> bool:
    return not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1])


def _polyline_crosses(
    first: Sequence[Point],
    second: Sequence[Point],
    eps: float,
) -> bool:
    for a, b in zip(first, first[1:]):
        for c, d in zip(second, second[1:]):
            if _segments_cross(a, b, c, d, eps):
                return True
    return False


def validate_disk_embedding(
    data: Mapping[str, object],
    *,
    check_geodesic_edges: bool = True,
    geodesic_samples: int = 36,
    eps: float = 1e-12,
    sample_limit: int = 5,
) -> Dict[str, object]:
    vertices_raw = data.get("vertices")
    if not isinstance(vertices_raw, list) or not vertices_raw:
        raise ValueError("Input JSON needs a non-empty 'vertices' list.")
    vertices = [
        (float(point[0]), float(point[1]))
        for point in vertices_raw
    ]

    faces_raw = data.get("faces", [])
    if faces_raw is None:
        faces_raw = []
    if not isinstance(faces_raw, list):
        raise ValueError("'faces' must be a list if provided.")
    faces: List[Face] = [
        (int(face[0]), int(face[1]), int(face[2]))
        for face in faces_raw
    ]

    edges_raw = data.get("edges")
    if edges_raw is None:
        edges = edges_from_faces(faces)
    else:
        if not isinstance(edges_raw, list):
            raise ValueError("'edges' must be a list if provided.")
        edges = [
            (min(int(edge[0]), int(edge[1])), max(int(edge[0]), int(edge[1])))
            for edge in edges_raw
            if int(edge[0]) != int(edge[1])
        ]

    outside_disk = [
        idx for idx, point in enumerate(vertices)
        if point_norm2(point) >= 1.0
    ]

    face_areas = [
        signed_area2(vertices[i], vertices[j], vertices[k])
        for i, j, k in faces
    ]
    positive_faces = sum(area > eps for area in face_areas)
    negative_faces = sum(area < -eps for area in face_areas)
    degenerate_faces = sum(abs(area) <= eps for area in face_areas)

    chord_crossings = 0
    chord_samples: List[Tuple[Edge, Edge]] = []
    for idx, (a, b) in enumerate(edges):
        for c, d in edges[idx + 1:]:
            if len({a, b, c, d}) < 4:
                continue
            if _segments_cross(vertices[a], vertices[b], vertices[c], vertices[d], eps):
                chord_crossings += 1
                if len(chord_samples) < sample_limit:
                    chord_samples.append(((a, b), (c, d)))

    geodesic_crossings = 0
    geodesic_samples_out: List[Tuple[Edge, Edge]] = []
    if check_geodesic_edges:
        polylines = [
            poincare_geodesic_points(
                vertices[a],
                vertices[b],
                samples=geodesic_samples,
            )
            for a, b in edges
        ]
        boxes = [_bbox(polyline) for polyline in polylines]
        for idx, (a, b) in enumerate(edges):
            for jdx, (c, d) in enumerate(edges[idx + 1:], idx + 1):
                if len({a, b, c, d}) < 4:
                    continue
                if not _bbox_overlaps(boxes[idx], boxes[jdx]):
                    continue
                if _polyline_crosses(polylines[idx], polylines[jdx], eps):
                    geodesic_crossings += 1
                    if len(geodesic_samples_out) < sample_limit:
                        geodesic_samples_out.append(((a, b), (c, d)))

    is_embedding = (
        not outside_disk
        and degenerate_faces == 0
        and (not faces or positive_faces == len(faces) or negative_faces == len(faces))
        and chord_crossings == 0
        and (not check_geodesic_edges or geodesic_crossings == 0)
    )

    return {
        "is_embedding": is_embedding,
        "vertex_count": len(vertices),
        "face_count": len(faces),
        "edge_count": len(edges),
        "outside_disk_count": len(outside_disk),
        "outside_disk_samples": outside_disk[:sample_limit],
        "positive_face_count": positive_faces,
        "negative_face_count": negative_faces,
        "degenerate_face_count": degenerate_faces,
        "min_abs_face_area2": min((abs(area) for area in face_areas), default=None),
        "euclidean_chord_crossing_count": chord_crossings,
        "euclidean_chord_crossing_samples": chord_samples,
        "poincare_geodesic_crossing_count": geodesic_crossings,
        "poincare_geodesic_crossing_samples": geodesic_samples_out,
    }


def require_disk_embedding(data: Mapping[str, object]) -> Dict[str, object]:
    report = validate_disk_embedding(data)
    if not report["is_embedding"]:
        raise ValueError(f"Input vertices/faces are not a valid embedding: {report}")
    return report
