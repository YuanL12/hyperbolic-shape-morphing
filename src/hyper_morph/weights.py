#!/usr/bin/env python3
"""Directed hyperbolic mean-value edge weights.

The current solver stores one scalar for each undirected edge.  This module
computes the two endpoint-centered weights needed for a directed discretization:
for edge ``(i, j)``, one weight from the star around ``i`` and one from the star
around ``j``.

When a cut mesh has slave vertices constrained to master vertices by Mobius
isometries, the star around an endpoint is assembled from every identified copy
of that endpoint.  Incident neighbors from another copy are transported into the
endpoint's local disk chart before the cyclic mean-value formula is evaluated.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .solver import (
    EPS,
    MobiusIsometry,
    build_edges,
    ensure_in_disk,
    from_complex,
    to_complex,
)


DirectedWeight = Tuple[float, float]
Edge = Tuple[int, int]
PointInput = Sequence[float] | complex
NormalizationMode = str
ANGLE_TOL = 1e-12


def compose_mobius(first: MobiusIsometry, second: MobiusIsometry) -> MobiusIsometry:
    """Return ``first after second`` for SU(1,1) disk isometries."""

    return MobiusIsometry(
        first.a * second.a + first.b * second.b.conjugate(),
        first.a * second.b + first.b * second.a.conjugate(),
    )


def star_angle_between(angle_a: float, angle_b: float) -> float:
    """Return the smaller unoriented angle between two geodesic directions."""

    angle = (angle_b - angle_a) % (2.0 * math.pi)
    if angle > math.pi:
        angle = 2.0 * math.pi - angle
    if abs(angle) <= ANGLE_TOL:
        return 0.0
    return angle


def positive_tan_half_angle(angle: float) -> float:
    """Return tan(angle / 2) with stable positive handling near pi."""

    if not (0.0 <= angle <= math.pi + ANGLE_TOL):
        raise ValueError(
            "Mean-value star angle must be in [0, pi]; "
            f"got {angle:.16g}."
        )
    if angle <= ANGLE_TOL:
        return 0.0
    if angle >= math.pi - ANGLE_TOL:
        return 1.0 / math.tan(0.5 * ANGLE_TOL)
    return math.tan(0.5 * angle)


IDENTITY_ISOMETRY = MobiusIsometry(1.0 + 0.0j, 0.0j)


@dataclass(frozen=True)
class StarOccurrence:
    """One neighbor occurrence in an attached vertex star."""

    source_center: int
    neighbor: int
    neighbor_position: complex
    edge: Edge


@dataclass(frozen=True)
class Constraint:
    slave: int
    master: int
    isometry: MobiusIsometry


class DirectedEdgeWeightCalculator:
    """Compute endpoint-centered hyperbolic mean-value edge weights.

    Parameters
    ----------
    data:
        Mesh JSON using the same schema as :mod:`hyper_morph.solver`.
    positions:
        Optional vertex positions.  If omitted, ``data["vertices"]`` is used.
        Full vertex positions and root-only positions are both accepted.  Slave
        vertices are resolved from their masters before weights are computed.
    """

    def __init__(
        self,
        data: Mapping[str, object],
        positions: Optional[Sequence[PointInput]] = None,
        normalization: NormalizationMode = "unnormalized",
        low_valence_policy: str = "error",
    ) -> None:
        if normalization not in {"unnormalized", "normalized"}:
            raise ValueError(
                "normalization must be either 'unnormalized' or 'normalized'."
            )
        if low_valence_policy not in {"error", "unit"}:
            raise ValueError("low_valence_policy must be either 'error' or 'unit'.")
        self.normalization = normalization
        self.low_valence_policy = low_valence_policy

        vertices_raw = data.get("vertices")
        if not isinstance(vertices_raw, list) or not vertices_raw:
            raise ValueError("Input JSON needs a non-empty 'vertices' list.")

        self.n = len(vertices_raw)
        self.data = data
        self.initial_positions = [to_complex(v) for v in vertices_raw]
        for idx, z in enumerate(self.initial_positions):
            ensure_in_disk(z, f"vertex {idx}")

        self.edges = build_edges(
            faces=data.get("faces"),
            edges=data.get("edges"),
            n_vertices=self.n,
        )
        if not self.edges:
            raise ValueError("Mesh has no edges.")

        self.parent: Dict[int, int] = {}
        self.parent_iso: Dict[int, MobiusIsometry] = {}
        self.constraints: List[Constraint] = []
        constraints_raw = data.get("constraints", [])
        if not isinstance(constraints_raw, list):
            raise ValueError("'constraints' must be a list if provided.")

        seen_slaves = set()
        for item in constraints_raw:
            if not isinstance(item, dict):
                raise ValueError(f"Constraint must be an object, got {item!r}")
            slave = int(item["slave"])
            master = int(item["master"])
            if slave == master:
                raise ValueError("A constrained vertex cannot reference itself.")
            if slave in seen_slaves:
                raise ValueError(f"Vertex {slave} has more than one constraint.")
            seen_slaves.add(slave)
            if not (0 <= slave < self.n and 0 <= master < self.n):
                raise ValueError(f"Constraint indices out of range: {item!r}")
            isometry = MobiusIsometry(
                to_complex(item["a"]),
                to_complex(item["b"]),
            )
            self.parent[slave] = master
            self.parent_iso[slave] = isometry
            self.constraints.append(Constraint(slave, master, isometry))

        self._check_constraint_forest()
        self.roots = [i for i in range(self.n) if i not in self.parent]
        self.root_index = {v: k for k, v in enumerate(self.roots)}
        self.root_for_vertex: Dict[int, int] = {}
        self.root_to_vertex: Dict[int, MobiusIsometry] = {}
        for vertex in range(self.n):
            self.root_for_vertex[vertex], self.root_to_vertex[vertex] = (
                self._root_and_transform(vertex)
            )

        self.copies_by_root: Dict[int, List[int]] = {root: [] for root in self.roots}
        for vertex, root in self.root_for_vertex.items():
            self.copies_by_root[root].append(vertex)

        (
            self.attachment_root_for_vertex,
            self.attachment_root_to_vertex,
            self.attachment_copies_by_root,
        ) = self._build_attachment_components()

        self.incident_edges: Dict[int, List[Tuple[int, Edge]]] = {
            i: [] for i in range(self.n)
        }
        for edge in self.edges:
            i, j = edge
            self.incident_edges[i].append((j, edge))
            self.incident_edges[j].append((i, edge))

        self.positions = self.resolve_positions(positions)

    def resolve_positions(
        self, positions: Optional[Sequence[PointInput]] = None
    ) -> List[complex]:
        """Resolve full vertex positions, applying slave/master constraints."""

        if positions is None:
            seed = self.initial_positions
        else:
            seed = [
                z if isinstance(z, complex) else to_complex(z)
                for z in positions
            ]

        if len(seed) == self.n:
            root_positions = [seed[root] for root in self.roots]
        elif len(seed) == len(self.roots):
            root_positions = list(seed)
        else:
            raise ValueError(
                "positions must contain either all vertices or only root vertices."
            )

        resolved: List[complex] = []
        for vertex in range(self.n):
            root = self.root_for_vertex[vertex]
            z = self.root_to_vertex[vertex].apply(root_positions[self.root_index[root]])
            ensure_in_disk(z, f"vertex {vertex}")
            resolved.append(z)
        return resolved

    def compute(
        self,
        positions: Optional[Sequence[PointInput]] = None,
        normalization: Optional[NormalizationMode] = None,
    ) -> List[DirectedWeight]:
        """Return ``[(weight_i_to_j, weight_j_to_i), ...]`` for ``self.edges``."""

        if positions is not None:
            self.positions = self.resolve_positions(positions)
        if normalization is None:
            normalization = self.normalization
        elif normalization not in {"unnormalized", "normalized"}:
            raise ValueError(
                "normalization must be either 'unnormalized' or 'normalized'."
            )

        directed: Dict[Tuple[int, int], float] = {}
        for center in range(self.n):
            directed.update(self._weights_for_center(center, normalization))

        out: List[DirectedWeight] = []
        missing: List[Tuple[int, int]] = []
        for i, j in self.edges:
            wij = directed.get((i, j))
            wji = directed.get((j, i))
            if wij is None:
                missing.append((i, j))
            if wji is None:
                missing.append((j, i))
            if wij is not None and wji is not None:
                out.append((wij, wji))

        if missing:
            sample = ", ".join(f"{i}->{j}" for i, j in missing[:8])
            raise ValueError(f"Could not compute directed weights for: {sample}")
        return out

    def as_json_data(
        self,
        positions: Optional[Sequence[PointInput]] = None,
        normalization: Optional[NormalizationMode] = None,
    ) -> Dict[str, object]:
        """Return a compact JSON-serializable directed-weight payload."""

        if normalization is None:
            normalization = self.normalization
        weights = self.compute(positions, normalization)
        return {
            "edges": [list(edge) for edge in self.edges],
            "directed_edge_weights": [[w0, w1] for w0, w1 in weights],
            "directed_edge_weight_schema": {
                "description": (
                    "For each edge [i, j], weights are "
                    "[weight centered at i, weight centered at j]."
                ),
                "normalization": normalization,
                "formula": "hyperbolic mean value coordinates on attached vertex stars",
            },
            "resolved_vertices": [from_complex(z) for z in self.positions],
        }

    def _weights_for_center(
        self,
        center: int,
        normalization: NormalizationMode,
    ) -> Dict[Tuple[int, int], float]:
        center_position = self.positions[center]
        occurrences = self.attached_star(center)
        if len(occurrences) < 3:
            if self.low_valence_policy == "error":
                raise ValueError(
                    f"Attached star for vertex {center} has valence {len(occurrences)}; "
                    "need at least 3 neighbor occurrences."
                )
            return {
                (center, occurrence.neighbor): 1.0
                for occurrence in occurrences
                if occurrence.source_center == center
            }

        polar_raw: List[Tuple[float, float, StarOccurrence]] = []
        for occurrence in occurrences:
            u = self._move_center_to_origin(
                center_position,
                occurrence.neighbor_position,
            )
            radius = abs(u)
            if radius < EPS:
                raise ValueError(
                    f"Neighbor {occurrence.neighbor} coincides with center {center} "
                    "after attachment."
                )
            angle = math.atan2(u.imag, u.real)
            polar_raw.append((angle, radius, occurrence))

        polar_raw.sort(key=lambda item: item[0])
        angle_groups: List[Tuple[float, List[Tuple[float, StarOccurrence]]]] = []
        for angle, radius, occurrence in polar_raw:
            if angle_groups and abs(angle - angle_groups[-1][0]) <= ANGLE_TOL:
                angle_groups[-1][1].append((radius, occurrence))
            else:
                angle_groups.append((angle, [(radius, occurrence)]))

        if (
            len(angle_groups) > 1
            and star_angle_between(angle_groups[-1][0], angle_groups[0][0]) <= ANGLE_TOL
        ):
            first_angle, first_items = angle_groups[0]
            _, last_items = angle_groups.pop()
            angle_groups[0] = (first_angle, last_items + first_items)

        count = len(angle_groups)
        alphas = []
        for k in range(count):
            next_k = (k + 1) % count
            alphas.append(
                star_angle_between(angle_groups[k][0], angle_groups[next_k][0])
            )

        occurrence_weights: List[Tuple[StarOccurrence, float]] = []
        for k, (_, group_items) in enumerate(angle_groups):
            prev_alpha = alphas[k - 1]
            next_alpha = alphas[k]
            numerator = positive_tan_half_angle(prev_alpha) + positive_tan_half_angle(
                next_alpha
            )
            numerator_share = numerator / len(group_items)
            for radius, occurrence in group_items:
                hyperbolic_radius = 2.0 * math.atanh(min(radius, 1.0 - 1e-15))
                denominator = math.sinh(hyperbolic_radius)
                if abs(denominator) < EPS:
                    raise ValueError(
                        f"Neighbor {occurrence.neighbor} is too close to center {center}."
                    )
                occurrence_weights.append((occurrence, numerator_share / denominator))

        if normalization == "normalized":
            weight_sum = sum(weight for _, weight in occurrence_weights)
            if abs(weight_sum) < EPS:
                raise ValueError(
                    f"Cannot normalize weights for vertex {center}: star sum is zero."
                )
            occurrence_weights = [
                (occurrence, weight / weight_sum)
                for occurrence, weight in occurrence_weights
            ]

        weights: Dict[Tuple[int, int], float] = {}
        for occurrence, weight in occurrence_weights:
            if occurrence.source_center == center:
                weights[(center, occurrence.neighbor)] = weight
        return weights

    def attached_star(self, center: int) -> List[StarOccurrence]:
        """Return all neighbor occurrences attached around ``center``."""

        root = self.attachment_root_for_vertex[center]
        target_from_root = self.attachment_root_to_vertex[center]
        occurrences: List[StarOccurrence] = []

        for source_center in self.attachment_copies_by_root[root]:
            source_from_root = self.attachment_root_to_vertex[source_center]
            target_from_source = compose_mobius(
                target_from_root,
                source_from_root.inverse(),
            )
            for neighbor, edge in self.incident_edges[source_center]:
                occurrences.append(
                    StarOccurrence(
                        source_center=source_center,
                        neighbor=neighbor,
                        neighbor_position=target_from_source.apply(
                            self.positions[neighbor]
                        ),
                        edge=edge,
                    )
                )
        return occurrences

    @staticmethod
    def _move_center_to_origin(center: complex, neighbor: complex) -> complex:
        denominator = 1.0 - center.conjugate() * neighbor
        if abs(denominator) < EPS:
            raise ValueError("Poincare disk isometry denominator is too small.")
        return (neighbor - center) / denominator

    def _check_constraint_forest(self) -> None:
        color = [0] * self.n

        def dfs(vertex: int) -> None:
            color[vertex] = 1
            if vertex in self.parent:
                parent = self.parent[vertex]
                if color[parent] == 1:
                    raise ValueError("Constraint graph contains a cycle.")
                if color[parent] == 0:
                    dfs(parent)
            color[vertex] = 2

        for vertex in range(self.n):
            if color[vertex] == 0:
                dfs(vertex)

    def _root_and_transform(self, vertex: int) -> Tuple[int, MobiusIsometry]:
        if vertex not in self.parent:
            return vertex, IDENTITY_ISOMETRY

        chain = []
        current = vertex
        while current in self.parent:
            chain.append(current)
            current = self.parent[current]
        root = current

        transform = IDENTITY_ISOMETRY
        for child in reversed(chain):
            transform = compose_mobius(self.parent_iso[child], transform)
        return root, transform

    def _build_attachment_components(
        self,
    ) -> Tuple[Dict[int, int], Dict[int, MobiusIsometry], Dict[int, List[int]]]:
        adjacency: Dict[int, List[Tuple[int, MobiusIsometry]]] = {
            vertex: [] for vertex in range(self.n)
        }

        def add_identification(
            source: int,
            target: int,
            source_from_target: MobiusIsometry,
        ) -> None:
            adjacency[target].append((source, source_from_target))
            adjacency[source].append((target, source_from_target.inverse()))

        for constraint in self.constraints:
            add_identification(
                constraint.slave,
                constraint.master,
                constraint.isometry,
            )

        for constraint in self._metadata_corner_constraints():
            add_identification(
                constraint.slave,
                constraint.master,
                constraint.isometry,
            )

        root_for_vertex: Dict[int, int] = {}
        root_to_vertex: Dict[int, MobiusIsometry] = {}
        copies_by_root: Dict[int, List[int]] = {}

        for start in range(self.n):
            if start in root_for_vertex:
                continue
            root = start
            copies_by_root[root] = []
            stack = [(start, IDENTITY_ISOMETRY)]
            while stack:
                vertex, transform = stack.pop()
                if vertex in root_for_vertex:
                    continue
                root_for_vertex[vertex] = root
                root_to_vertex[vertex] = transform
                copies_by_root[root].append(vertex)
                for neighbor, neighbor_from_vertex in adjacency[vertex]:
                    stack.append(
                        (
                            neighbor,
                            compose_mobius(neighbor_from_vertex, transform),
                        )
                    )

        return root_for_vertex, root_to_vertex, copies_by_root

    def _metadata_corner_constraints(self) -> List[Constraint]:
        metadata = self.data.get("metadata")
        if not isinstance(metadata, dict):
            return []

        boundary_meta = metadata.get("boundary_meta")
        pairings = metadata.get("pairings")
        side_subdivisions = metadata.get("side_subdivisions")
        if (
            not isinstance(boundary_meta, list)
            or not isinstance(pairings, list)
            or side_subdivisions is None
        ):
            return []

        boundary_count = len(boundary_meta)
        if boundary_count == 0:
            return []
        boundary_start = self.n - boundary_count
        side_count = sum(
            1 for item in boundary_meta if isinstance(item, list) and item[:1] == ["corner"]
        )
        if side_count <= 0:
            return []
        side_step = boundary_count // side_count
        if side_count * side_step != boundary_count:
            return []

        constraints_by_pair: Dict[Tuple[int, int], MobiusIsometry] = {}
        for constraint in self.constraints:
            slave_meta = boundary_meta[constraint.slave - boundary_start]
            master_meta = boundary_meta[constraint.master - boundary_start]
            if (
                isinstance(slave_meta, list)
                and isinstance(master_meta, list)
                and len(slave_meta) >= 3
                and len(master_meta) >= 3
                and slave_meta[0] == "side"
                and master_meta[0] == "side"
            ):
                constraints_by_pair[(int(slave_meta[1]), int(master_meta[1]))] = (
                    constraint.isometry
                )

        def boundary_vertex(side: int, offset: int) -> int:
            return boundary_start + ((side * side_step + offset) % boundary_count)

        out: List[Constraint] = []
        for item in pairings:
            if not isinstance(item, dict):
                continue
            master_side = int(item["master_side"])
            slave_side = int(item["slave_side"])
            isometry = constraints_by_pair.get((slave_side, master_side))
            if isometry is None:
                continue

            out.append(
                Constraint(
                    boundary_vertex(slave_side, 0),
                    boundary_vertex(master_side, side_step),
                    isometry,
                )
            )
            out.append(
                Constraint(
                    boundary_vertex(slave_side, side_step),
                    boundary_vertex(master_side, 0),
                    isometry,
                )
            )
        return out


def load_json(path: str) -> Dict[str, object]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError("Top-level JSON value must be an object.")
    return data


def write_json(path: str, data: Mapping[str, object]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Path to mesh JSON.")
    parser.add_argument(
        "--output",
        help="Optional JSON output path for edges and directed_edge_weights.",
    )
    parser.add_argument(
        "--normalization",
        choices=("unnormalized", "normalized"),
        default="unnormalized",
        help="Use raw mean-value weights or normalize each attached star to sum to 1.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    calculator = DirectedEdgeWeightCalculator(
        load_json(args.input),
        normalization=args.normalization,
    )
    result = calculator.as_json_data()
    weights = result["directed_edge_weights"]
    flat_weights = [w for pair in weights for w in pair]

    print(f"Normalization: {args.normalization}")
    print(f"Edges: {len(weights)}")
    print(f"Directed weights: {len(flat_weights)}")
    print(f"Min directed weight: {min(flat_weights):.12g}")
    print(f"Max directed weight: {max(flat_weights):.12g}")

    if args.output:
        write_json(args.output, result)
        print(f"Wrote directed weights to {args.output}")


if __name__ == "__main__":
    main()
