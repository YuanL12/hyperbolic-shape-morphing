#!/usr/bin/env python3
"""Discrete harmonic map descent in the Poincare disk.

This is a small, dependency-free implementation of the algorithmic core from
the paper section you pointed to.

What this script does:
1. Reads a mesh embedded in the Poincare disk from JSON.
2. Optionally enforces cut-mesh identifications using disk isometries.
3. Minimizes the discrete Dirichlet energy

       E(f) = 1/2 * sum_{ij in E} c_ij * d_H(f_i, f_j)^2

   by Riemannian gradient descent using the exponential map of the Poincare
   disk.

Input JSON schema
-----------------
{
  "vertices": [[x0, y0], [x1, y1], ...],
  "faces": [[i, j, k], ...],               // optional if "edges" is given
  "edges": [[i, j], ...],                  // optional if "faces" is given
  "edge_weights": [1.0, 2.0, ...],         // optional, defaults to all 1
  "directed_edge_weights": [[1.0, 1.5], ...],
                                             // optional gradient weights
                                             // for edge [i, j]: [at i, at j]
  "fixed": [0, 3, 4],                      // optional, only for root vertices
  "constraints": [                         // optional
    {
      "slave": 7,
      "master": 2,
      "a": [a_re, a_im],
      "b": [b_re, b_im]
    }
  ],
  "iterations": 200,                       // optional
  "step_size": 0.2,                        // optional
  "tolerance": 1e-8                        // optional
}

Each constraint represents the orientation-preserving disk isometry

    gamma(z) = (a z + b) / (conj(b) z + conj(a)),

with |a|^2 - |b|^2 = 1, and enforces

    f_slave = gamma(f_master).

Constraints may form a forest of parent-child relations. Cycles are rejected.

Example
-------
python3 poincare_harmonic_map.py example_disk_mesh.json --output result.json


Explanation of the code
-------------------------
Vertices Partitioned by:
interior/boundary = geometric location
root/slave        = optimization dependency
fixed/free        = whether a root is allowed to move


"""

from __future__ import annotations

import argparse
import cmath
import json
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


EPS = 1e-12


@dataclass(frozen=True)
class MobiusIsometry:
    """Orientation-preserving Poincare disk isometry."""

    a: complex
    b: complex

    def __post_init__(self) -> None:
        det = abs(self.a) ** 2 - abs(self.b) ** 2
        if abs(det - 1.0) > 1e-7:
            raise ValueError(
                "Invalid SU(1,1) parameters: expected |a|^2 - |b|^2 = 1, "
                f"got {det:.12f}"
            )

    def apply(self, z: complex) -> complex:
        # apply the mobius isometry to the point z
        den = self.b.conjugate() * z + self.a.conjugate()
        if abs(den) < EPS:
            raise ValueError("Mobius transform denominator is too small.")
        return (self.a * z + self.b) / den

    def derivative(self, z: complex) -> complex:
        # derivative of the mobius isometry at the point z
        den = self.b.conjugate() * z + self.a.conjugate()
        if abs(den) < EPS:
            raise ValueError("Mobius derivative denominator is too small.")
        return 1.0 / (den * den)

    def push_forward(self, z: complex, v: complex) -> complex:
        # push the vector v forward by the mobius isometry at the point z
        return self.derivative(z) * v

    def inverse(self) -> "MobiusIsometry":
        return MobiusIsometry(self.a.conjugate(), -self.b)

    def pull_back_from_image(self, image_z: complex, image_v: complex) -> complex:
        inv = self.inverse()
        return inv.push_forward(image_z, image_v)


def to_complex(xy: Sequence[float]) -> complex:
    if len(xy) != 2:
        raise ValueError(f"Expected a 2-vector, got {xy!r}")
    return complex(float(xy[0]), float(xy[1]))


def from_complex(z: complex) -> List[float]:
    return [float(z.real), float(z.imag)]


def ensure_in_disk(z: complex, name: str = "point") -> None:
    if abs(z) >= 1.0:
        raise ValueError(f"{name} must lie strictly inside the unit disk: {z!r}")


def mobius_add(x: complex, y: complex) -> complex:
    den = 1.0 + x.conjugate() * y
    if abs(den) < EPS:
        raise ValueError("Mobius addition denominator is too small.")
    return (x + y) / den


def lambda_x(x: complex) -> float:
    denom = 1.0 - abs(x) ** 2
    if denom <= 0.0:
        raise ValueError("Point left the Poincare disk.")
    return 2.0 / denom


def hyperbolic_distance(x: complex, y: complex) -> float:
    delta = abs(mobius_add(-x, y))
    delta = min(max(delta, 0.0), 1.0 - 1e-15)
    return 2.0 * math.atanh(delta)


def log_map(x: complex, y: complex) -> complex:
    delta = mobius_add(-x, y)
    norm = abs(delta)
    if norm < EPS:
        return 0.0j
    scale = (2.0 / lambda_x(x)) * math.atanh(min(norm, 1.0 - 1e-15)) / norm
    return scale * delta


def exp_map(x: complex, v: complex) -> complex:
    norm_v = abs(v)
    if norm_v < EPS:
        return x
    factor = math.tanh(lambda_x(x) * norm_v / 2.0) / norm_v
    y = mobius_add(x, factor * v)
    if abs(y) >= 1.0:
        # Numerical guard near the boundary.
        y *= (1.0 - 1e-12) / abs(y)
    return y


def tangent_norm(x: complex, v: complex) -> float:
    return lambda_x(x) * abs(v)


def build_edges(
    faces: Optional[Sequence[Sequence[int]]],
    edges: Optional[Sequence[Sequence[int]]],
    n_vertices: int,
) -> List[Tuple[int, int]]:
    if edges is not None:
        out: Set[Tuple[int, int]] = set()
        for edge in edges:
            if len(edge) != 2:
                raise ValueError(f"Edge must have two endpoints, got {edge!r}")
            i, j = int(edge[0]), int(edge[1])
            if i == j:
                continue
            if not (0 <= i < n_vertices and 0 <= j < n_vertices):
                raise ValueError(f"Edge index out of range: {edge!r}")
            out.add((min(i, j), max(i, j)))
        return sorted(out)

    if faces is None:
        raise ValueError("Need either 'faces' or 'edges' in the input JSON.")

    out: Set[Tuple[int, int]] = set()
    for face in faces:
        if len(face) != 3:
            raise ValueError(f"Face must have three vertices, got {face!r}")
        i, j, k = (int(face[0]), int(face[1]), int(face[2]))
        for u, v in ((i, j), (j, k), (k, i)):
            if u == v:
                continue
            if not (0 <= u < n_vertices and 0 <= v < n_vertices):
                raise ValueError(f"Face index out of range: {face!r}")
            out.add((min(u, v), max(u, v)))
    return sorted(out)


class HarmonicMapSolver:
    def __init__(self, data: Dict[str, object]) -> None:
        vertices_raw = data.get("vertices")
        if not isinstance(vertices_raw, list) or not vertices_raw:
            raise ValueError("Input JSON needs a non-empty 'vertices' list.")

        self.n = len(vertices_raw)
        self.root_positions: Dict[int, complex] = {}
        self.parent: Dict[int, int] = {}
        self.parent_iso: Dict[int, MobiusIsometry] = {}

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

        edge_weights_raw = data.get("edge_weights")
        if edge_weights_raw is None:
            self.edge_weights = [1.0] * len(self.edges)
        else:
            if not isinstance(edge_weights_raw, list):
                raise ValueError("'edge_weights' must be a list if provided.")
            if len(edge_weights_raw) != len(self.edges):
                raise ValueError(
                    f"'edge_weights' has length {len(edge_weights_raw)} but there are "
                    f"{len(self.edges)} edges."
                )
            self.edge_weights = [float(w) for w in edge_weights_raw]

        directed_edge_weights_raw = data.get("directed_edge_weights")
        if directed_edge_weights_raw is None:
            self.directed_edge_weights = [
                (weight, weight) for weight in self.edge_weights
            ]
            self.uses_directed_edge_weights = False
        else:
            if not isinstance(directed_edge_weights_raw, list):
                raise ValueError(
                    "'directed_edge_weights' must be a list if provided."
                )
            if len(directed_edge_weights_raw) != len(self.edges):
                raise ValueError(
                    "'directed_edge_weights' has length "
                    f"{len(directed_edge_weights_raw)} but there are "
                    f"{len(self.edges)} edges."
                )
            self.directed_edge_weights = []
            for edge, weights in zip(self.edges, directed_edge_weights_raw):
                if not isinstance(weights, list) or len(weights) != 2:
                    raise ValueError(
                        "Each directed edge weight must be a two-item list "
                        f"[weight_at_i, weight_at_j] for edge {edge}, got "
                        f"{weights!r}."
                    )
                self.directed_edge_weights.append(
                    (float(weights[0]), float(weights[1]))
                )
            self.uses_directed_edge_weights = True

        self.edge_energy_weights = [
            0.5 * (weight_i + weight_j)
            for weight_i, weight_j in self.directed_edge_weights
        ]

        self.fixed: Set[int] = set(int(v) for v in data.get("fixed", []))
        self.iterations = int(data.get("iterations", 200))
        self.step_size = float(data.get("step_size", 0.2))
        self.tolerance = float(data.get("tolerance", 1e-8))

        constraints_raw = data.get("constraints", [])
        if not isinstance(constraints_raw, list):
            raise ValueError("'constraints' must be a list if provided.")

        seen_slaves: Set[int] = set()
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
            a = to_complex(item["a"])
            b = to_complex(item["b"])
            self.parent[slave] = master
            self.parent_iso[slave] = MobiusIsometry(a, b)

        self.roots = [i for i in range(self.n) if i not in self.parent]
        if not self.roots:
            raise ValueError("At least one unconstrained root vertex is required.")

        for idx in self.fixed:
            if idx in self.parent:
                raise ValueError(
                    f"Fixed vertex {idx} is constrained. Fix only root vertices."
                )

        self._check_constraint_forest()
        self.positions = self._resolve_all_positions(self.initial_positions)
        self.energy_history: List[float] = []

    def _check_constraint_forest(self) -> None:
        color = [0] * self.n

        def dfs(v: int) -> None:
            color[v] = 1
            if v in self.parent:
                u = self.parent[v]
                if color[u] == 1:
                    raise ValueError("Constraint graph contains a cycle.")
                if color[u] == 0:
                    dfs(u)
            color[v] = 2

        for v in range(self.n):
            if color[v] == 0:
                dfs(v)

    def _resolve_vertex(
        self,
        idx: int,
        roots: Sequence[complex],
        cache: Dict[int, complex],
        root_index: Dict[int, int],
    ) -> complex:
        if idx in cache:
            return cache[idx]
        # recursively resolve the parent vertex
        if idx in self.parent:
            parent = self.parent[idx]
            parent_z = self._resolve_vertex(parent, roots, cache, root_index)
            z = self.parent_iso[idx].apply(parent_z)
        else:
            # if the vertex is not constrained, use the root position
            z = roots[root_index[idx]]
        ensure_in_disk(z, f"vertex {idx}")
        # cache the resolved position
        cache[idx] = z
        return z

    def _resolve_all_positions(self, seed_positions: Sequence[complex]) -> List[complex]:
        """
        Resolve all vertex positions from the seed/root positions.
        If any vertex is a slave, its position is resolved from the master's position by the mobius isometry.
        """
        if len(seed_positions) == self.n:
            root_values = [seed_positions[v] for v in self.roots]
        elif len(seed_positions) == len(self.roots):
            root_values = list(seed_positions)
        else:
            raise ValueError(
                "seed_positions must contain either all vertex positions or only "
                "the root-vertex positions."
            )
        root_index = {v: k for k, v in enumerate(self.roots)}
        cache: Dict[int, complex] = {}
        return [self._resolve_vertex(i, root_values, cache, root_index) for i in range(self.n)]

    def current_root_positions(self) -> List[complex]:
        return [self.positions[v] for v in self.roots]

    def energy_and_raw_gradient(
        self, positions: Sequence[complex]
    ) -> Tuple[float, List[complex]]:
        energy = 0.0
        grad = [0.0j for _ in range(self.n)]
        for (i, j), energy_weight, (weight_i, weight_j) in zip(
            self.edges,
            self.edge_energy_weights,
            self.directed_edge_weights,
        ):
            zi = positions[i]
            zj = positions[j]
            d = hyperbolic_distance(zi, zj)
            energy += 0.5 * energy_weight * d * d
            grad[i] -= weight_i * log_map(zi, zj)
            grad[j] -= weight_j * log_map(zj, zi)
        return energy, grad

    def _pull_gradient_to_roots(
        self, positions: Sequence[complex], raw_grad: Sequence[complex]
    ) -> Dict[int, complex]:
        acc = {v: 0.0j for v in self.roots}

        def push_to_root(idx: int, grad_v: complex) -> None:
            if idx in self.parent: # not a root
                pulled = self.parent_iso[idx].pull_back_from_image(
                    positions[idx], grad_v
                )
                push_to_root(self.parent[idx], pulled)
            else: # root
                acc[idx] += grad_v

        for idx, g in enumerate(raw_grad):
            push_to_root(idx, g)
        return acc

    def _max_hyperbolic_step(
        self, old_positions: Sequence[complex], new_positions: Sequence[complex]
    ) -> float:
        return max(
            hyperbolic_distance(z_old, z_new)
            for z_old, z_new in zip(old_positions, new_positions)
        )

    def solve(self) -> Dict[str, object]:
        root_positions = self.current_root_positions()
        last_grad_norm = None
        # free roots are the roots that are not fixed
        free_roots = [root for root in self.roots if root not in self.fixed]

        for _ in range(self.iterations):
            positions = self._resolve_all_positions(root_positions)
            energy, raw_grad = self.energy_and_raw_gradient(positions)
            root_grad = self._pull_gradient_to_roots(positions, raw_grad)
            self.energy_history.append(energy)

            if free_roots:
                grad_norm_sq = 0.0
                for root in free_roots:
                    grad_norm_sq += tangent_norm(positions[root], root_grad[root]) ** 2
                last_grad_norm = math.sqrt(grad_norm_sq / len(free_roots))
            else:
                last_grad_norm = 0.0

            if last_grad_norm < self.tolerance:
                root_positions = [positions[v] for v in self.roots]
                break

            step = self.step_size
            accepted = False
            for _ in range(25):
                candidate_roots: List[complex] = []
                for root in self.roots:
                    z = positions[root]
                    if root in self.fixed:
                        candidate_roots.append(z)
                    else:
                        candidate_roots.append(exp_map(z, -step * root_grad[root]))

                candidate_positions = self._resolve_all_positions(candidate_roots)
                candidate_energy, _ = self.energy_and_raw_gradient(candidate_positions)
                if candidate_energy <= energy:
                    root_positions = candidate_roots
                    max_step = self._max_hyperbolic_step(positions, candidate_positions)
                    accepted = True
                    if max_step < self.tolerance:
                        self.positions = candidate_positions
                        self.energy_history.append(candidate_energy)
                        return self._result(candidate_energy, last_grad_norm)
                    break
                step *= 0.5

            if not accepted:
                root_positions = [positions[v] for v in self.roots]
                break

        self.positions = self._resolve_all_positions(root_positions)
        final_energy, _ = self.energy_and_raw_gradient(self.positions)
        return self._result(final_energy, last_grad_norm)

    def _result(self, final_energy: float, final_grad_norm: Optional[float]) -> Dict[str, object]:
        return {
            "vertices": [from_complex(z) for z in self.positions],
            "energy": final_energy,
            "mean_free_gradient_norm": final_grad_norm,
            "iterations_recorded": len(self.energy_history),
            "energy_history": self.energy_history,
            "edges": [list(edge) for edge in self.edges],
            "edge_weights": self.edge_weights,
            "directed_edge_weights": [
                [weight_i, weight_j]
                for weight_i, weight_j in self.directed_edge_weights
            ],
            "uses_directed_edge_weights": self.uses_directed_edge_weights,
            "roots": self.roots,
            "fixed": sorted(self.fixed),
        }


def load_json(path: str) -> Dict[str, object]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError("Top-level JSON value must be an object.")
    return data


def write_json(path: str, data: Dict[str, object]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Path to the input JSON file.")
    parser.add_argument(
        "--output",
        help="Optional output JSON path. Defaults to printing a short summary only.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = load_json(args.input)
    solver = HarmonicMapSolver(data)
    result = solver.solve()

    if args.output:
        write_json(args.output, result)

    print(f"Final energy: {result['energy']:.12f}")
    print(
        "Mean free gradient norm: "
        f"{float(result['mean_free_gradient_norm'] or 0.0):.12e}"
    )
    print(f"Recorded iterations: {result['iterations_recorded']}")
    if args.output:
        print(f"Wrote result to {args.output}")


if __name__ == "__main__":
    main()
