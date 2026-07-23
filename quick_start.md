## Quick start: genus-2 morph

The project uses [uv](https://docs.astral.sh/uv/) to provide a reproducible
environment:

```bash
uv sync
```

Generate the source embedding $Z^0$, a symmetric genus-2 triangulation in a
regular octagonal fundamental domain:

```bash
uv run python examples/generate_genus2_bolza_symmetric_embedding.py \
  --output examples/input/Z0_genus2_bolza_symmetric_center_fan.json
```

Create a weighted endpoint input and solve for its harmonic embedding $Z^1$:

```bash
uv run python examples/make_sin8_endpoint_input.py \
  examples/input/Z0_genus2_bolza_symmetric_center_fan.json \
  --output examples/input/Z1_genus2_bolza_symmetric_center_fan_sin8_input.json \
  --min-weight 1 \
  --max-weight 50 \
  --frequency 8

uv run hyper-morph-solve \
  examples/input/Z1_genus2_bolza_symmetric_center_fan_sin8_input.json \
  --output examples/input/Z1_genus2_bolza_symmetric_center_fan_sin8_solution.json \
  --iterations 20000 \
  --step-size 0.004 \
  --tolerance 1e-10
```

Interpolate the directed mean-value weights of the two embeddings and solve 50
morph frames:

```bash
uv run hyper-morph \
  --source-embedding examples/input/Z0_genus2_bolza_symmetric_center_fan.json \
  --target-embedding examples/input/Z1_genus2_bolza_symmetric_center_fan_sin8_solution.json \
  --target-directed-weights mean_value \
  --output-dir examples/output/frames/genus2_sin8 \
  --frames 50 \
  --iterations 5000 \
  --step-size 0.01 \
  --tolerance 1e-9 \
  --reference-fundamental-domain

uv run python utils/visualize_morph_frames.py \
  examples/output/frames/genus2_sin8 \
  --output examples/output_html/genus2_sin8.html
```

## Input JSON

The generated genus-2 input has the following structure. This is abridged:
`...` marks omitted values and is not literal JSON.

```json
{
  "name": "genus2_bolza_symmetric_center_fan",
  "vertices": [
    [0.0, 0.0],
    [0.2802988050845715, 0.0],
    ...
  ],
  "faces": [
    [0, 1, 2],
    ...
  ],
  "edges": [
    [0, 1],
    [0, 2],
    ...
  ],
  "fixed": [81, 86, 91, 96, 101, 106, 111, 116],
  "constraints": [
    {
      "slave": 92,
      "master": 85,
      "a": [-1.7071067811865488, 1.7071067811865488],
      "b": [2.0301035302564374, -0.8408964152537153]
    },
    ...
  ]
}
```

### Fields

| Field | Meaning |
| --- | --- |
| `name`, `description` | Optional human-readable labels. The solver ignores them; put data that should be returned under `metadata`. |
| `vertices` | Required list of 2D `[x, y]` coordinates. Every point must lie strictly inside the unit disk. Vertex indices are positions in this list. |
| `faces` | Face-boundary cycles. If `edges` is absent, every face must be a triangle `[i, j, k]` so the solver can derive the edges. Quadrilateral and pentagonal faces are only supported when an explicit `edges` list is also present. |
| `edges` | Undirected pairs `[i, j]`. Required when `faces` is absent or contains non-triangular faces. The solver removes duplicates, writes each pair as `(min(i, j), max(i, j))`, and sorts the result lexicographically. |
| `edge_weights` | One scalar per canonical edge. Defaults to `1.0`. |
| `directed_edge_weights` | One pair `[w_ij, w_ji]` for every canonical edge `[i, j]`: the first weight acts at `i` toward `j`, and the second acts at `j` toward `i`. When absent, the scalar `edge_weights` value is used in both directions. |
| `fixed` | Indices of unconstrained root vertices that must not move. Defaults to `[]`. Fixed vertices do not have to be boundary vertices. A constraint slave cannot also be fixed. |
| `constraints` | Optional slave/master identifications. Each item enforces `z_slave = (a z_master + b) / (conj(b) z_master + conj(a))`, where `a` and `b` are complex numbers stored as `[real, imaginary]`. Constraints must form a forest. |
| `metadata` | Optional JSON object copied into the result without affecting the solve. |

If both `faces` and `edges` are present, `edges` defines the solver topology.
Weight arrays must follow the solver's canonical sorted edge order.

### Solver and morph arguments

There are also arguments for optimization and rendering. The quick-start commands above show where to specify them.

| Argument | Meaning |
| --- | --- |
| `--iterations` | Maximum descent iterations. |
| `--step-size` | Initial Riemannian descent step. |
| `--tolerance` | Stopping tolerance for the selected convergence criterion. |
| `--line-search-objective` | `"energy"`, `"gradient_norm"`, or `"none"`. |
| `--convergence-criterion` | `"gradient_norm"` or `"relative_step"`. |
| `--edge-force-model` | Edge potential and force used by `hyper-morph-solve`: `"squared_distance"` $\sum \frac{1}{2} w d^2$ or `"hyperbolic_mean_value"` $\sum w (\cosh(d) - 1)$. |
| `--reference-fundamental-domain` | Tells `hyper-morph` to include a fixed drawing reference in its output frames. The reference is derived from the source embedding's fixed corner positions; it is not an input JSON field. |

Run `uv run hyper-morph-solve --help` or `uv run hyper-morph --help` for
command-specific defaults and the complete argument list.

### Quadrilateral and pentagonal faces

The repository includes polygonal examples for both surfaces:

- [`examples/generate_genus2_quad_pent_graph.py`](examples/generate_genus2_quad_pent_graph.py)
- [`examples/generate_genus3_quad_pent_graph.py`](examples/generate_genus3_quad_pent_graph.py)

These generators coarsen a triangulation by removing construction edges. Their
output stores quadrilateral and pentagonal boundary cycles under `faces` and
also stores the complete canonical graph under `edges`. The solver uses the
explicit edges; polygonal faces remain available for validation, rendering, and
the returned result. Run either script with `--help` to see its base-embedding
and output arguments.

### Directed weights on undirected edges

The topology does not need two stored arrows for each edge. Direction lives in
the two weights associated with one undirected edge. For example:

```json
{
  "edges": [[2, 5], [3, 5]],
  "directed_edge_weights": [
    [0.3, 0.7],
    [0.4, 0.6]
  ]
}
```

For edge `[2, 5]`, `0.3` is \(w_{2,5}\), used in the equation at vertex `2`,
and `0.7` is \(w_{5,2}\), used in the equation at vertex `5`. In code, their
contributions are:

```text
gradient[2] -= 0.3 * edge_force(2, 5)
gradient[5] -= 0.7 * edge_force(5, 2)
```

The solver first canonicalizes all input edges and then pairs the weight rows
with that canonical list by index. It does not reorder the weight rows or swap
the two values inside a row. Therefore, when weights are supplied manually:

1. store every edge as `[min(i, j), max(i, j)]`;
2. sort the edge list lexicographically; and
3. construct `edge_weights` and `directed_edge_weights` in exactly that order.

`hyper-morph-weights` avoids this bookkeeping: its output contains an `edges`
list and a matching `directed_edge_weights` list that should be kept together.

### What `--edge-force-model` changes

The selected edge-force model is used during every energy and raw-gradient
evaluation. For a canonical edge `[i, j]`, let

```text
d = hyperbolic_distance(z_i, z_j)
w_bar = (w_ij + w_ji) / 2
```

With `"squared_distance"`, the edge contribution is

```text
energy       += 1/2 * w_bar * d^2
gradient[i]  -= w_ij * log_z_i(z_j)
gradient[j]  -= w_ji * log_z_j(z_i)
```

With `"hyperbolic_mean_value"`, it becomes

```text
energy       += w_bar * (cosh(d) - 1)
scale         = sinh(d) / d
gradient[i]  -= w_ij * scale * log_z_i(z_j)
gradient[j]  -= w_ji * scale * log_z_j(z_i)
```

The limiting scale is `1` when `d` is zero. Because the hyperbolic norm of the
log-map vector is `d`, the first model has force magnitude proportional to
`d`, while the second has magnitude proportional to `sinh(d)`. The models are
similar for short edges but the hyperbolic mean-value force grows much faster
for long edges.

The scalar energy uses the average `w_bar`, while the two endpoint equations
use `w_ij` and `w_ji` separately. If those two directed weights differ, the
reported symmetric energy is a diagnostic and the directed residual is not, in
general, its exact gradient. This is why the morph command sets
`edge_force_model` to `"hyperbolic_mean_value"` and defaults its line-search
objective to `"none"`. The standalone solver defaults to
`"squared_distance"` and an energy-based line search unless the CLI overrides
them.

### Boundary vertices may move

“Boundary” describes a vertex's position in the cut-open fundamental domain;
it does not imply that the solver fixes it. For the same genus-2 mesh, this is
a valid alternative input:

```json
{
  "vertices": [
    [0.0, 0.0],
    [0.2802988050845715, 0.0],
    ...
  ],
  "edges": [
    [0, 1],
    ...
  ],
  "fixed": [0, 1],
  "constraints": [
    {
      "slave": 92,
      "master": 85,
      "a": [-1.7071067811865488, 1.7071067811865488],
      "b": [2.0301035302564374, -0.8408964152537153]
    },
    ...
  ]
}
```

Here vertices `0` and `1` are interior roots used to anchor the embedding.
None of the boundary vertices `81` through `120` are fixed: boundary roots
are optimized, and boundary slaves follow their masters through the Möbius
constraints. `fixed` may also be empty, although an unanchored problem can
retain a global disk-isometry ambiguity.

## Python API and solver result

```python
import json

from hyper_morph import HarmonicMapSolver

with open(
    "examples/input/Z1_genus2_bolza_symmetric_center_fan_sin8_input.json",
    encoding="utf-8",
) as file:
    mesh = json.load(file)

result = HarmonicMapSolver(mesh).solve()

print(result["stop_reason"])
print(result["energy"])
print(result["mean_free_gradient_norm"])
print(result["vertices"])
```

`solve()` does not mutate `mesh`. It returns a JSON-serializable dictionary:

| Result field | Meaning |
| --- | --- |
| `vertices` | Final positions for every vertex, including slave positions resolved from their masters. |
| `energy` | Final scalar energy evaluated at `vertices`. |
| `mean_free_gradient_norm` | Mean hyperbolic gradient norm over free root vertices at the last recorded iteration; it is `null` when no iteration ran. |
| `iterations_recorded` | Number of entries stored in `energy_history`. |
| `energy_history` | Energy measured during descent, useful for plotting or diagnosing convergence. |
| `stop_reason` | Why the solve ended: `gradient_tolerance`, `step_tolerance`, `max_iterations`, `line_search_failed`, or `zero_iterations`. |
| `roots` | Vertices optimized independently; constraint slaves are omitted. |
| `fixed`, `constraints` | The fixed-root set and side-pairing constraints used by the solve. |
| `edges`, `faces` | The resolved edge list and, when supplied, the original faces. |
| `edge_weights`, `directed_edge_weights` | The weights actually used by the solver. |
| `edge_force_model`, `line_search_objective`, `convergence_criterion` | The effective solver options, including defaults. |
| `metadata` | The input metadata copied unchanged. |

In particular, check `stop_reason` before treating the output as converged:
`max_iterations` and `line_search_failed` return the best final state reached
but do not indicate that the requested tolerance was met.

To compute directed mean-value weights without solving:

```bash
uv run hyper-morph-weights \
  examples/input/Z0_genus2_bolza_symmetric_center_fan.json \
  --normalization normalized \
  --output examples/output/Z0_genus2_directed_weights.json
```
