# Graph Morphing Problem on Hyperbolic Surface

## Theorem
Tutte Embedding works for nonsymmetric weights on non-triangulation(just planar graph). 


## Current Experiment

The current experiment produces four mean-value morph viewers:

1. genus-2 Bolza triangulation,
2. genus-2 Bolza quad/pent cell graph,
3. genus-3 Klein quartic triangulation,
4. genus-3 Klein quartic quad/pent cell graph.

Each morph starts from two convex endpoint embeddings with the same topology.
The endpoints are produced by minimizing with scalar edge weights:

- all edges have weight `1`,
- half the edges have weight `10` and half have weight `1`.

The morph script then computes normalized directed mean-value weights from each
endpoint embedding, interpolates those directed weights, and solves one
hyperbolic harmonic map problem per frame. The local update uses the
squared-distance/log-map gradient for

```text
F_i(y) = 1/2 sum_j omega_ij d_D(y, z_j)^2.
```

## Reproduce the Four Morph Viewers

Run the commands below from the repository root. They assume the endpoint
solution JSON files already exist in `input/`.

### Morph Script Arguments

The four commands below use the same main arguments:

- `--source-embedding`: the source endpoint JSON. The script takes the graph
  topology, constraints, fixed vertices, and frame-0 vertex positions from this
  file, so this must be an embedding.
- `--target-embedding`: the target endpoint JSON. Its vertices are used to
  compute the target mean-value directed weights. It must have the same vertex
  order and topology as `--source-embedding`.
- `--output-dir`: directory where `frame_*.json` morph frames are written.
- `--frames`: number of morph frames.
- `--start-directed-weights mean_value`: compute directed mean-value weights
  from `--source-embedding` for frame 0.
- `--target-directed-weights mean_value`: compute directed mean-value weights
  from `--target-embedding` for the final frame.
- `--normalization normalized`: normalize the outgoing mean-value weights at
  each vertex.
- `--iterations`, `--step-size`, `--tolerance`: solver iteration budget,
  gradient step size, and stopping tolerance for each frame.
- `--line-search-objective none`: disable line-search energy minimization
  during morphing.
- `--convergence-criterion relative_step`: stop each frame solve using the
  relative update size rather than gradient norm.
- `--skip-face-orientation-check`: validate edge crossings without requiring
  every face orientation check, which is useful for the quad/pent cell graphs.

The `visualize_morph_frames.py` command takes the frame directory as its first
argument and writes a self-contained HTML viewer to `--output`.

### Genus-2 Triangulation

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --source-embedding input/example_genus2_bolza_delaunay_all_1_solution.json \
  --target-embedding input/example_genus2_bolza_delaunay_half_10_solution.json \
  --output-dir boundary_morph_genus2_bolza_triangulation_normalized_F_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --normalization normalized \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check

.venv/bin/python visualize_morph_frames.py \
  boundary_morph_genus2_bolza_triangulation_normalized_F_frames \
  --output morph_viewer_genus2_bolza_triangulation_normalized_F.html
```

### Genus-2 Cells

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --source-embedding input/example_genus2_bolza_convex_quad_pent_graph_all_1_solution.json \
  --target-embedding input/example_genus2_bolza_convex_quad_pent_graph_half_10_solution.json \
  --output-dir boundary_morph_genus2_bolza_cells_normalized_F_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --normalization normalized \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check

.venv/bin/python visualize_morph_frames.py \
  boundary_morph_genus2_bolza_cells_normalized_F_frames \
  --output morph_viewer_genus2_bolza_cells_normalized_F.html
```

### Genus-3 Triangulation

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --source-embedding input/example_genus3_klein_quartic_triangulation_all_1_solution.json \
  --target-embedding input/example_genus3_klein_quartic_triangulation_half_10_solution.json \
  --output-dir boundary_morph_genus3_triangulation_normalized_F_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --normalization normalized \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check

.venv/bin/python visualize_morph_frames.py \
  boundary_morph_genus3_triangulation_normalized_F_frames \
  --output morph_viewer_genus3_triangulation_normalized_F.html
```

### Genus-3 Cells

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --source-embedding input/example_genus3_klein_quartic_convex_quad_pent_graph_all_1_solution.json \
  --target-embedding input/example_genus3_klein_quartic_convex_quad_pent_graph_half_10_solution.json \
  --output-dir boundary_morph_genus3_cells_normalized_F_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --normalization normalized \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check

.venv/bin/python visualize_morph_frames.py \
  boundary_morph_genus3_cells_normalized_F_frames \
  --output morph_viewer_genus3_cells_normalized_F.html
```

## Generate Inputs From Code

The four morph commands above use endpoint solution files in `input/`. A fresh
clone can regenerate those files from code. More detailed notes are recorded in:

- `notes/genus2_bolza_mean_value_morph_generation.md`
- `notes/genus3_mean_value_morph_generation.md`

The four graph embeddings used by the morphs are generated first:

```bash
.venv/bin/python generate_genus2_bolza_embedding.py

.venv/bin/python poincare_harmonic_map.py \
  input/example_genus2_bolza_delaunay_irregular.json \
  --output input/example_genus2_bolza_delaunay_all_1_solution.json \
  --iterations 8000 \
  --step-size 0.006 \
  --tolerance 1e-9

.venv/bin/python generate_genus2_quad_pent_graph.py \
  --quads 25 \
  --pentagons 6 \
  --seed-limit 700

.venv/bin/python generate_genus3_klein_embedding.py

.venv/bin/python generate_genus3_quad_pent_graph.py
```

These commands write:

```text
input/example_genus2_bolza_delaunay_irregular.json
input/example_genus2_bolza_convex_quad_pent_graph.json
input/example_genus3_klein_quartic_embedded_irregular.json
input/example_genus3_klein_quartic_convex_quad_pent_graph.json
```

For each graph, solve the all-one endpoint:

```bash
.venv/bin/python poincare_harmonic_map.py \
  input/example_genus2_bolza_convex_quad_pent_graph.json \
  --output input/example_genus2_bolza_convex_quad_pent_graph_all_1_solution.json \
  --iterations 8000 \
  --step-size 0.006 \
  --tolerance 1e-9

.venv/bin/python poincare_harmonic_map.py \
  input/example_genus3_klein_quartic_embedded_irregular.json \
  --output input/example_genus3_klein_quartic_triangulation_all_1_solution.json \
  --iterations 8000 \
  --step-size 0.006 \
  --tolerance 1e-9

.venv/bin/python poincare_harmonic_map.py \
  input/example_genus3_klein_quartic_convex_quad_pent_graph.json \
  --output input/example_genus3_klein_quartic_convex_quad_pent_graph_all_1_solution.json \
  --iterations 8000 \
  --step-size 0.006 \
  --tolerance 1e-9
```

Then create the half-10 endpoint inputs and solve them:

```bash
.venv/bin/python make_half10_endpoint_input.py \
  input/example_genus2_bolza_delaunay_irregular.json \
  --output input/example_genus2_bolza_delaunay_half_10.json

.venv/bin/python poincare_harmonic_map.py \
  input/example_genus2_bolza_delaunay_half_10.json \
  --output input/example_genus2_bolza_delaunay_half_10_solution.json \
  --iterations 8000 \
  --step-size 0.006 \
  --tolerance 1e-9

.venv/bin/python make_half10_endpoint_input.py \
  input/example_genus2_bolza_convex_quad_pent_graph.json \
  --output input/example_genus2_bolza_convex_quad_pent_graph_half_10.json

.venv/bin/python poincare_harmonic_map.py \
  input/example_genus2_bolza_convex_quad_pent_graph_half_10.json \
  --output input/example_genus2_bolza_convex_quad_pent_graph_half_10_solution.json \
  --iterations 8000 \
  --step-size 0.006 \
  --tolerance 1e-9

.venv/bin/python make_half10_endpoint_input.py \
  input/example_genus3_klein_quartic_embedded_irregular.json \
  --output input/example_genus3_klein_quartic_triangulation_half_10.json

.venv/bin/python poincare_harmonic_map.py \
  input/example_genus3_klein_quartic_triangulation_half_10.json \
  --output input/example_genus3_klein_quartic_triangulation_half_10_solution.json \
  --iterations 8000 \
  --step-size 0.006 \
  --tolerance 1e-9

.venv/bin/python make_half10_endpoint_input.py \
  input/example_genus3_klein_quartic_convex_quad_pent_graph.json \
  --output input/example_genus3_klein_quartic_convex_quad_pent_graph_half_10.json

.venv/bin/python poincare_harmonic_map.py \
  input/example_genus3_klein_quartic_convex_quad_pent_graph_half_10.json \
  --output input/example_genus3_klein_quartic_convex_quad_pent_graph_half_10_solution.json \
  --iterations 8000 \
  --step-size 0.006 \
  --tolerance 1e-9
```

After these endpoint files exist, run the four morph-viewer commands above.

## Paper Figures

After the four frame directories exist, regenerate the SVG grids used by the
TeX file with:

```bash
.venv/bin/python generate_morph_paper_figures.py \
  --variant normalized-F \
  --output-dir notes/pictures
```

The `notes/pictures` directory is reserved for SVG files included by the TeX
source.
