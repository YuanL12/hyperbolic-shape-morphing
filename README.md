# Graph Morphing Problem on Hyperbolic Surface

## Theorem
Tutte Embedding works for nonsymmetric weights on non-triangulation(just planar graph). 


## Current Experiment
We consider 
- genus 2: regular octagonal fundamental polygon for the Bolza surface
- genus 3: regular 14-gonal fundamental polygon for the Klein quartic surface


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
- Directed mean-value weights are normalized at each vertex by default. Pass
  `--unnormalized` to use raw outgoing mean-value weights instead.
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
  --output-dir output/frames/boundary_morph_genus2_bolza_triangulation_mean_value_fixed_corner_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check \
  --reference-fundamental-domain

.venv/bin/python visualize_morph_frames.py \
  output/frames/boundary_morph_genus2_bolza_triangulation_mean_value_fixed_corner_frames \
  --output output_html/morph_viewer_genus2_bolza_triangulation_mean_value_fixed_corner.html
```

### Genus-2 Cells

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --source-embedding input/example_genus2_bolza_convex_quad_pent_graph_all_1_solution.json \
  --target-embedding input/example_genus2_bolza_convex_quad_pent_graph_half_10_solution.json \
  --output-dir output/frames/boundary_morph_genus2_bolza_cells_mean_value_fixed_corner_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check \
  --reference-fundamental-domain

.venv/bin/python visualize_morph_frames.py \
  output/frames/boundary_morph_genus2_bolza_cells_mean_value_fixed_corner_frames \
  --output output_html/morph_viewer_genus2_bolza_cells_mean_value_fixed_corner.html
```

### Genus-3 Triangulation

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --source-embedding input/example_genus3_klein_quartic_triangulation_all_1_solution.json \
  --target-embedding input/example_genus3_klein_quartic_triangulation_half_10_solution.json \
  --output-dir output/frames/boundary_morph_genus3_triangulation_mean_value_fixed_corner_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check \
  --reference-fundamental-domain

.venv/bin/python visualize_morph_frames.py \
  output/frames/boundary_morph_genus3_triangulation_mean_value_fixed_corner_frames \
  --output output_html/morph_viewer_genus3_triangulation_mean_value_fixed_corner.html
```

### Genus-3 Cells

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --source-embedding input/example_genus3_klein_quartic_convex_quad_pent_graph_all_1_solution.json \
  --target-embedding input/example_genus3_klein_quartic_convex_quad_pent_graph_half_10_solution.json \
  --output-dir output/frames/boundary_morph_genus3_cells_mean_value_fixed_corner_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check \
  --reference-fundamental-domain

.venv/bin/python visualize_morph_frames.py \
  output/frames/boundary_morph_genus3_cells_mean_value_fixed_corner_frames \
  --output output_html/morph_viewer_genus3_cells_mean_value_fixed_corner.html
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

After the four frame directories exist, regenerate the PDF grids used by the
Overleaf paper with:

```bash
.venv/bin/python generate_morph_paper_figures.py \
  --variant mean-value-fixed_corner \
  --output-dir paper/shape-morphing-overleaf/pictures
```

## Relax Corner Orbits

To try the same four normalized mean-value morphs with the corner orbits relaxed, first
make copied endpoint inputs. This does not modify the original `input/*.json`
files used above. Passing corner positions `0 1` relaxes the one Bolza corner
orbit in genus 2 and both Klein-quartic corner orbits in genus 3.

```bash
.venv/bin/python make_relaxed_corner_input_copies.py \
  --corner-positions 0 1 \
  --suffix _relaxed_corner_orbits
```

The copied endpoint files remove the selected corner orbit(s) from `fixed` and
add the corresponding corner Mobius constraints, so the selected corners are
free and their paired/slave corner copies follow by transformation. The original
reference polygon is preserved in `reference_fundamental_domain`, a drawing-only
JSON object with `corner_indices` and `vertices`. Thus `fixed` means solver-fixed
vertices only, while `reference_fundamental_domain` controls the fixed black
fundamental polygon used as visual context.

Here "fundamental domain" means the cut-open polygon in the Poincare disk. The
genus-2 examples use a regular octagonal fundamental polygon for the Bolza
surface. The genus-3 examples use a regular 14-gonal fundamental polygon for the
Klein quartic surface; the Klein quartic is the surface, not the polygon itself.

Then rerun the four morph commands with the copied inputs, suffixed output
directories, and `--reference-fundamental-domain`. For example, the genus-3 cell
variant is:

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --source-embedding input/example_genus3_klein_quartic_convex_quad_pent_graph_all_1_relaxed_corner_orbits_solution.json \
  --target-embedding input/example_genus3_klein_quartic_convex_quad_pent_graph_half_10_relaxed_corner_orbits_solution.json \
  --output-dir output/frames/boundary_morph_genus3_cells_mean_value_relaxed_corner_frames \
  --frames 25 \
  --start-directed-weights mean_value \
  --target-directed-weights mean_value \
  --iterations 3000 \
  --step-size 0.001 \
  --tolerance 1e-8 \
  --line-search-objective none \
  --convergence-criterion relative_step \
  --skip-face-orientation-check \
  --reference-fundamental-domain

.venv/bin/python visualize_morph_frames.py \
  output/frames/boundary_morph_genus3_cells_mean_value_relaxed_corner_frames \
  --output output_html/morph_viewer_genus3_cells_mean_value_relaxed_corner.html
```

The PDF grids for all four relaxed-corner-orbit frame directories can be regenerated
with:

```bash
.venv/bin/python generate_morph_paper_figures.py \
  --variant mean-value-relaxed_corner \
  --output-dir paper/shape-morphing-overleaf/pictures
```
