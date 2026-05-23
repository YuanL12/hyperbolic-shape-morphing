# Graph Morphing Problem on Hyperbolic Surface
Theorem: Tutte Embedding works for nonsymmetric weights on non-triangulation(just planar graph). 


## Workflow

```bash

# generate several frames as .json files
python make_boundary_weight_morph_mean_value.py --base input/example_genus2_boundary_vertices.json  --output-dir morph_frames --frame 25

# visualzie it as a html file
python visualize_morph_frames.py morph_frames --output morph_viewer.html
```


## Experiment Result 
1. `morph_viewer_genus_2_regular.html`
``` bash 
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --base input/example_genus2_boundary_vertices_embedded_1_20.json \
  --output-dir morph_frames

.venv/bin/python visualize_morph_frames.py \
  morph_frames --output morph_viewer_genus_2_regular.html
```

2. `morph_viewer_genus_2_irregular.html`
``` bash 
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --base input/example_genus2_boundary_vertices_embedded_irregular_1_20.json \
  --output-dir boundary_morph_irregular_1_20_frames

.venv/bin/python visualize_morph_frames.py \
  boundary_morph_irregular_1_20_frames \
  --output morph_viewer_genus_2_irregular.html
```

3. `morph_viewer_genus_3_1`
- Base:  `input/example_genus3_klein_quartic_embedded_irregular.json`
- Initial shape: center vertex with high degree
- weight: mean-value coordinates -> half-20-half-1


4. `morph_viewer_genus_3_2`
- Base:  `input/example_genus3_klein_quartic_embedded_irregular.json`
- Initial shape: random points + Delaunay triangulation
- weight: mean-value coordinates -> half-20-half-1 

Generator:
- `generate_genus3_klein_embedding.py`

5. `morph_viewer_genus3_quad_pent_1_to_half_10.html`
`input/example_genus3_klein_quartic_quad_pent_graph.json`: planar graph containing quads and pentagons, but not convex
`input/example_genus3_klein_quartic_quad_pent_graph_solution.json` planar graph containing quads and pentagons, and convex. Generated from the above one by using `poincare_harmonic_map` to minizize energy with weights all one

Weights: all ones to half-10-half-one

6. `morph_viewer_genus3_quad_pent_mean_value_frame0_to_frame24.html`
Weights: mean-valued weights to mean-valued weights


ok, now we are close to the final goal. Use our previous chat. 
For Genus-2 or Genus-3:
Each of them should contain two morphing using mean-value coordinates for weights
1. triangulation (using random samling points then and Delaunay triangulation used above).
2. cells (containing quads and pentagons by removing edges like above). Since mean-value coordinates needs convex embedding, i.e. every cell should be convex, we need to minimize its Dirichlet Energy to make it convex. 

In other words, the input I want is two convex initial embedding, then use mean-value weights to generate the morphing. To generate the two initial embeddings, as we have discussed previously, one can use all ones edge weight in [poincare_harmonic_map.py](poincare_harmonic_map.py) and the other can be half-one-half-10.  

Note: 
1. Another way is to add edges to make its a valid triangulation. The chocie of those two methods depend on the purpose/application, sometimes we want to keep the cells. 

2. When computing the Ensure that computing mean-value directed weights always can see the attached vertices when checking the valence. 

- 

## Two different stages/workflows.
### Harmonic map from a given edge weight
[poincare_harmonic_map.py](/Users/yluo/Documents/Codes/shape-morphing/poincare_harmonic_map.py:1) is the general solver.

It takes one problem JSON:

```text
vertices + edges/faces + weights + constraints + fixed vertices
```

and solves one harmonic map minimization:

```text
f = argmin E(f)
```

Example:

```bash
.venv/bin/python poincare_harmonic_map.py \
  example_genus2_boundary_vertices_weighted_1_100.json \
  --output example_genus2_boundary_vertices_weighted_1_100_solution.json
```

That computes the weighted solution `f_0`.

### Morphing by interpolating two harmonic maps
`make_boundary_weight_morph.py` is the morph/continuation workflow.

It assumes you already have:

```text
example_genus2_boundary_vertices_weighted_1_100_solution.json
```

Then it gradually changes the weights:

```text
1/100-weighted case -> all weights 1
```

and solves many harmonic map problems, one per frame:

```text
f_0, f_1, f_2, ..., f_24
```


## Generate html view
```bash
.venv/bin/python visualize_morph_frames.py boundary_morph_1_20_frames --output morph_viewer_20.html
```

## Mean value coordinates
An embedding (realized as 2D coordinates) gives mean value coordinates(edge weights). 


## Term
1. Poincare disk 
2. Hyperbolic Disck (same as the above)
3. Bolza Surface (regular 8-gons)
