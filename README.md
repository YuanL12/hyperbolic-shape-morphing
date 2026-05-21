## Workflow

```bash

# generate several frames as .json files
python make_boundary_weight_morph_mean_value.py --base input/example_genus2_boundary_vertices.json  --output-dir morph_frames --frame 25

# visualzie it as a html file
python visualize_morph_frames.py morph_frames --output morph_viewer.html
```



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
