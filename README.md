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
[make_boundary_weight_1_100_morph.py](/Users/yluo/Documents/Codes/shape-morphing/make_boundary_weight_1_100_morph.py:1) is the morph/continuation workflow.

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

So the relationship is:

```text
poincare_harmonic_map.py
  solves one problem
  can generate the starting solution JSON

make_boundary_weight_1_100_morph.py
  solves a sequence of problems
  uses the starting solution JSON as frame 0's initial condition
```

In pipeline form:

```text
example_genus2_boundary_vertices_weighted_1_100.json
        ↓ poincare_harmonic_map.py
example_genus2_boundary_vertices_weighted_1_100_solution.json
        ↓ make_boundary_weight_1_100_morph.py
boundary_morph_1_100_frames/frame_000.json
boundary_morph_1_100_frames/frame_001.json
...
```