# Irregular Genus-3 Klein Quartic Embedding

This note summarizes the current method for generating the less-regular
genus-3 initial embedding used before any morphing.

## Goal

Generate a Klein-quartic-style genus-3 cut surface in the Poincare disk with:

- a regular 14-gon fundamental domain,
- side identifications matching the standard Klein quartic picture,
- fixed 14 corner vertices,
- no high-degree center vertex,
- an irregular, non-ring triangulation,
- and a validated planar disk embedding before running morphing.

The generator is:

```bash
.venv/bin/python generate_genus3_klein_embedding.py
```

Default outputs:

```text
input/example_genus3_klein_quartic_embedded_irregular.json
output/genus3_klein_initial_embedding.svg
```

## Side Pairing Model

The fundamental polygon uses 14 sides. The side labels are ordered as:

```text
[1, 7, 3, 2, 5, 4, 7, 6, 2, 1, 4, 3, 6, 5]
```

This follows the standard Klein quartic 14-gon figure: sides with equal labels
are paired with reversed boundary orientation.

For each paired side, the script computes a Poincare disk Mobius isometry
mapping the master side to the slave side. Boundary interior vertices on slave
sides are constrained to the corresponding reversed-order master-side vertices.
The 14 polygon corners are fixed.

## Triangulation

The previous version used concentric rings and a center fan, which was too
regular and created a high-degree central vertex.

The current version avoids that by:

1. Sampling boundary vertices on the 14-gon sides.
2. Sampling random interior points inside the Euclidean 14-gon.
3. Rejecting points that are too close to existing samples.
4. Running a local Bowyer-Watson Delaunay triangulation.
5. Keeping triangles whose centroids lie inside the 14-gon.

This produces a less symmetric triangulation and removes the center vertex.

Current default statistics:

```text
vertices=250
faces=428
edges=677
constraints=28
fixed corners=14
max vertex degree=10
average degree=5.416
vertices with r < 0.05: 0
```

## Relaxation and Irregular Embedding

The raw triangulation plus side constraints is not used directly. The script
first solves a constrained harmonic map problem using `HarmonicMapSolver`.
This relaxes the constrained cut mesh into a valid disk embedding.

After relaxation, the generator applies a small smooth perturbation to root
vertices only. Slave vertices are then resolved through the existing side
pairing constraints, so the paired boundary data remains consistent.

## Validation

The generated embedding is checked with `validate_disk_embedding`.

The validation confirms:

```text
outside disk vertices=0
positive faces=428
negative faces=0
degenerate faces=0
Euclidean chord crossings=0
Poincare geodesic crossings=0
```

The output metadata also records how many vertices lie outside the regular
14-gon fundamental domain while still staying inside the Poincare disk. The
current default has:

```text
outside regular 14-gon domain vertices=29
```

## Preview

The script writes an SVG preview before morphing:

```text
output/genus3_klein_initial_embedding.svg
```

The preview draws:

- green Poincare geodesic mesh edges,
- dashed black fundamental-domain sides,
- red fixed corner vertices,
- black free vertices.

This preview should be inspected before running the morph script.

## Morphing Command

After accepting the initial embedding, frames can be generated with:

```bash
.venv/bin/python make_boundary_weight_morph_mean_value.py \
  --base input/example_genus3_klein_quartic_embedded_irregular.json \
  --output-dir boundary_morph_genus3_klein_irregular_frames
```

The morphing step was intentionally not run while iterating on the initial embedding preview.
