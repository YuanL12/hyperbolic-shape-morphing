# Sin(8 Theta) Checkerboard Morph

This note records the smooth angular edge-weight experiment for the symmetric
genus-2 center-fan triangulation.

## Base Embedding

The source embedding is
`input/Z0_genus2_bolza_symmetric_center_fan.json`.

It is the symmetric genus-2 Bolza triangulation with:

- one center vertex,
- three concentric rings,
- `40` sectors on each ring,
- `200` triangular faces,
- boundary vertices placed on the Poincare geodesic sides of the regular
  Bolza octagon.

This source embedding is used as `Z^0`.

## Smooth Endpoint Weights

For the endpoint input, each undirected edge receives a weight determined by
the angular position of its Euclidean midpoint in the Poincare disk.

For edge $ij$, let

$$\theta_{ij} = \operatorname{atan2}((y_i + y_j) / 2, (x_i + x_j) / 2).$$

The scalar edge weight is

$$
w_{ij} = w_{\min} + (w_{\max} - w_{\min}) * (1/2 + 1/2 * \sin(8 \theta_{ij})).
$$

For the current run:

```text
w_min = 1
w_max = 50
frequency = 8
phase = 0
```

The realized finite-edge range is approximately:

```text
min weight = 1.0466196024372703
max weight = 49.953380397562725
```

The endpoint input was generated with:

```bash
python3 make_sin8_endpoint_input.py \
  input/Z0_genus2_bolza_symmetric_center_fan.json \
  --output input/Z1_genus2_bolza_symmetric_center_fan_sin8_1_50_input.json \
  --min-weight 1 \
  --max-weight 50 \
  --frequency 8
```

This produces symmetric directed weights by storing the same scalar value in
both edge directions:

```text
directed_edge_weights[e] = [w_e, w_e]
```

## Harmonic Endpoint

The endpoint embedding `Z^1` is obtained by minimizing the Poincare harmonic
energy with the sin(8 theta) edge weights, keeping the quotient-side
constraints and the reference fundamental domain metadata.

Command:

```bash
python3 poincare_harmonic_map.py \
  input/Z1_genus2_bolza_symmetric_center_fan_sin8_1_50_input.json \
  --output input/Z1_genus2_bolza_symmetric_center_fan_sin8_1_50_solution.json \
  --iterations 20000 \
  --step-size 0.004 \
  --tolerance 1e-10
```

Observed solve summary:

```text
Final energy: 1046.383129688146
Mean free gradient norm: 8.480458026992e-07
Recorded iterations: 145
```

After solving, the output metadata was updated to preserve the
`reference_fundamental_domain` from the input.

## Mean-Value Morph Frames

The visible morph uses mean-value directed edge weights for both endpoints.
The source endpoint is `Z^0`; the target endpoint is the sin(8 theta)
harmonic solution `Z^1`.

Command:

```bash
python3 make_boundary_weight_morph_mean_value.py \
  --source-embedding input/Z0_genus2_bolza_symmetric_center_fan.json \
  --target-directed-weights mean_value \
  --target-embedding input/Z1_genus2_bolza_symmetric_center_fan_sin8_1_50_solution.json \
  --output-dir output/frames/Z0_to_Z1_sin8_1_50_symmetric_mean_value_frames_50 \
  --frames 50 \
  --iterations 5000 \
  --step-size 0.01 \
  --tolerance 1e-9 \
  --reference-fundamental-domain
```

Validation checked all `50` morph frames plus the endpoint:

```text
checked embeddings = 51
min_abs_face_area2 = 0.004228885437926522
max_crossings = 0
bad_count = 0
```

## Ring-Alternating Checkerboard Rendering

The black-white art uses the face IDs of `Z^0` as the color reference, so the
checkerboard regions move with the same faces during the morph.

Each depth was rendered with:

```bash
python3 generate_poincare_checkerboard.py \
  --input output/frames/Z0_to_Z1_sin8_1_50_symmetric_mean_value_frames_50/frame_000.json \
  --color-reference input/Z0_genus2_bolza_symmetric_center_fan.json \
  --pattern ring_alternating \
  --tile-depth 2 \
  --output output/checkerboard_morph/Z0_to_Z1_sin8_1_50_ring_alternating/depth2/frame_000.svg
```

The actual run repeats this command for all frames `000` through `049` and for
tile depths `0`, `1`, and `2`.

The HTML viewer is:

```text
output_html/checkerboard_morph_Z0_Z1_sin8_1_50_ring_alternating.html
```

## Boundary Gaps

The gray gaps near the disk boundary in finite-depth renderings are expected.
They are not holes in the quotient surface. They occur because only finitely
many side-pairing copies are drawn. In the infinite tiling, deeper copies fill
those regions.

Small hairline seams can also appear from SVG antialiasing and clipping; those
are rendering artifacts, not embedding validity failures.
