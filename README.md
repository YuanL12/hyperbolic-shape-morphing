# Hyperbolic Shape Morphing in the Poincaré Disk

`hyper_morph` provides a
frame generator for morphing between two embeddings. 

It also supports 
-  computing discrete harmonic maps and directed hyperbolic
- mean-value weights for graphs embedded in the Poincaré disk. 


## Example

The same 50-frame genus-2 morph is shown at three tiling depths.

| Depth 0 | Depth 1 | Depth 2 |
| --- | --- | --- |
| ![Morph at tile depth 0](examples/output/checkerboard_morph/Z0_to_Z1_sin8_1_50_ring_alternating/sin8_ring_alternating_depth0.gif) | ![Morph at tile depth 1](examples/output/checkerboard_morph/Z0_to_Z1_sin8_1_50_ring_alternating/sin8_ring_alternating_depth1.gif) | ![Morph at tile depth 2](examples/output/checkerboard_morph/Z0_to_Z1_sin8_1_50_ring_alternating/sin8_ring_alternating_depth2.gif) |

MP4 versions: [depth 0](examples/output/checkerboard_morph/Z0_to_Z1_sin8_1_50_ring_alternating/sin8_ring_alternating_depth0.mp4),
[depth 1](examples/output/checkerboard_morph/Z0_to_Z1_sin8_1_50_ring_alternating/sin8_ring_alternating_depth1.mp4),
[depth 2](examples/output/checkerboard_morph/Z0_to_Z1_sin8_1_50_ring_alternating/sin8_ring_alternating_depth2.mp4).

## Quick start
We use uv to manage the environment and dependencies, but it is not required.

Install the locked environment:

```bash
uv sync
```

Generate a small input mesh, then solve it:

```bash
uv run python examples/generate_input.py
uv run hyper-morph-solve \
  examples/input/simple_disk.json \
  --output examples/output/simple_disk_solution.json
```

The generator writes a five-vertex disk mesh with four fixed boundary vertices and one off-center interior vertex.

An input JSON file needs `vertices` plus either `faces` or `edges`:

```json
{
  "vertices": [[0.18, 0.12], [-0.45, -0.45], [0.45, -0.45]],
  "edges": [[0, 1], [1, 2], [2, 0]],
  "fixed": [1, 2],
  "iterations": 300,
  "step_size": 0.05,
  "tolerance": 1e-10
}
```

All vertices are `[x, y]` points strictly inside the unit disk. `fixed`,
`iterations`, `step_size`, and `tolerance` are optional. Cut meshes can also
provide Möbius `constraints`; see the generators in `examples/` for complete
surface examples.

## Python API

```python
import json

from hyper_morph import HarmonicMapSolver

with open("examples/input/simple_disk.json", encoding="utf-8") as file:
    mesh = json.load(file)

result = HarmonicMapSolver(mesh).solve()
print(result["vertices"])
```

Compute normalized directed mean-value weights:

```bash
uv run hyper-morph-weights \
  examples/input/simple_disk.json \
  --normalization normalized \
  --output examples/output/simple_disk_weights.json
```

Generate a multi-frame morph between compatible source and target embeddings:

```bash
uv run hyper-morph \
  --source-embedding path/to/source.json \
  --target-embedding path/to/target.json \
  --target-directed-weights mean_value \
  --output-dir examples/output/frames/demo
```

Create a self-contained HTML viewer for those frames:

```bash
uv run python utils/visualize_morph_frames.py \
  examples/output/frames/demo \
  --output examples/output_html/demo.html
```
