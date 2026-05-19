#!/usr/bin/env python3
"""Build a standalone HTML viewer for morph frame JSON files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List


def load_frames(input_dir: Path) -> List[Dict[str, object]]:
    frames = []
    for path in sorted(input_dir.glob("frame_*.json")):
        with open(path, "r", encoding="utf-8") as fh:
            frame = json.load(fh)
        if not isinstance(frame, dict):
            raise ValueError(f"{path} does not contain a JSON object")
        for key in ("vertices", "edges"):
            if key not in frame:
                raise ValueError(f"{path} is missing required key {key!r}")
        frame["_name"] = path.name
        frames.append(frame)
    if not frames:
        raise ValueError(f"No frame_*.json files found in {input_dir}")
    return frames


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_dir",
        nargs="?",
        default="boundary_morph_1_100_frames",
        help="Directory containing frame_*.json files.",
    )
    parser.add_argument(
        "--output",
        default="morph_viewer.html",
        help="Output HTML path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output = Path(args.output)
    frames = load_frames(input_dir)

    # Keep the viewer self-contained, but avoid embedding long per-iteration logs.
    compact_frames = []
    for frame in frames:
        compact_frames.append(
            {
                "name": frame.get("_name", ""),
                "vertices": frame["vertices"],
                "edges": frame["edges"],
                "edge_weights": frame.get("edge_weights", []),
                "fixed": frame.get("fixed", []),
                "energy": frame.get("energy"),
                "morph_t": frame.get("morph_t"),
            }
        )

    payload = json.dumps(compact_frames, separators=(",", ":"))
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Morph Frame Viewer</title>
<style>
  :root {{
    color-scheme: light;
    font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  }}
  body {{
    margin: 0;
    background: #f6f7f8;
    color: #14181f;
  }}
  main {{
    max-width: 980px;
    margin: 0 auto;
    padding: 24px;
  }}
  h1 {{
    font-size: 22px;
    font-weight: 650;
    margin: 0 0 16px;
  }}
  .toolbar {{
    display: grid;
    grid-template-columns: auto 1fr auto;
    gap: 12px;
    align-items: center;
    margin-bottom: 14px;
  }}
  button {{
    border: 1px solid #9aa4b2;
    background: #ffffff;
    border-radius: 6px;
    padding: 7px 12px;
    font: inherit;
    cursor: pointer;
  }}
  input[type="range"] {{
    width: 100%;
  }}
  .stats {{
    font-size: 13px;
    color: #4b5563;
    white-space: nowrap;
  }}
  canvas {{
    display: block;
    width: min(100%, 900px);
    aspect-ratio: 1;
    background: #ffffff;
    border: 1px solid #d4d8de;
    border-radius: 8px;
  }}
</style>
</head>
<body>
<main>
  <h1>Morph Frame Viewer</h1>
  <div class="toolbar">
    <button id="play" type="button">Play</button>
    <input id="slider" type="range" min="0" max="{len(compact_frames) - 1}" value="0" step="1">
    <div id="stats" class="stats"></div>
  </div>
  <canvas id="canvas" width="900" height="900"></canvas>
</main>
<script>
const frames = {payload};
const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
const slider = document.getElementById("slider");
const play = document.getElementById("play");
const stats = document.getElementById("stats");
let timer = null;

function project(p) {{
  const margin = 52;
  const r = (canvas.width - 2 * margin) / 2;
  return [
    canvas.width / 2 + p[0] * r,
    canvas.height / 2 - p[1] * r
  ];
}}

function edgeColor(weight, maxWeight) {{
  if (!Number.isFinite(weight) || maxWeight <= 1) return "rgba(55, 65, 81, 0.32)";
  const t = Math.max(0, Math.min(1, (weight - 1) / (maxWeight - 1)));
  const red = Math.round(80 + 175 * t);
  const green = Math.round(92 - 42 * t);
  const blue = Math.round(108 - 68 * t);
  return `rgba(${{red}}, ${{green}}, ${{blue}}, ${{0.28 + 0.52 * t}})`;
}}

function draw(index) {{
  const frame = frames[index];
  const vertices = frame.vertices;
  const edges = frame.edges;
  const weights = frame.edge_weights || [];
  const fixed = new Set(frame.fixed || []);
  const maxWeight = Math.max(1, ...weights);

  ctx.clearRect(0, 0, canvas.width, canvas.height);

  const unitCenter = project([0, 0]);
  const unitEdge = project([1, 0]);
  const unitRadius = Math.abs(unitEdge[0] - unitCenter[0]);
  ctx.beginPath();
  ctx.arc(unitCenter[0], unitCenter[1], unitRadius, 0, Math.PI * 2);
  ctx.fillStyle = "#fbfcfd";
  ctx.fill();
  ctx.strokeStyle = "#1f2937";
  ctx.lineWidth = 2;
  ctx.stroke();

  ctx.lineCap = "round";
  for (let e = 0; e < edges.length; e++) {{
    const [i, j] = edges[e];
    const a = project(vertices[i]);
    const b = project(vertices[j]);
    ctx.beginPath();
    ctx.moveTo(a[0], a[1]);
    ctx.lineTo(b[0], b[1]);
    ctx.strokeStyle = edgeColor(weights[e], maxWeight);
    ctx.lineWidth = 0.8 + 1.8 * Math.min(1, ((weights[e] || 1) - 1) / Math.max(1, maxWeight - 1));
    ctx.stroke();
  }}

  for (let i = 0; i < vertices.length; i++) {{
    const p = project(vertices[i]);
    ctx.beginPath();
    ctx.arc(p[0], p[1], fixed.has(i) ? 4.2 : 2.2, 0, Math.PI * 2);
    ctx.fillStyle = fixed.has(i) ? "#0f766e" : "#111827";
    ctx.fill();
  }}

  const energy = Number(frame.energy);
  const morphT = Number(frame.morph_t);
  stats.textContent = `${{frame.name}}  t=${{Number.isFinite(morphT) ? morphT.toFixed(3) : "n/a"}}  energy=${{Number.isFinite(energy) ? energy.toFixed(6) : "n/a"}}`;
}}

function setFrame(index) {{
  slider.value = String(index);
  draw(index);
}}

slider.addEventListener("input", () => draw(Number(slider.value)));
play.addEventListener("click", () => {{
  if (timer) {{
    clearInterval(timer);
    timer = null;
    play.textContent = "Play";
    return;
  }}
  play.textContent = "Pause";
  timer = setInterval(() => {{
    const next = (Number(slider.value) + 1) % frames.length;
    setFrame(next);
  }}, 180);
}});

draw(0);
</script>
</body>
</html>
"""
    output.write_text(html, encoding="utf-8")
    print(f"Wrote {output} with {len(frames)} frame(s)")


if __name__ == "__main__":
    main()
