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
        default="examples/output/frames/morph_frames",
        help="Directory containing frame_*.json files.",
    )
    parser.add_argument(
        "--output",
        default="examples/output_html/morph_viewer.html",
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
                "reference_fundamental_domain": frame.get(
                    "reference_fundamental_domain"
                ),
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

function norm2(p) {{
  return p[0] * p[0] + p[1] * p[1];
}}

function normalizeAngle(delta) {{
  while (delta <= -Math.PI) delta += Math.PI * 2;
  while (delta > Math.PI) delta -= Math.PI * 2;
  return delta;
}}

function geodesicPoints(a, b) {{
  const det = a[0] * b[1] - a[1] * b[0];
  const samples = 36;

  if (Math.abs(det) < 1e-10) {{
    return [a, b];
  }}

  // Poincare disk geodesics are circles orthogonal to the unit circle.
  // Their Euclidean center c satisfies 2 c.p = |p|^2 + 1 for both endpoints.
  const rhsA = (norm2(a) + 1) * 0.5;
  const rhsB = (norm2(b) + 1) * 0.5;
  const center = [
    (rhsA * b[1] - rhsB * a[1]) / det,
    (a[0] * rhsB - b[0] * rhsA) / det
  ];
  const radius = Math.hypot(a[0] - center[0], a[1] - center[1]);
  if (!Number.isFinite(radius) || radius <= 1e-10) {{
    return [a, b];
  }}

  const start = Math.atan2(a[1] - center[1], a[0] - center[0]);
  const end = Math.atan2(b[1] - center[1], b[0] - center[0]);
  let delta = normalizeAngle(end - start);

  const midShort = [
    center[0] + radius * Math.cos(start + delta * 0.5),
    center[1] + radius * Math.sin(start + delta * 0.5)
  ];
  if (norm2(midShort) > 1 + 1e-7) {{
    delta += delta > 0 ? -Math.PI * 2 : Math.PI * 2;
  }}

  const points = [];
  for (let k = 0; k <= samples; k++) {{
    const t = k / samples;
    const angle = start + delta * t;
    points.push([
      center[0] + radius * Math.cos(angle),
      center[1] + radius * Math.sin(angle)
    ]);
  }}
  return points;
}}

function drawFundamentalDomain(domainVertices) {{
  if (!domainVertices || domainVertices.length < 3) {{
    return;
  }}

  ctx.save();
  ctx.setLineDash([]);
  ctx.lineCap = "round";
  ctx.lineJoin = "round";
  ctx.globalAlpha = 0.62;
  ctx.strokeStyle = "#111827";
  ctx.lineWidth = 1.2;

  for (let k = 0; k < domainVertices.length; k++) {{
    const a = domainVertices[k];
    const b = domainVertices[(k + 1) % domainVertices.length];
    const points = geodesicPoints(a, b).map(project);
    ctx.beginPath();
    ctx.moveTo(points[0][0], points[0][1]);
    for (let p = 1; p < points.length; p++) {{
      ctx.lineTo(points[p][0], points[p][1]);
    }}
    ctx.stroke();
  }}

  ctx.restore();
}}

function draw(index) {{
  const frame = frames[index];
  const vertices = frame.vertices;
  const edges = frame.edges;
  const referenceDomain = frame.reference_fundamental_domain || {{}};
  const referenceCornerIndices = referenceDomain.corner_indices || [];
  const domainVertices = referenceDomain.vertices || [];
  const boundaryCount = referenceCornerIndices.length >= 2
    ? referenceCornerIndices.length * (referenceCornerIndices[1] - referenceCornerIndices[0])
    : 0;

  ctx.clearRect(0, 0, canvas.width, canvas.height);

  const unitCenter = project([0, 0]);
  const unitEdge = project([1, 0]);
  const unitRadius = Math.abs(unitEdge[0] - unitCenter[0]);
  ctx.beginPath();
  ctx.arc(unitCenter[0], unitCenter[1], unitRadius, 0, Math.PI * 2);
  ctx.fillStyle = "#fbfcfd";
  ctx.fill();
  ctx.strokeStyle = "#1f2937";
  ctx.globalAlpha = 0.8;
  ctx.lineWidth = 1.2;
  ctx.stroke();
  ctx.globalAlpha = 1;

  drawFundamentalDomain(domainVertices);

  ctx.lineCap = "round";
  ctx.lineJoin = "round";
  ctx.strokeStyle = "#138a2e";
  ctx.globalAlpha = 0.85;
  ctx.lineWidth = 1.5;
  for (let e = 0; e < edges.length; e++) {{
    const [i, j] = edges[e];
    const points = geodesicPoints(vertices[i], vertices[j]).map(project);
    ctx.beginPath();
    ctx.moveTo(points[0][0], points[0][1]);
    for (let k = 1; k < points.length; k++) {{
      ctx.lineTo(points[k][0], points[k][1]);
    }}
    ctx.stroke();
  }}
  ctx.globalAlpha = 1;

  for (let i = 0; i < vertices.length; i++) {{
    const p = project(vertices[i]);
    ctx.beginPath();
    ctx.arc(p[0], p[1], 2.2, 0, Math.PI * 2);
    ctx.fillStyle = i < boundaryCount ? "#dc2626" : "#111827";
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
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    print(f"Wrote {output} with {len(frames)} frame(s)")


if __name__ == "__main__":
    main()
