#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TRACKS = {"CRUSH v2 1 kb": "#b2182b", "eigenvector 1 kb": "#2166ac", "eigenvector 25 kb": "#67a9cf"}

r = json.loads((HERE / "domain_classes.json").read_text())["tracks"]
fig, ax = plt.subplots(figsize=(2.8, 3.0))
vals = [r[t]["loop_domains_with_internal_switch"] * 100 for t in TRACKS]
ax.bar(range(len(TRACKS)), vals, color=list(TRACKS.values()), width=0.65)
for i, v in enumerate(vals):
    ax.text(i, v + 0.5, f"{v:.1f}", ha="center", va="bottom", fontsize=7)
ax.set_xticks(range(len(TRACKS)))
ax.set_xticklabels([t.replace(" ", "\n", 1) for t in TRACKS], fontsize=7)
ax.set_ylabel("corner-loop domains with an\ninternal compartment switch (%)", fontsize=8)
ax.tick_params(labelsize=7)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(HERE / f"loop_subdivision.{ext}", dpi=300, bbox_inches="tight")
