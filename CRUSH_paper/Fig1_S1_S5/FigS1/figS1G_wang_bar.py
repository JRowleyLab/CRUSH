#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"

SWEEP = OUT / "figS1G_wang_sweep_prom.json"
if not SWEEP.exists():
    import subprocess
    for step in ("figS1G_prom_termlists.R", "figS1G_wang_sweep_prom.R"):
        need = OUT / ("top200_goids_prom.json" if "termlists" in step
                      else "figS1G_wang_sweep_prom.json")
        if need.exists():
            continue
        print(f"{need.name} not found -- running {step}", flush=True)
        subprocess.run(["Rscript", str(HERE / step)], cwd=str(FIG), check=True)
    print()
D = json.loads(SWEEP.read_text())
INK, INK2 = "#16181c", "#5a5f68"
LABELS = False
NS = [50, 100]
SERIES = [("crush_reps", "CRUSH, rep 1 vs rep 2", "#8B0000"),
          ("cross_rep1", "CRUSH vs eigenvector, rep 1", "#2E77B5"),
          ("cross_rep2", "CRUSH vs eigenvector, rep 2", "#7FB0DC")]

fig, ax = plt.subplots(figsize=(6.4, 4.2))
fig.patch.set_facecolor("white"); ax.set_facecolor("white")

ax.axhspan(-3, 2, color="#ececea", lw=0, zorder=0)
ax.axhline(2, color=INK2, lw=1.0, linestyle=(0, (4, 2.5)), zorder=1)

w = 0.26
x = np.arange(len(NS))
for j, (key, lab, col) in enumerate(SERIES):
    ys = [D[f"{n}|{key}"]["z"] for n in NS]
    ws = [D[f"{n}|{key}"]["wang"] for n in NS]
    xs = x + (j - 1) * w
    ax.bar(xs, ys, w * 0.88, color=col, lw=0, label=lab, zorder=3)
    for xi, yi in zip(xs, ys):
        if abs(yi) < 0.15:
            ax.plot([xi - w * 0.44, xi + w * 0.44], [0, 0], color=col,
                    lw=3.0, solid_capstyle="butt", zorder=4)
    if LABELS:
        for xi, yi, wi in zip(xs, ys, ws):
            ax.text(xi, yi + 0.35, f"{wi:.2f}", ha="center", va="bottom",
                    fontsize=8, color=INK2, zorder=4)

ax.set_xticks(x)
ax.set_xticklabels([f"top {n}" for n in NS], fontsize=10.5, color=INK)
ax.set_xlim(-0.55, len(NS) - 0.45)
ax.set_ylim(-3, 17)
if LABELS:
    ax.set_ylabel("Wang similarity, z vs a null matched to N",
                  fontsize=10, color=INK, labelpad=8)
    ax.legend(frameon=False, fontsize=9, loc="upper left", handlelength=1.3)
ax.tick_params(colors=INK2, labelsize=9.5, length=0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color("#d7d7d3")
ax.set_axisbelow(True)
ax.yaxis.grid(True, color="#f3f3f1", lw=0.8)

fig.tight_layout()
out = OUT / "figS1G_wang_bar.png"
for ext in ("png", "pdf"):
    fig.savefig(out.with_suffix(f".{ext}"), dpi=300, facecolor="white")
print(f"wrote {out} (+ .pdf)")

import csv
rows = []
for n in NS:
    for key, lab, _ in SERIES:
        d = D[f"{n}|{key}"]
        rows.append({"top_N": n, "comparison": lab, "wang": d["wang"],
                     "null_mean": d["null_mean"], "null_sd": d["null_sd"],
                     "z": d["z"], "null_q95": d["null_q95"]})
csvp = OUT / "figS1G_wang_bar_sourcedata.csv"
with csvp.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print(f"wrote {csvp}")
for n in NS:
    for key, lab, _ in SERIES:
        d = D[f"{n}|{key}"]
        print(f"  top {n:<4}{lab:<30}Wang {d['wang']:.3f}   z {d['z']:+5.1f}")
