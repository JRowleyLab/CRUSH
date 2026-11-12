#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
DATA = FIG.parent / "reviewer" / "validation"
MAXSEP = 40_000_000
RED, BLUE, INK, INK2 = "#C6303A", "#2E77B5", "#16181c", "#5a5f68"


def main() -> None:
    a = np.load(DATA / "strength_vs_dist_chr20.npy")
    a = a[a[:, 0] <= MAXSEP]
    x, cr, ei = a[:, 0] / 1e6, a[:, 1], a[:, 2]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    fig.patch.set_facecolor("white")
    ax.plot(x, cr, color=RED, linewidth=5.0, solid_capstyle="round",
            clip_on=False, zorder=3)
    ax.axhline(1.0, color="#8a8f98", linewidth=1.6, linestyle=(0, (5, 3)), zorder=1)
    ax.text(x.max(), 1.0, "  no compartment preference", fontsize=9.5,
            color=INK2, va="center")
    ax.set_xlim(0, x.max())
    ax.set_ylim(0, cr.max() * 1.18)
    ax.set_xlabel("minimum separation considered (Mb)", fontsize=10.5,
                  color=INK, labelpad=8)
    ax.set_ylabel("compartment strength  (AA+BB)/2AB", fontsize=10.5,
                  color=INK, labelpad=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d7d7d3")
    ax.tick_params(colors=INK2, labelsize=10, length=0)
    ax.grid(True, color="#ececE8", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.text(0.0, 1.04, "chr20, 1 kb, whole chromosome", transform=ax.transAxes,
            fontsize=9.5, color=INK2)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "figS2B_strength_vs_dist.png"
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, facecolor="white")
    print(f"wrote {out} (+ .pdf)")
    import csv
    q = OUT / "figS2B_strength_vs_dist_sourcedata.csv"
    with q.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["min_separation_Mb", "CRUSH_strength", "eigen_strength"])
        for xi, c, e in zip(x, cr, ei):
            w.writerow([f"{xi:g}", f"{c:.6f}", f"{e:.6f}"])
    print(f"wrote {q}")
    print(f"  CRUSH       {cr[0]:.2f} at 0 Mb -> {cr[-1]:.2f} at {x[-1]:.0f} Mb"
          f"   (min {cr.min():.2f}, max {cr.max():.2f})")
    print(f"  (eigenvector not plotted; {ei[0]:.2f} -> {ei[-1]:.2f}, "
          f"retained in the .npy)")


if __name__ == "__main__":
    main()
