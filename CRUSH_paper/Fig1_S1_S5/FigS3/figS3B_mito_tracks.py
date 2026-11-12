#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
MITO = FIG.parent / "reviewer" / "mitotic"
sys.path.insert(0, str(FIG.parent))
import MSD_tracks as msd

RES, CHROM = 25_000, "1"
LABELS = False
SCORES = True
FLAGS = FIG.parent / "reviewer" / "diagnostic" / "mito_flag.json"
A_COL, B_COL, INK, INK2 = "#8B0000", "#000000", "#16181c", "#5a5f68"
SAMPLES = [("WT-90G2-R1", "G2"),
           ("WT-NS-R1", "non-synchronised"),
           ("WT-92Pphase-10m-R1", "prophase 10 min"),
           ("WT-87PMphase-30m-R1", "prometaphase 30 min"),
           ("WT-92PMphase-60m-R1", "prometaphase 60 min")]


def load(sample):
    p = MITO / "crush" / f"{sample}_r{RES}" / f"{sample}_mergedCrush_{RES}.bedgraph"
    if not p.exists():
        return None
    d = msd.load(str(p)).dropna(subset=["Score"])
    c = msd.norm_chrom(d["Chromosome"]).to_numpy()
    d = d[c == CHROM]
    if d.empty:
        return None
    n = int(d["Start"].max()) // RES + 1
    v = np.full(n, np.nan)
    idx = (d["Start"].to_numpy() // RES).astype(int)
    ok = idx < n
    v[idx[ok]] = d["Score"].to_numpy()[ok]
    return v


def main() -> None:
    rows = [(load(s), lab, s) for s, lab in SAMPLES]
    rows = [r for r in rows if r[0] is not None]
    if not rows:
        print("no tracks found"); return
    n = max(len(v) for v, _, _ in rows)
    lim = CRUSH_LIM
    x = np.arange(n) * RES / 1e6
    scores = json.loads(FLAGS.read_text()) if (SCORES and FLAGS.exists()) else {}

    fig, axes = plt.subplots(len(rows), 1, figsize=(11.5, 0.95 * len(rows)),
                             sharex=True)
    fig.patch.set_facecolor("white")
    for ax, (v, lab, s) in zip(axes, rows):
        vv = np.full(n, np.nan); vv[:len(v)] = v
        ax.fill_between(x, 0, np.nan_to_num(vv), where=vv > 0, color=A_COL, lw=0)
        ax.fill_between(x, 0, np.nan_to_num(vv), where=vv <= 0, color=B_COL, lw=0)
        ax.axhline(0, color="#c9c9c4", lw=0.8, zorder=1)
        ax.set_ylim(-lim, lim); ax.set_yticks([])
        if LABELS:
            ax.set_ylabel(lab, fontsize=9.5, color=INK, rotation=0,
                          ha="right", va="center", labelpad=8)
        key = f"{s}|{CHROM}"
        if key in scores:
            ax.text(1.004, 0.5, f"{scores[key]['strength']:.2f}",
                    transform=ax.transAxes, fontsize=9, color=INK,
                    va="center", ha="left")
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color("#d7d7d3")
        ax.tick_params(colors=INK2, labelsize=9.5, length=0)
    if LABELS:
        axes[-1].set_xlabel(f"chr{CHROM} (Mb)", fontsize=10.5, color=INK, labelpad=6)
    fig.subplots_adjust(left=0.13 if LABELS else 0.02,
                        right=0.955 if scores else 0.99,
                        top=0.97, bottom=0.11, hspace=0.16)
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "figS3B_mito_tracks.png"
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, facecolor="white")
    if scores:
        import csv
        q = OUT / "figS3B_mito_tracks_sourcedata.csv"
        shown = {s_ for s_, _ in SAMPLES}
        with q.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sample", "stage", "plotted", "strength_AABB_2AB",
                        "null_mean", "null_sd", "z", "scored_at_bp", "norm"])
            label = dict(SAMPLES)
            for key, v in scores.items():
                samp = key.split("|", 1)[0]
                w.writerow([samp, label.get(samp, "replicate, not plotted"),
                            "yes" if samp in shown else "no",
                            f"{v['strength']:.6f}", f"{v['null']:.6f}",
                            f"{v.get('null_sd', float('nan')):.6f}",
                            f"{v['z']:.3f}", v.get("scale", ""), v.get("norm", "")])
        print(f"wrote {q}")
    print(f"wrote {out}")
    for v, lab, s in rows:
        cov = np.isfinite(v)
        sc = scores.get(f"{s}|{CHROM}", {}).get("strength", float("nan"))
        print(f"  {lab:<28}{int(cov.sum()):>7,} bins  {100*np.mean(v[cov]>0):>5.1f}% A"
              f"   score {sc:5.3f}")


if __name__ == "__main__":
    main()
