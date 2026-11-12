#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
XDATA = FIG.parent / "reviewer" / "inactive_x"
sys.path.insert(0, str(FIG.parent))
import MSD_tracks as msd

RES, CRES, ERES = 5000, 500_000, 500_000
EIG = Path("GSE99991")
ARMS = [("WT_Xa", "WT active X", "WT1.merged.cas"),
        ("WT_Xi", "WT inactive X", "WT1.merged.mus")]
DXZ4 = 72_500_000
A_COL, B_COL, INK, INK2 = "#8B0000", "#000000", "#16181c", "#5a5f68"
E_POS, E_NEG = "#14532D", "#12386B"
LABELS = False
SCORES = True
CRUSH_LIM = 200.0
FLAGS = FIG.parent / "reviewer" / "diagnostic" / "x_flag.json"


GENES_BED = Path("mm9.refGene_nochr.bed")


def gene_density(starts, res):
    """Genes per bin on chrX, for orienting a track by the CRUSH convention."""
    if not GENES_BED.exists():
        return None
    g = pd.read_csv(GENES_BED, sep="\t", header=None, usecols=[0, 1],
                    names=["c", "s"], dtype={"c": str})
    g = g[g["c"].str.removeprefix("chr") == "X"]
    idx = (g["s"].to_numpy() // res).astype(int)
    n = int(starts.max()) // res + 2
    d = np.bincount(idx[idx < n], minlength=n).astype(float)
    return d[(starts // res).astype(int)]


def load_eigen(stem, n_target):
    """Eigenvector at 500 kb, oriented to gene density, expanded for display."""
    p = EIG / f"{stem}.eigen_{ERES//1000}kb.bedgraph"
    if not p.exists():
        return None
    d = msd.load(str(p)).dropna(subset=["Score"])
    c = msd.norm_chrom(d["Chromosome"]).to_numpy()
    d = d[c == "X"]
    gd = gene_density(d["Start"].to_numpy(), ERES)
    if gd is not None:
        ok = np.isfinite(d["Score"].to_numpy()) & np.isfinite(gd)
        if ok.sum() > 10:
            r = np.corrcoef(d["Score"].to_numpy()[ok], gd[ok])[0, 1]
            if np.isfinite(r) and r < 0:
                d = d.assign(Score=-d["Score"])
                print(f"    oriented {stem}: flipped (r to gene density {r:+.3f})")
    v = np.full(n_target, np.nan)
    for s_, val in zip(d["Start"].to_numpy(), d["Score"].to_numpy()):
        lo = int(s_) // RES
        v[lo:lo + ERES // RES] = val
    return v


def load(arm, res=RES, n_target=None):
    p = XDATA / "haplotype" / arm / f"{arm}_mergedCrush_{res}.bedgraph"
    if not p.exists():
        return None
    d = msd.load(str(p)).dropna(subset=["Score"])
    c = msd.norm_chrom(d["Chromosome"]).to_numpy()
    d = d[c == "X"]
    n = n_target if n_target else int(d["Start"].max()) // RES + 1
    v = np.full(n, np.nan)
    if res == RES:
        idx = (d["Start"].to_numpy() // RES).astype(int)
        ok = idx < n
        v[idx[ok]] = d["Score"].to_numpy()[ok]
    else:
        for s_, val in zip(d["Start"].to_numpy(), d["Score"].to_numpy()):
            lo = int(s_) // RES
            v[lo:lo + res // RES] = val
    return v


def main() -> None:
    cr = [(load(a), lab) for a, lab, _ in ARMS]
    n = max(len(v) for v, _ in cr)
    lim = CRUSH_LIM

   
    c5_all = {a: load(a, CRES, n) for a, _, _ in ARMS}
    ev_all = {stem: load_eigen(stem, n) for _, _, stem in ARMS}
    c5_lim = CRUSH_LIM
    ev_lim = max(np.nanmax(np.abs(v)) for v in ev_all.values()
                 if v is not None) * 1.05
    x = np.arange(n) * RES / 1e6
    rows = []
    for (v, lab), (arm, _, stem) in zip(cr, ARMS):
        ev = ev_all.get(stem)
        if ev is not None:
            rows.append((ev, f"{lab}\neigenvector 500 kb", ev_lim, E_POS, E_NEG,
                         f"{arm}|eigen 500 kb"))
        c5 = c5_all.get(arm)
        if c5 is not None:
            rows.append((c5, f"{lab}\nCRUSH 500 kb", c5_lim, A_COL, B_COL,
                         f"{arm}|CRUSH 500 kb"))
        rows.append((v, f"{lab}\nCRUSH 5 kb", lim, A_COL, B_COL,
                     f"{arm}|CRUSH 5 kb"))

    fig, axes = plt.subplots(len(rows), 1, figsize=(11.5, 0.95 * len(rows)),
                             sharex=True)
    fig.patch.set_facecolor("white")
    scores = json.loads(FLAGS.read_text()) if (SCORES and FLAGS.exists()) else {}
    for ax, (v, lab, ylim, cpos, cneg, skey) in zip(axes, rows):
        vv = np.full(n, np.nan); vv[:len(v)] = v
        ax.fill_between(x, 0, np.nan_to_num(vv), where=vv > 0, color=cpos, lw=0)
        ax.fill_between(x, 0, np.nan_to_num(vv), where=vv <= 0, color=cneg, lw=0)
        ax.axhline(0, color="#c9c9c4", lw=0.8, zorder=1)
        ax.axvline(DXZ4 / 1e6, color="#2E77B5", lw=1.4,
                   linestyle=(0, (4, 2.5)), zorder=4)
        ax.set_ylim(-ylim, ylim); ax.set_yticks([])
        if skey in scores:
            ax.text(1.004, 0.5, f"{scores[skey]['strength']:.2f}",
                    transform=ax.transAxes, fontsize=9, color=INK,
                    va="center", ha="left")
        if LABELS:
            ax.set_ylabel(lab, fontsize=9.5, color=INK, rotation=0,
                          ha="right", va="center", labelpad=8)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color("#d7d7d3")
        ax.tick_params(colors=INK2, labelsize=9.5, length=0)
    if LABELS:
        axes[-1].set_xlabel("chrX (Mb)", fontsize=10.5, color=INK, labelpad=6)
        axes[0].text(DXZ4 / 1e6, rows[0][2] * 0.98, " DXZ4 / megadomain boundary",
                     fontsize=9, color="#2E77B5", va="top")
    fig.subplots_adjust(left=0.13 if LABELS else 0.02,
                        right=0.955 if scores else 0.99,
                        top=0.97, bottom=0.11, hspace=0.16)
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "figS3A_x_tracks.png"
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, facecolor="white")
    if scores:
        import csv
        q = OUT / "figS3A_x_tracks_sourcedata.csv"
        with q.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["arm", "track", "AA", "BB", "AB", "strength_AABB_2AB",
                        "null_mean", "null_sd", "z", "scored_at_bp", "norm"])
            plotted = {a for a, _, _ in ARMS}
            for key, v in scores.items():
                arm, track = key.split("|", 1)
                if arm not in plotted:
                    continue
                w.writerow([arm, track, f"{v['AA']:.6f}", f"{v['BB']:.6f}",
                            f"{v['AB']:.6f}", f"{v['strength']:.6f}",
                            f"{v['null']:.6f}", f"{v['null_sd']:.6f}",
                            f"{v['z']:.3f}", v["scale"], v["norm"]])
        print(f"wrote {q}")
    print(f"wrote {out}")
    for (v, lab, *_rest) in rows:
        cov = np.isfinite(v)
        print(f"  {lab.replace(chr(10),' / '):<40}{int(cov.sum()):>8,} bins"
              f"  {100*np.mean(v[cov] > 0):>5.1f}% A")


if __name__ == "__main__":
    main()
