#!/usr/bin/env python3
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
REV = FIG.parent
sys.path.insert(0, str(REV))
import MSD_tracks as msd

RES = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
RMIN = 0.05
THRESH = 50.0
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
K27 = Path("tracks")
GENES = "hg19genes.bed"
AB, BA, ST, INK, INK2 = "#45BFDE", "#E9483D", "#b8b8b2", "#16181c", "#5a5f68"
CLASSES = [("A to B", "AB", AB), ("Stable", "ST", ST), ("B to A", "BA", BA)]


def gene_density():
    g = pd.read_csv(GENES, sep="\t", header=None, usecols=[0, 1],
                    names=["c", "s"], dtype={"c": str})
    g["c"] = g["c"].str.removeprefix("chr")
    g["b"] = (g["s"] // RES).astype(np.int64)
    return g.groupby(["c", "b"]).size().rename("n").reset_index()


def crush(line, dens):
    p = REV / "tracks" / "tnbc" / f"{line}_r{RES}" / f"{line}_mergedCrush_{RES}.bedgraph"
    d = msd.load(str(p)).dropna(subset=["Score"])
    t = pd.DataFrame({"c": msd.norm_chrom(d["Chromosome"]).to_numpy(),
                      "b": (d["Start"].to_numpy() // RES).astype(np.int64),
                      "v": d["Score"].to_numpy(float)})
    t = t[(t["v"] != 0) & t["c"].isin(CHROMS)].reset_index(drop=True)
    m = t.merge(dens, on=["c", "b"], how="left").fillna({"n": 0.0})
    flipped = []
    for c, sub in m.groupby("c"):
        if len(sub) < 50 or sub["n"].std() == 0:
            continue
        r = np.corrcoef(sub["v"], sub["n"])[0, 1]
        if np.isfinite(r) and abs(r) >= RMIN and r < 0:
            t.loc[t["c"] == c, "v"] *= -1
            flipped.append(c)
    print(f"    {line:<7} {len(t):>9,} bins, {len(flipped)} chroms flipped", flush=True)
    return t


def k27(path):
    bw = pyBigWig.open(str(path))
    sz = bw.chroms()
    rows = []
    for c in CHROMS:
        k = c if c in sz else "chr" + c
        if k not in sz:
            continue
        n = sz[k] // RES + 1
        v = bw.stats(k, 0, min(n * RES, sz[k]), type="mean", nBins=n)
        v = np.array([np.nan if x is None else x for x in v], float)
        rows.append(pd.DataFrame({"c": c, "b": np.arange(len(v)), "k": v}))
    bw.close()
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    print(f"\n  {RES // 1000} kb")
    dens = gene_density()
    h, b = crush("hmec", dens), crush("bt549", dens)
    m = h.merge(b, on=["c", "b"], suffixes=("_h", "_b"))
    kh = k27(K27 / "HMEC_K27ac.CPM.bw").rename(columns={"k": "kh"})
    kb = k27(K27 / "BT549_K27ac.CPM.bw").rename(columns={"k": "kb"})
    m = m.merge(kh, on=["c", "b"]).merge(kb, on=["c", "b"]).dropna()
    m["delta"] = m["kb"] - m["kh"]
    h_a, h_b = m["v_h"] >= THRESH, m["v_h"] <= -THRESH
    b_a, b_b = m["v_b"] >= THRESH, m["v_b"] <= -THRESH
    m["cls"] = np.where(h_a & b_b, "AB",
                np.where(h_b & b_a, "BA",
                np.where((h_a & b_a) | (h_b & b_b), "ST", "AMB")))
    n_amb = int((m["cls"] == "AMB").sum())
    print(f"\n  {n_amb:,} of {len(m):,} bins ambiguous (|score| < {THRESH:g} in "
          f"one or both lines) -- dropped")
    m = m[m["cls"] != "AMB"]

    groups = {k: m.loc[m["cls"] == k, "delta"].to_numpy() for _, k, _ in CLASSES}
    rec = {}
    print(f"\n  {'class':<8}{'n':>9}{'median':>10}{'mean':>10}")
    for name, key, _ in CLASSES:
        g = groups[key]
        rec[key] = {"n": int(len(g)), "median": float(np.median(g)),
                    "mean": float(g.mean())}
        print(f"  {name:<8}{len(g):>9,}{np.median(g):>10.3f}{g.mean():>10.3f}")
    for a, bkey in (("AB", "ST"), ("BA", "ST"), ("AB", "BA")):
        u, p = mannwhitneyu(groups[a], groups[bkey])
        rec[f"{a}_vs_{bkey}"] = {"p": float(p)}
        print(f"  {a} vs {bkey}: Mann-Whitney p = {p:.3g}")

    fig, ax = plt.subplots(figsize=(3.6, 4.0))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    data = [groups[k] for _, k, _ in CLASSES]
    bp = ax.boxplot(data, widths=0.6, showfliers=False, patch_artist=True,
                    medianprops=dict(color=INK, lw=1.4),
                    whiskerprops=dict(color=INK2, lw=1.0),
                    capprops=dict(color=INK2, lw=1.0),
                    boxprops=dict(lw=1.0, edgecolor=INK2))
    for patch, (_, _, col) in zip(bp["boxes"], CLASSES):
        patch.set_facecolor(col)
    ax.axhline(0, color=INK2, lw=0.8, ls="--", zorder=1)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([n for n, _, _ in CLASSES], fontsize=9, color=INK)
    ax.set_ylabel("Δ H3K27ac (BT549 − HMEC, CPM)", fontsize=9, color=INK,
                  labelpad=6)
    ax.tick_params(colors=INK2, labelsize=8.5, length=3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d7d7d3")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    png = OUT / f"figS5D_switch_k27ac_{RES // 1000}kb.png"
    fig.savefig(png, dpi=300, facecolor="white")

    ax.set_ylabel(""); ax.set_xticklabels([""] * 3); ax.set_yticklabels([])
    fig.savefig(OUT / f"figS5D_switch_k27ac_{RES // 1000}kb_nolabels.png",
                dpi=300, facecolor="white")
    plt.close(fig)
    (OUT / f"figS5D_switch_k27ac_{RES // 1000}kb.json").write_text(
        json.dumps(rec, indent=1))
    print(f"\nwrote {png.name} (+ _nolabels, .json)")


if __name__ == "__main__":
    main()
