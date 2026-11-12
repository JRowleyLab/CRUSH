#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results"
LIB = OUT / "GO_BP_2025.gmt"
CONTRASTS = [("vs_rpos", "TNBC vs receptor-positive"),
             ("vs_hmec", "TNBC vs HMEC"),
             ("bt549_vs_hcc70", "BT549 vs HCC70")]
PANEL = "vs_rpos"
TOPN = 8
PLOT_TOPN = 4
A_COL, B_COL, INK, INK2 = "#8B0000", "#000000", "#16181c", "#5a5f68"


def build() -> None:
    """Score table -> ranks -> fgsea -> rotation, skipping what already exists."""
    if not (OUT / "figS5E_geneorder.tsv").exists():
        subprocess.run([sys.executable, str(HERE / "figS5E_genescores.py")], check=True)
    for tag, _ in CONTRASTS:
        rnk, gsea = OUT / f"figS5E_tnbc_{tag}.rnk", OUT / f"figS5E_gsea_{tag}.tsv"
        rot = OUT / f"figS5E_rotation_{tag}.tsv"
        if not gsea.exists():
            subprocess.run(["Rscript", str(HERE / "figS4_gsea.R"), str(rnk),
                            str(LIB), str(gsea)], check=True)
        if not rot.exists():
            subprocess.run([sys.executable, str(HERE / "figS4_rotation.py"),
                            str(rnk), str(gsea), str(LIB), str(rot),
                            str(OUT / "figS5E_geneorder.tsv")], check=True)


def plot(tag: str, labels: bool = True) -> None:
    d = pd.read_csv(OUT / f"figS5E_rotation_{tag}.tsv", sep="\t")
    d["short"] = d["term"].str.split(" (GO:", regex=False).str[0]
    pos = d[d["NES"] > 0].nlargest(PLOT_TOPN, "NES")
    neg = d[d["NES"] < 0].nsmallest(PLOT_TOPN, "NES")
    d = pd.concat([pos, neg]).sort_values("NES").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7.6, 0.42 * 2 * PLOT_TOPN + 1.0))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    y = np.arange(len(d))
    cols = [A_COL if v > 0 else B_COL for v in d["NES"]]
    ax.barh(y, d["NES"], 0.72, color=cols, lw=0, zorder=3)
    ax.axvline(0, color=INK, lw=0.9, zorder=4)
    if labels:
        for yi, (nes, term, p, n) in enumerate(zip(d["NES"], d["short"],
                                                   d["p_rotation"], d["n_genes"])):
            side = 0.06 if nes < 0 else -0.06
            ax.text(side, yi, term, ha="left" if nes < 0 else "right",
                    va="center", fontsize=7.5, color=INK, zorder=5)
            ax.text(nes + (-0.06 if nes < 0 else 0.06), yi,
                    f"p={p:.3f}  n={n}", ha="right" if nes < 0 else "left",
                    va="center", fontsize=6.5, color=INK2, zorder=5)
        who = "BT549" if tag == "bt549_vs_hcc70" else "TNBC"
        ax.set_xlabel(f"NES   (toward B in {who}  <-  ->  toward A in {who})",
                      fontsize=8.5, color=INK, labelpad=7)
    ax.set_yticks([]); ax.set_ylim(-0.8, len(d) - 0.2)
    lim = float(np.abs(d["NES"]).max()) * 1.55
    ax.set_xlim(-lim, lim)
    ax.tick_params(colors=INK2, labelsize=8, length=3)
    if not labels:
        ax.set_xticklabels([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(INK)
    fig.tight_layout()
    suffix = "" if labels else "_nolabels"
    panel = "figS5F" if tag == "bt549_vs_hcc70" else "figS5E"
    png = OUT / f"{panel}_gsea_{tag}{suffix}.png"
    fig.savefig(png, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {png.name}")


def main() -> None:
    build()
    for tag, name in CONTRASTS:
        d = pd.read_csv(OUT / f"figS5E_rotation_{tag}.tsv", sep="\t")
        sig = int((d["p_rotation"] <= 0.05).sum())
        print(f"\n  {name}: {sig}/{len(d)} of the top terms at p(rot) <= 0.05")
        for _, r in d.sort_values("NES", ascending=False).iterrows():
            mark = "*" if r["p_rotation"] <= 0.05 else " "
            print(f"   {mark} {r['term'].split(' (GO:')[0][:52]:<52}"
                  f"{r['NES']:>7.2f}{int(r['n_genes']):>5}{r['p_rotation']:>8.3f}")
    for tag in (PANEL, "bt549_vs_hcc70"):
        plot(tag, labels=True)
        plot(tag, labels=False)


if __name__ == "__main__":
    main()
