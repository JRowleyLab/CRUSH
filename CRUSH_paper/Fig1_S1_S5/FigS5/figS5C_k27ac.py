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
from sklearn.metrics import roc_auc_score, roc_curve

REV = Path("revision")
sys.path.insert(0, str(REV))
import MSD_tracks as msd

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
RES = 100_000
RMIN = 0.05
QCUT = 0.75
CR_COL, EV_COL, INK, INK2 = "#2E5C8A", "#9a9a94", "#16181c", "#5a5f68"
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
CHROM = "10"
UNIFORM_AS_A = True

LINES = [("hmec",    "HMEC",    "GSM5098062_HMEC",    "GSM5098095"),
         ("t47D",    "T47D",    "GSM5098064_T47D",    "GSM5098097"),
         ("zr7530",  "ZR7530",  "GSM5098066_ZR7530",  "GSM5098099"),
         ("hcc1954", "HCC1954", "GSM5098068_HCC1954", "GSM5098101"),
         ("hcc70",   "HCC70",   "GSM5098070_HCC70",   "GSM5098103"),
         ("bt549",   "BT549",   "GSM5098072_BT549",   "GSM5098105")]
EIGDIR = Path("eigen_vcsqrt")
K27DIR = Path("ChIP")
K27NAME = {"HMEC": "HMEC", "T47D": "T47D", "ZR7530": "ZR7530",
           "HCC1954": "HCC1954", "HCC70": "HCC70", "BT549": "BT549"}
GENES = ("hg19genes.bed")


def gene_density():
    g = pd.read_csv(GENES, sep="\t", header=None, usecols=[0, 1],
                    names=["c", "s"], dtype={"c": str})
    g["c"] = g["c"].str.removeprefix("chr")
    g["b"] = (g["s"] // RES).astype(np.int64)
    return g.groupby(["c", "b"]).size().rename("n").reset_index()


def keyed(path):
    d = msd.load(str(path)).dropna(subset=["Score"])
    t = pd.DataFrame({"c": msd.norm_chrom(d["Chromosome"]).to_numpy(),
                      "b": (d["Start"].to_numpy() // RES).astype(np.int64),
                      "v": d["Score"].to_numpy(float)})
    return t[(t["v"] != 0) & t["c"].isin(CHROMS)].reset_index(drop=True)


def orient(t, dens):
    """Flip each chromosome so positive score tracks gene density. Chromosomes
    where |r| < RMIN are left alone -- the rule cannot decide there."""
    t = t.copy()
    m = t.merge(dens, on=["c", "b"], how="left").fillna({"n": 0.0})
    flipped, undet = [], []
    for c, sub in m.groupby("c"):
        if len(sub) < 50 or sub["n"].std() == 0:
            continue
        r = np.corrcoef(sub["v"], sub["n"])[0, 1]
        if not np.isfinite(r) or abs(r) < RMIN:
            undet.append(c)
            if UNIFORM_AS_A and (sub["v"] > 0).mean() < 0.5:
                t.loc[t["c"] == c, "v"] *= -1
            continue
        if r < 0:
            t.loc[t["c"] == c, "v"] *= -1
            flipped.append(c)
    return t, flipped, undet


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


def score(t, kk):
    d = t.merge(kk, on=["c", "b"], how="inner").dropna(subset=["v", "k"])
    y = (d["k"] >= np.quantile(d["k"], QCUT)).astype(int).to_numpy()
    return y, d["v"].to_numpy(), len(d)


BOXED = {"T47D", "BT549"}


def plot_bars(out, png, labels=True):
    import matplotlib.transforms as mtrans
    labs = [l for _, l, _, _ in LINES]
    fig, ax = plt.subplots(figsize=(5.4, 2.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    x = np.arange(len(labs))
    w = 0.40
    for j, (key, col) in enumerate((("CRUSH", CR_COL), ("Eigenvector", EV_COL))):
        ys = [out[l][key]["auroc"] for l in labs]
        ax.bar(x + (j - 0.5) * w, ys, w * 0.92, color=col, lw=0, zorder=3)
        if labels:
            for xi, yi in zip(x + (j - 0.5) * w, ys):
                ax.text(xi, yi + 0.02, f"{yi:.2f}", ha="center", va="bottom",
                        fontsize=6.5, color=INK2, zorder=4)
    tr = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
    for i, l in enumerate(labs):
        if l not in BOXED:
            continue
        ax.add_patch(plt.Rectangle((i - 0.50, -0.30), 1.00, 1.32, transform=tr,
                                   fill=False, ec=INK, lw=1.0, ls=(0, (1.6, 1.6)),
                                   clip_on=False, zorder=5))
    ax.set_xticks(x)
    ax.set_xticklabels(labs if labels else [""] * len(labs), fontsize=8,
                       color=INK, rotation=45, ha="right")
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "0.5", "1"] if labels else ["", "", ""])
    ax.set_ylim(0, 1.06)
    ax.set_xlim(-0.7, len(labs) - 0.3)
    if labels:
        ax.set_ylabel("Predictive of\nH3K27ac (AUROC)", fontsize=8.5, color=INK,
                      labelpad=6)
    ax.tick_params(colors=INK2, labelsize=8, length=3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK)
    if labels:
        for k, (txt, col) in enumerate(((f"chr{CHROM}", INK), ("CRUSH", CR_COL),
                                        ("Eigenvector", EV_COL))):
            ax.text(1.03, 0.98 - 0.13 * k, txt, transform=ax.transAxes,
                    fontsize=8.5, color=col, ha="left", va="top")
    fig.tight_layout()
    fig.savefig(png, dpi=300, facecolor="white", bbox_inches="tight")
    print(f"wrote {png}")


def main() -> None:
    dens = gene_density()
    out = {}
    roc = {}
    print(f"\n  positives = top {100*(1-QCUT):.0f}% of 100 kb bins by H3K27ac")
    print(f"  chr{CHROM} is the scored scope; genome-wide shown as a check\n")
    print(f"  {'line':<9}{'n bins':>8}{'CRUSH':>8}{'eigen':>8}"
          f"{'CRUSH gw':>10}{'eigen gw':>10}   chr{CHROM} eigen")
    for line, lab, elab, gsm in LINES:
        kk = k27(K27DIR / f"{gsm}_ChIP-seq_H3K27ac_{K27NAME[lab]}.bw")
        cr_raw = keyed(REV / "tracks" / "tnbc" / f"{line}_r100000"
                       / f"{line}_mergedCrush_{RES}.bedgraph")
        ev_raw = keyed(EIGDIR / elab / f"Eigen_{elab}_100kb.bedgraph")
        cr, crf, cru = orient(cr_raw[cr_raw["c"] == CHROM], dens)
        ev, evf, evu = orient(ev_raw[ev_raw["c"] == CHROM], dens)
        cr_gw, _, _ = orient(cr_raw, dens)
        ev_gw, _, _ = orient(ev_raw, dens)
        rec, curves = {}, {}
        for name, t in (("CRUSH", cr), ("Eigenvector", ev)):
            y, v, n = score(t, kk)
            rec[name] = {"auroc": float(roc_auc_score(y, v)), "n": int(n)}
            fpr, tpr, _ = roc_curve(y, v)
            curves[name] = (fpr, tpr)
        for name, t in (("CRUSH", cr_gw), ("Eigenvector", ev_gw)):
            y, v, _ = score(t, kk)
            rec[name]["auroc_genomewide"] = float(roc_auc_score(y, v))
        for name, t in (("CRUSH", cr_raw[cr_raw["c"] == CHROM]),
                        ("Eigenvector", ev_raw[ev_raw["c"] == CHROM])):
            y, v, _ = score(t, kk)
            rec[name]["auroc_unoriented"] = float(roc_auc_score(y, v))
        rec["oriented"] = {"crush_flipped": len(crf), "crush_undetermined": len(cru),
                           "eigen_flipped": len(evf), "eigen_undetermined": len(evu)}
        out[lab], roc[lab] = rec, curves
        print(f"  {lab:<9}{rec['CRUSH']['n']:>8,}"
              f"{rec['CRUSH']['auroc']:>8.3f}{rec['Eigenvector']['auroc']:>8.3f}"
              f"{rec['CRUSH']['auroc_genomewide']:>10.3f}"
              f"{rec['Eigenvector']['auroc_genomewide']:>10.3f}"
              f"   {'UNDETERMINED, set A' if evu else 'oriented'}")

    res = OUT
    res.mkdir(exist_ok=True)
    plot_bars(out, res / "figS5C_k27ac_auroc.png")
    plot_bars(out, res / "figS5C_k27ac_auroc_nolabels.png", labels=False)

    fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.6))
    fig.patch.set_facecolor("white")
    for ax, (_, lab, _, _) in zip(axes.ravel(), LINES):
        ax.set_facecolor("white")
        for name, col in (("CRUSH", CR_COL), ("Eigenvector", EV_COL)):
            fpr, tpr = roc[lab][name]
            ax.plot(fpr, tpr, color=col, lw=1.9, zorder=3,
                    label=f"{name} (AUROC = {out[lab][name]['auroc']:.3f})")
        ax.plot([0, 1], [0, 1], color="#d7d7d3", lw=0.9, ls="--", zorder=2)
        ax.set_title(lab, fontsize=11, color=INK, pad=6)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("False Positive Rate", fontsize=9.5, color=INK)
        ax.set_ylabel("True Positive Rate", fontsize=9.5, color=INK)
        ax.legend(frameon=False, fontsize=8.5, loc="lower right", handlelength=1.2)
        ax.tick_params(colors=INK2, labelsize=9, length=3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color("#d7d7d3")
    fig.tight_layout()
    fig.savefig(res / "figS5C_k27ac_roc.png", dpi=300, facecolor="white")
    (res / "figS5C_k27ac.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {res / 'figS5C_k27ac_roc.png'}")
    print(f"wrote {res / 'figS5C_k27ac.json'}")


if __name__ == "__main__":
    main()
