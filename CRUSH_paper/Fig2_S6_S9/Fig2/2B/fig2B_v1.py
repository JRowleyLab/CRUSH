#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fig2B

ORIG = Path("normalized_FM_Subs_E_CR_Cs")
HERE = fig2B.HERE
fig2B.WORK = HERE / "work_v1"
CLASSES, COLORS = fig2B.CLASSES, fig2B.COLORS


def load_calls():
    rd = lambda f, v: pd.read_csv(ORIG / f, sep="\t", header=None, usecols=[0, 1, 3],
                                  names=["c", "s", v], dtype={"c": str})
    e, cr = rd("Eigen_1kb_fm_N.bedgraph", "eig"), rd("Crush_1kb_fm_N.bedgraph", "cr")
    for t in (e, cr):
        t["c"] = t["c"].str.replace("chr", "", regex=False)
    cr = cr[cr["c"].isin(fig2B.CHROMS) & (cr["cr"] != 0)]
    d = e[e["c"].isin(fig2B.CHROMS)].merge(cr, on=["c", "s"])
    d = d[d["eig"] != 0].reset_index(drop=True)
    d["cls"] = (np.where(d["eig"] > 0, "A", "B").astype(object) + "->" +
                np.where(d["cr"] > 0, "A", "B").astype(object))
    return d, cr


def zscore(x):
    return (x - np.nanmean(x)) / np.nanstd(x)


def boxes(ax, d, col):
    data = [d.loc[d["cls"] == k, col].dropna().to_numpy() for k in CLASSES]
    bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black"))
    for b, k in zip(bp["boxes"], CLASSES):
        b.set_facecolor(COLORS[k])


def style(ax):
    ax.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    d, cr = load_calls()
    g = pd.read_csv(ORIG / "histones_signal" / "fmCR1kb_gmgro", sep="\t", usecols=["chr", "start", "means1"],
                    dtype={"chr": str}).rename(columns={"chr": "c", "start": "s"})
    d = d.merge(g, on=["c", "s"], how="left")
    d["gro"] = zscore(d["means1"].to_numpy())
    c1 = d["c"] == "1"
    d["gro_chr1"] = np.nan
    d.loc[c1, "gro_chr1"] = zscore(d.loc[c1, "gro"].to_numpy())
    ab = fig2B.contact_sums(cr)
    d = d.merge(ab, on=["c", "s"], how="left")

    res = {}
    for k in CLASSES:
        s = d[d["cls"] == k]
        q = s["gro"].quantile([0.25, 0.5, 0.75])
        q1 = s["gro_chr1"].quantile([0.25, 0.5, 0.75])
        res[k] = {"n_bins": int(len(s)), "n_bins_chr1": int((s["c"] == "1").sum()),
                  "gro_z_q25_median_q75": [float(x) for x in q],
                  "gro_z_chr1_q25_median_q75": [float(x) for x in q1],
                  "contacts_to_A": float(s["toA"].sum()), "contacts_to_B": float(s["toB"].sum()),
                  "contact_ratio_log2": float(np.log2(s["toA"].sum() / s["toB"].sum()))}
    (HERE / "fig2B_v1.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))

    r = [res[k]["contact_ratio_log2"] for k in CLASSES]
    for notext in (False, True):
        fig, ax = plt.subplots(2, 1, figsize=(2.6, 4.4), gridspec_kw={"hspace": 0.35})
        boxes(ax[0], d, "gro")
        ax[0].set_xticks([])
        ax[1].bar(range(1, 5), r, color=[COLORS[k] for k in CLASSES], width=0.6)
        ax[1].axhline(0, color="black", lw=0.6)
        ax[1].set_xticks(range(1, 5))
        for a in ax:
            style(a)
        if notext:
            ax[1].set_xticklabels([])
        else:
            ax[0].set_ylabel("GRO-seq\n(z-scores)", fontsize=8)
            ax[1].set_ylabel("Contact ratio\n(log2 A/B)", fontsize=8)
            ax[1].set_xticklabels(CLASSES, fontsize=7)
            ax[1].set_xlabel("Eigen. -> CRUSH", fontsize=8)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"fig2B_v1{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(2.6, 2.2))
    boxes(ax, d[c1], "gro_chr1")
    ax.set_xticks(range(1, 5))
    ax.set_xticklabels(CLASSES, fontsize=7)
    ax.set_ylabel("GRO-seq (z-scores)\nchr1 only", fontsize=8)
    style(ax)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"fig2B_v1_chr1.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
