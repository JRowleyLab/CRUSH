#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
MATRIX = Path("comp_bru_matrix.tsv")
TRACKS = Path("encode")
LINES = ["hct", "hepg2", "imr90", "k562"]
LABEL = {"hct": "HCT-116", "hepg2": "HEPG2", "imr90": "IMR90", "k562": "K562"}
COLOR = {"hct": "tab:red", "hepg2": "tab:green", "imr90": "tab:purple", "k562": "tab:orange"}
RES = 1000
BRU_PCT = 95.0


def crush(line: str, chrom: pd.Series, start: pd.Series) -> np.ndarray:
    bw = pyBigWig.open(str(TRACKS / f"{line}_r1000" / f"{line}_mergedCrush_{RES}.bw"))
    out = np.full(len(chrom), np.nan)
    for c, idx in chrom.groupby(chrom).groups.items():
        if c not in bw.chroms():
            continue
        v = np.full(bw.chroms()[c] // RES + 1, np.nan)
        iv = np.array(bw.intervals(c) or [], dtype=float).reshape(-1, 3)
        v[(iv[:, 0] // RES).astype(int)] = iv[:, 2]
        b = start.loc[idx].to_numpy() // RES
        ok = b < len(v)
        vals = np.full(len(b), np.nan)
        vals[ok] = v[b[ok]]
        out[np.asarray(idx)] = vals
    bw.close()
    return out


def main() -> None:
    cols = ["chr", "start", "end", "hct_c", "hela_c", "hepg2_c", "imr90_c", "k562_c",
            "chr2", "start2", "end2", "hct", "hela", "hepg2", "imr90", "k562"]
    df = pd.read_csv(MATRIX, sep="\t", header=None, names=cols,
                     usecols=["chr", "start"] + LINES, dtype={"chr": str})
    bru = df[LINES].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    comp = np.column_stack([crush(l, df["chr"], df["start"]) for l in LINES])
    comp[comp == 0] = np.nan

    res, uniq_all, other_all = {}, [], []
    for j, l in enumerate(LINES):
        thr = np.nanpercentile(bru[:, j], BRU_PCT)
        others = np.delete(np.arange(len(LINES)), j)
        uniq = (bru[:, j] >= thr) & np.all(bru[:, others] <= 0, axis=1)
        delta = comp[:, j] - np.nanmean(comp[:, others], axis=1)
        ok = np.isfinite(delta) & np.all(np.isfinite(comp), axis=1)
        du, do = delta[uniq & ok], delta[~uniq & ok]
        uniq_all.append(du)
        other_all.append(do)
        res[l] = {"n_unique": int(du.size), "n_other": int(do.size),
                  "median_unique": float(np.median(du)), "median_other": float(np.median(do))}
    u, o = np.concatenate(uniq_all), np.concatenate(other_all)
    p = mannwhitneyu(u, o, alternative="greater").pvalue
    res["pooled"] = {"median_unique": float(np.median(u)), "median_other": float(np.median(o)),
                     "p_one_sided": float(p)}
    (HERE / "fig2D.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))

    fig, ax = plt.subplots(figsize=(3.4, 3.4))
    ax.boxplot([o, u], showfliers=False, widths=0.55,
               medianprops=dict(color="black"), boxprops=dict(color="black"),
               whiskerprops=dict(color="black"), capprops=dict(color="black"))
    ax.axhline(0, ls="--", lw=0.8, color="black")
    for off, l in zip(np.linspace(-0.12, 0.12, len(LINES)), LINES):
        ax.scatter([1 + off, 2 + off], [res[l]["median_other"], res[l]["median_unique"]],
                   s=36, color=COLOR[l], edgecolor="black", linewidth=0.4,
                   zorder=3, label=LABEL[l])
    top = ax.get_ylim()[1]
    ax.plot([1, 1, 2, 2], [top * 0.92, top * 0.95, top * 0.95, top * 0.92], color="black", lw=0.8)
    ax.text(1.5, top * 0.96, "***", ha="center", va="bottom")
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Commonly\nTranscribed", "Uniquely\nTranscribed"], fontsize=8)
    ax.set_ylabel("ΔCRUSH 1 kb\n(Unique Line - Average)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"fig2D.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
