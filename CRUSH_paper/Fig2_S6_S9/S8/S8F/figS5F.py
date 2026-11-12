#!/usr/bin/env python3
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import fig2D

HERE = Path(__file__).resolve().parent
RECALL = HERE.parent / "recall"
TITLE = {"v1rust": "CRUSH --v1 (Rust port of the published pipeline)",
         "fixed": "CRUSH v2, --seed-select both --orient 1"}


def stars(p: float) -> str:
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def main():
    cols = ["chr", "start", "end", "hct_c", "hela_c", "hepg2_c", "imr90_c", "k562_c",
            "chr2", "start2", "end2", "hct", "hela", "hepg2", "imr90", "k562"]
    df = pd.read_csv(fig2D.MATRIX, sep="\t", header=None, names=cols,
                     usecols=["chr", "start"] + fig2D.LINES, dtype={"chr": str})
    bru = df[fig2D.LINES].apply(pd.to_numeric, errors="coerce").to_numpy(float)

    for mode in sys.argv[1:] or ("highpass", "coarse"):
        fig2D.TRACKS = RECALL / mode
        comp = np.column_stack([fig2D.crush(l, df["chr"], df["start"]) for l in fig2D.LINES])
        comp[comp == 0] = np.nan
        res = {}
        fig, ax = plt.subplots(2, 4, figsize=(8.0, 4.6), sharey="row")
        for j, l in enumerate(fig2D.LINES):
            thr = np.nanpercentile(bru[:, j], fig2D.BRU_PCT)
            others = np.delete(np.arange(len(fig2D.LINES)), j)
            uniq = (bru[:, j] >= thr) & np.all(bru[:, others] <= 0, axis=1)
            ok = np.all(np.isfinite(comp), axis=1) & np.isfinite(bru[:, j])
            dbru = bru[:, j] - np.nanmean(bru[:, others], axis=1)
            dcr = comp[:, j] - np.nanmean(comp[:, others], axis=1)
            res[l] = {}
            for r, (name, delta) in enumerate((("ΔBRU-seq", dbru), ("ΔCRUSH 1 kb", dcr))):
                keep = ok & np.isfinite(delta)
                o, u = delta[~uniq & keep], delta[uniq & keep]
                p = mannwhitneyu(u, o, alternative="greater").pvalue
                res[l][name] = {"median_other": float(np.median(o)),
                                "median_unique": float(np.median(u)),
                                "n_unique": int(u.size), "p_one_sided": float(p)}
                bp = ax[r, j].boxplot([o, u], showfliers=False, widths=0.6, patch_artist=True,
                                      medianprops=dict(color="black"))
                bp["boxes"][0].set_facecolor("#cccccc")
                bp["boxes"][1].set_facecolor(fig2D.COLOR[l])
                ax[r, j].axhline(0, ls="--", lw=0.6, color="black")
                ax[r, j].set_title(f"{stars(p)}", fontsize=8, pad=2)
                ax[r, j].set_xticks([1, 2])
                ax[r, j].set_xticklabels(["Other", "Unique"] if r else ["", ""], fontsize=7)
                ax[r, j].tick_params(labelsize=7)
                for s in ("top", "right"):
                    ax[r, j].spines[s].set_visible(False)
                if j == 0:
                    ax[r, j].set_ylabel(f"{name}\n(line - mean of others)", fontsize=8)
            ax[0, j].text(0.5, 1.18, f"{fig2D.LABEL[l]} (n={res[l]['ΔCRUSH 1 kb']['n_unique']:,})",
                          transform=ax[0, j].transAxes, ha="center", fontsize=8)
        fig.suptitle(TITLE.get(mode, f"CRUSH v0.2, -M {mode}"), fontsize=8, y=1.02)
        fig.tight_layout()
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"figS5F_{mode}.{ext}", dpi=300, bbox_inches="tight")
        (HERE / f"figS5F_{mode}.json").write_text(json.dumps(res, indent=1))
        print(mode, {l: round(res[l]["ΔCRUSH 1 kb"]["median_unique"], 1) for l in res})


if __name__ == "__main__":
    main()
