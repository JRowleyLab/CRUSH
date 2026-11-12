#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
TFS = ["BHLHE40", "ADNP", "GMEB2", "POU6F2", "FOXR2", "OTX2"]
E_COL, M_COL = "#b5307a", "#2f6fe0"
AMIN = 0.25


def main():
    r = pd.read_csv(HERE / "fig4E_rank.tsv", sep="\t").set_index("TF")
    n = len(r)
    t = r.loc[TFS].copy()
    t["rank_from_end"] = np.where(t.combined < 0, t["rank"], n + 1 - t["rank"])
    elep = t[t.combined < 0].sort_values("rank_from_end")
    mamm = t[t.combined > 0].sort_values("rank_from_end", ascending=False)
    t = pd.concat([elep, mamm])
    t[["rank_from_end", "elephant_CRUSH", "mammoth_CRUSH", "n_targets", "n_more_A_elephant", "n_more_A_mammoth",
       "rank_biserial", "FDR_rotation", "combined"]].to_csv(HERE / "fig4E_new.tsv", sep="\t", float_format="%.4g")
    print(t[["rank_from_end", "elephant_CRUSH", "mammoth_CRUSH", "rank_biserial", "FDR_rotation"]])
    t["mlog10_FDR"] = -np.log10(t.FDR_rotation)
    vmax = np.ceil(t.mlog10_FDR.max() * 2) / 2

    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(3.4, 2.4))
        x = np.arange(len(t))
        for xi, (_, s) in zip(x, t.iterrows()):
            a = AMIN + (1 - AMIN) * min(s.mlog10_FDR / vmax, 1)
            ax.bar(xi - 0.19, s.elephant_CRUSH, width=0.38, color=E_COL, alpha=a, lw=0)
            ax.bar(xi + 0.19, s.mammoth_CRUSH, width=0.38, color=M_COL, alpha=a, lw=0)
        ax.axhline(0, color="black", lw=0.8)
        ax.axvline(len(elep) - 0.5, color="grey", lw=0.5, ls=":")
        ax.set_xticks(x)
        ax.set_xlim(-0.6, len(t) - 0.4)
        ax.tick_params(labelsize=7)
        for sp in ("top", "right", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="x", length=0)
        cax = fig.add_axes([0.93, 0.25, 0.025, 0.45])
        grey = LinearSegmentedColormap.from_list("g", [(0, 0, 0, AMIN), (0, 0, 0, 1)])
        cb = fig.colorbar(plt.cm.ScalarMappable(cmap=grey, norm=plt.Normalize(0, vmax)), cax=cax)
        cb.outline.set_linewidth(0.4)
        cb.set_ticks([0, -np.log10(0.05), vmax])
        cb.ax.tick_params(labelsize=6, length=2)
        if notext:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            cb.set_ticklabels([])
        else:
            ax.set_xticklabels(t.index, rotation=60, ha="right", fontsize=7, style="italic")
            ax.set_ylabel("CRUSH @ 5 kb", fontsize=8)
            cb.ax.set_yticklabels(["0", "1.3", f"{vmax:g}"])
            cb.set_label("Targets, -log$_{10}$(FDR)", fontsize=7)
            ax.text(0.25, 1.02, "A in elephant", transform=ax.transAxes, ha="center", fontsize=7, color=E_COL)
            ax.text(0.75, 1.02, "A in mammoth", transform=ax.transAxes, ha="center", fontsize=7, color=M_COL)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"fig4E_new{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
