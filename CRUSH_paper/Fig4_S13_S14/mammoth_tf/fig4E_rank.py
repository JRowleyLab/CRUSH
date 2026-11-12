#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
CT = Path("comptracks")
OLD = ["NR1D2", "BHLHE40", "ADNP", "THAP5", "POU6F2", "ZNF215", "FOXR2", "GCM2"]
NTOP = 10
AMIN, VMAX = 0.25, 3
E_COL, M_COL = "#b5307a", "#2f6fe0"


def spct(x):
    return 2 * (x.rank() - 0.5) / len(x) - 1


def main():
    r = pd.read_csv(HERE / "regulon_rank.tsv", sep="\t")
    r = r[(r.regulon == "all") & r.TF_dCRUSH_norm.notna()].copy()
    r["TF_score"] = spct(r.TF_dCRUSH_norm)
    r["target_score"] = spct(r.rank_biserial)
    r["combined"] = (r.TF_score + r.target_score) / 2
    r["species"] = np.where(r.combined > 0, "mammoth", "elephant")
    r = r.sort_values("combined").reset_index(drop=True)
    r["rank"] = np.arange(1, len(r) + 1)
    g = pd.read_csv(CT / "genes_median", sep="\t", header=None, names=["gene", "e", "m"])
    g = g.drop_duplicates("gene").set_index("gene")
    r["elephant_CRUSH"], r["mammoth_CRUSH"] = r.TF.map(g.e), r.TF.map(g.m)
    r["published_fig4E"] = np.where(r.TF.isin(OLD), "yes", "")
    top = pd.concat([r.head(NTOP), r.tail(NTOP)])
    r.to_csv(HERE / "fig4E_rank.tsv", sep="\t", index=False, float_format="%.4g")
    print(top[["rank", "TF", "TF_dCRUSH_norm", "n_targets", "n_more_A_elephant", "n_more_A_mammoth",
               "rank_biserial", "FDR_rotation", "combined", "TableS9_set", "published_fig4E"]].to_string(index=False))
    print("\npublished Fig. 4E TFs:")
    print(r[r.TF.isin(OLD)][["rank", "TF", "TF_dCRUSH_norm", "rank_biserial", "FDR_rotation", "combined",
                             "concordant"]].to_string(index=False))

    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(3.2, 2.6))
        col = np.where(r.combined > 0, M_COL, E_COL)
        sig = (r.FDR_rotation < 0.05).to_numpy()
        ax.scatter(r["rank"][~sig], r.combined[~sig], s=4, facecolor="none", edgecolor=col[~sig], lw=0.3)
        ax.scatter(r["rank"][sig], r.combined[sig], s=4, color=col[sig], lw=0)
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xlim(0, len(r) + 1)
        ax.set_ylim(-1.05, 1.05)
        ax.tick_params(labelsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if notext:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            ax.set_xlabel("TF rank", fontsize=8)
            ax.set_ylabel("Combined score (TF + targets)", fontsize=8)
            for side, x in (("elephant", r.head(NTOP)), ("mammoth", r.tail(NTOP))):
                y0 = -0.95 if side == "elephant" else 0.95
                xt = len(r) * (0.62 if side == "elephant" else 0.04)
                for i, (_, t) in enumerate(x.iterrows()):
                    ax.text(xt, y0 + (0.09 * i if side == "elephant" else -0.09 * i), t.TF, fontsize=5,
                            style="italic", color=E_COL if side == "elephant" else M_COL,
                            va="center")
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"fig4E_rank{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    top = top.assign(mlog10_FDR=-np.log10(top.FDR_rotation))
    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(5.0, 2.2))
        x = np.arange(len(top))
        for k, (c, base) in enumerate((("elephant_CRUSH", E_COL), ("mammoth_CRUSH", M_COL))):
            for xi, (_, t) in zip(x, top.iterrows()):
                a = AMIN + (1 - AMIN) * min(t.mlog10_FDR / VMAX, 1)
                ax.bar(xi + (k - 0.5) * 0.38, t[c], width=0.38, color=base, alpha=a, lw=0)
        ax.axhline(0, color="black", lw=0.6)
        ax.axvline(NTOP - 0.5, color="grey", lw=0.5, ls=":")
        ax.set_xticks(x)
        ax.set_xlim(-0.6, len(top) - 0.4)
        ax.tick_params(labelsize=7)
        for sp in ("top", "right", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="x", length=0)
        cax = fig.add_axes([0.93, 0.25, 0.018, 0.45])
        grey = LinearSegmentedColormap.from_list("g", [(0, 0, 0, AMIN), (0, 0, 0, 1)])
        cb = fig.colorbar(plt.cm.ScalarMappable(cmap=grey, norm=plt.Normalize(0, VMAX)), cax=cax)
        cb.outline.set_linewidth(0.4)
        cb.set_ticks([0, -np.log10(0.05), VMAX])
        cb.ax.tick_params(labelsize=6, length=2)
        if notext:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            cb.set_ticklabels([])
        else:
            cb.ax.set_yticklabels(["0", "1.3", f"{VMAX:g}"])
            cb.set_label("Targets, -log$_{10}$(FDR)", fontsize=7)
            ax.set_xticklabels(top.TF, rotation=60, ha="right", fontsize=6, style="italic")
            ax.set_ylabel("CRUSH @ 5 kb", fontsize=8)
            ax.text(0.25, 1.02, "More A in elephant (TF + targets)", transform=ax.transAxes, ha="center",
                    fontsize=7, color=E_COL)
            ax.text(0.75, 1.02, "More A in mammoth (TF + targets)", transform=ax.transAxes, ha="center",
                    fontsize=7, color=M_COL)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"fig4E_top{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
