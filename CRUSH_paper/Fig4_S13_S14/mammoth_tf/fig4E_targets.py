#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fig4E_new as f

HERE = Path(__file__).resolve().parent


def main():
    r = pd.read_csv(HERE / "fig4E_rank.tsv", sep="\t").set_index("TF")
    n = len(r)
    t = r.loc[f.TFS].copy()
    t["rank_from_end"] = np.where(t.combined < 0, t["rank"], n + 1 - t["rank"])
    elep = t[t.combined < 0].sort_values("rank_from_end")
    mamm = t[t.combined > 0].sort_values("rank_from_end", ascending=False)
    t = pd.concat([elep, mamm])
    t["pct_A_elephant"] = 100 * t.n_more_A_elephant / t.n_targets
    t["pct_A_mammoth"] = 100 * t.n_more_A_mammoth / t.n_targets
    t["pct_concordant"] = np.where(t.combined < 0, t.pct_A_elephant, t.pct_A_mammoth)
    t[["n_targets", "n_more_A_elephant", "n_more_A_mammoth", "pct_A_elephant", "pct_A_mammoth", "pct_concordant",
       "rank_biserial", "FDR_rotation"]].to_csv(HERE / "fig4E_targets.tsv", sep="\t", float_format="%.4g")
    print(t[["n_targets", "n_more_A_elephant", "n_more_A_mammoth", "pct_concordant", "FDR_rotation"]])

    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(4.0, 2.0))
        x = np.arange(len(t))
        ax.bar(x, t.pct_A_elephant, width=0.7, color=f.E_COL, lw=0)
        ax.bar(x, t.pct_A_mammoth, bottom=t.pct_A_elephant, width=0.7, color=f.M_COL, lw=0)
        ax.axhline(50, color="black", lw=0.6, ls="--")
        ax.axvline(len(elep) - 0.5, color="grey", lw=0.5, ls=":")
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 50, 100])
        ax.set_xticks(x)
        ax.set_xlim(-0.6, len(t) - 0.4)
        ax.tick_params(labelsize=7)
        ax.tick_params(axis="x", length=0)
        for sp in ("top", "right", "bottom"):
            ax.spines[sp].set_visible(False)
        if notext:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            ax.set_xticklabels([f"{tf}\nn={int(s.n_targets):,}\nFDR {s.FDR_rotation:.1g}" for tf, s in t.iterrows()],
                               fontsize=5.5)
            for lab in ax.get_xticklabels():
                lab.set_fontstyle("normal")
            ax.set_ylabel("ChEA3 targets (%)", fontsize=8)
            for xi, (_, s) in zip(x, t.iterrows()):
                ax.text(xi, 101, f"{s.pct_concordant:.0f}%", ha="center", va="bottom", fontsize=6)
            ax.text(1.02, 0.25, "More A in\nelephant", transform=ax.transAxes, fontsize=6, color=f.E_COL,
                    va="center")
            ax.text(1.02, 0.75, "More A in\nmammoth", transform=ax.transAxes, fontsize=6, color=f.M_COL,
                    va="center")
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"fig4E_targets{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
