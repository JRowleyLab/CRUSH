#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
LINES = ["A549", "HCT116", "HeLa-S3", "HepG2", "HUVEC", "IMR90", "K562"]
SETS = {"enhancer": "enhancerDB", "active enhancer": "Active", "active enhancer, distal": "Active, distal"}
COL = ["#c6dbef", "#6baed6", "#08519c"]
W = 0.26


def main():
    t = pd.read_csv(HERE / "enh_predict.tsv", sep="\t")
    v = t.pivot(index="line", columns="label", values="AUROC").loc[LINES, list(SETS)]
    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(3.6, 2.2))
        for k, (lab, name) in enumerate(SETS.items()):
            x = [i + (k - 1) * W for i in range(len(v))]
            ax.bar(x, v[lab], width=W, color=COL[k], edgecolor="black", lw=0.4, label=name)
            if not notext:
                for xi, a in zip(x, v[lab]):
                    ax.text(xi, a + 0.01, f"{a:.2f}", ha="center", va="bottom", fontsize=4, rotation=90)
        ax.axhline(0.5, ls="--", color="grey", lw=0.8)
        ax.set_ylim(0, 1)
        ax.set_xticks(range(len(v)))
        ax.tick_params(labelsize=7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if notext:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            ax.set_xticklabels(LINES, rotation=40, ha="right")
            ax.set_ylabel("AUROC", fontsize=8)
            ax.legend(fontsize=6, frameon=False, loc="lower left", bbox_to_anchor=(1, 0))
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"active_enh_auroc_bar{'_notext' if notext else ''}.{ext}", dpi=300,
                        bbox_inches="tight")
        plt.close(fig)
    print(v.round(3).to_string())


if __name__ == "__main__":
    main()
