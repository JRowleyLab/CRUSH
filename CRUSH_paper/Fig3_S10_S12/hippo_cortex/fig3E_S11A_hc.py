#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
RUN = Path("DipC/"
           "hippoCortex_NewCombinedEigenWay/analysis_old_method")
GROUPS = {"Neuron": "#d62728", "Neonatal Astrocyte": "#98df8a", "Adult Astrocyte": "#2ca02c",
          "Oligodendrocyte Progenitor": "#e377c2", "Mature Oligodendrocyte": "#9467bd",
          "Microglia Etc.": "#1f77b4", "Unknown": "#c7c7c7"}
REGIONS = {"Cortex": "#e7298a", "Hippocampus": "#1b9e77"}


def umap(d, col, colors, stem, notext):
    fig, ax = plt.subplots(figsize=(2.8, 2.6))
    for k in reversed(list(colors)):
        s = d[d[col] == k]
        ax.scatter(s.UMAP_1, s.UMAP_2, s=3, color=colors[k], lw=0, label=f"{k} ({len(s)})")
    ax.set_xticks([])
    ax.set_yticks([])
    if not notext:
        ax.set_xlabel("UMAP 1", fontsize=8)
        ax.set_ylabel("UMAP 2", fontsize=8)
        h, l = ax.get_legend_handles_labels()
        ax.legend(h[::-1], l[::-1], frameon=False, fontsize=6, markerscale=2.5,
                  loc="upper left", bbox_to_anchor=(1.0, 1.0))
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    d = pd.read_csv(RUN / "hippo_cortex_oldcalls_25kb_chr1_umap_coords.csv")
    print(pd.read_csv(RUN / "hippo_cortex_oldcalls_purity_summary.csv").round(3).T.to_string())
    for notext in (False, True):
        sfx = "_notext" if notext else ""
        umap(d, "Group", GROUPS, "fig3E_hc" + sfx, notext)
        umap(d, "Region", REGIONS, "figS11A_hc_region" + sfx, notext)


if __name__ == "__main__":
    main()
