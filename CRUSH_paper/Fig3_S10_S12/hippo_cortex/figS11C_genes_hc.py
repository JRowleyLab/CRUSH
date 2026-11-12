#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

import fig3F_signatures_hc as f
from fig3G_genes_hc import stars

HERE = Path(__file__).resolve().parent
GENES = {"Hhat": "Neuron", "Brinp3": "Mature Oligodendrocyte", "Slc40a1": "Microglia Etc.",
         "Prex2": "Adult Astrocyte", "Gli2": "Neonatal Astrocyte", "Sntg1": "Oligodendrocyte Progenitor"}


def main():
    U = pd.read_csv(f.UMAP)
    G = f.gene_matrix()
    typed = (U.Group != "Unknown").to_numpy()
    grp = U.Group.to_numpy()
    stats = []
    for notext in (False, True):
        fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), gridspec_kw={"wspace": 0.35, "hspace": 0.6})
        for ax, (gene, t) in zip(axes.flat, GENES.items()):
            v = G.loc[gene].to_numpy()
            v = (v - v.mean()) / v.std()
            y = (grp[typed] == t).astype(int)
            vt = v[typed]
            auc = roc_auc_score(y, vt)
            p = mannwhitneyu(vt[y == 1], vt[y == 0], alternative="greater").pvalue
            if not notext:
                stats.append({"gene": gene, "cell_type": t, "n": int(y.sum()), "AUC": auc,
                              "median_z_type": np.median(vt[y == 1]), "median_z_rest": np.median(vt[y == 0]),
                              "p_wilcoxon_greater": p})
            data = [v[grp == k] for k in f.TYPES]
            bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.6,
                            medianprops=dict(color="black", lw=0.8))
            for b, k in zip(bp["boxes"], f.TYPES):
                b.set_facecolor(f.COL[k] if k == t else "#e0e0e0")
                b.set_linewidth(0.6)
            ax.axhline(0, color="grey", lw=0.5, ls=":")
            ax.set_xticks(range(1, len(f.TYPES) + 1))
            i = f.TYPES.index(t) + 1
            top = max(np.percentile(d, 75) + 1.5 * (np.percentile(d, 75) - np.percentile(d, 25)) for d in data)
            ax.text(i, top, stars(p), ha="center", va="bottom", fontsize=8)
            ax.set_ylim(top=top + 0.7)
            ax.tick_params(labelsize=6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if notext:
                ax.set_xticklabels([])
            else:
                ax.set_title(gene, fontsize=8, style="italic", loc="left")
                ax.set_title(f"AUC = {auc:.2f}", fontsize=7, loc="right")
                ax.set_xticklabels([f.SHORT[k] for k in f.TYPES], rotation=40, ha="right", fontsize=6)
        if not notext:
            for ax in axes[:, 0]:
                ax.set_ylabel("CRUSH (z-score)", fontsize=7)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"figS11C_hc{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)
    st = pd.DataFrame(stats)
    st.to_csv(HERE / "figS11C_hc_stats.tsv", sep="\t", index=False, float_format="%.4g")
    print(st.to_string(index=False))


if __name__ == "__main__":
    main()
