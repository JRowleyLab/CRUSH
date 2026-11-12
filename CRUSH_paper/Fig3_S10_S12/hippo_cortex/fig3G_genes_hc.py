#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score, roc_curve

import fig3F_signatures_hc as f

HERE = Path(__file__).resolve().parent
GENES = {"Kcnh1": "Neuron", "St18": "Mature Oligodendrocyte", "Dock10": "Microglia Etc.",
         "Gli2": "Neonatal Astrocyte"}
LIM = 2


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def main():
    U = pd.read_csv(f.UMAP)
    G = f.gene_matrix()
    typed = (U.Group != "Unknown").to_numpy()
    grp = U.Group.to_numpy()
    stats = []
    for notext in (False, True):
        fig, axes = plt.subplots(len(GENES), 3, figsize=(7.2, 2.5 * len(GENES)),
                                 gridspec_kw={"width_ratios": [1, 1.5, 0.9], "wspace": 0.6, "hspace": 0.75})
        for row, (gene, t) in zip(axes, GENES.items()):
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
            ax = row[0]
            o = np.argsort(v)
            sc = ax.scatter(U.UMAP_1.to_numpy()[o], U.UMAP_2.to_numpy()[o], c=v[o], s=2, lw=0,
                            cmap="BrBG_r", vmin=-LIM, vmax=LIM)
            ax.set_xticks([])
            ax.set_yticks([])
            cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
            cb.ax.tick_params(labelsize=5)

            ax = row[1]
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

            ax = row[2]
            fpr, tpr = roc_curve(y, vt)[:2]
            ax.plot(fpr, tpr, color=f.COL[t], lw=1.2)
            ax.plot([0, 1], [0, 1], ls="--", color="grey", lw=0.6)
            ax.set_aspect("equal")
            for a in row[1:]:
                a.tick_params(labelsize=6)
                for s in ("top", "right"):
                    a.spines[s].set_visible(False)
            if notext:
                row[1].set_xticklabels([])
            else:
                row[0].set_title(gene, fontsize=8, style="italic")
                row[1].set_xticklabels([f.SHORT[k] for k in f.TYPES], rotation=40, ha="right", fontsize=6)
                row[1].set_ylabel("CRUSH (z-score)", fontsize=7)
                row[2].set_title(f"AUC = {auc:.2f}", fontsize=7)
                row[2].set_xlabel("FPR", fontsize=7)
                row[2].set_ylabel("TPR", fontsize=7)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"fig3G_hc{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)
    st = pd.DataFrame(stats)
    st.to_csv(HERE / "fig3G_hc_stats.tsv", sep="\t", index=False, float_format="%.4g")
    print(st.to_string(index=False))


if __name__ == "__main__":
    main()
