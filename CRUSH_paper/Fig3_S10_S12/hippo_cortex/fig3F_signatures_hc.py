#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata
from sklearn.metrics import roc_auc_score, roc_curve

HERE = Path(__file__).resolve().parent
D = Path("hippoCortex_NewCombinedEigenWay")
MATRIX = D / "hippocortex_25kb_combined_Eigenway_NAed_startRemoved.txt"
STARTS = D / "crush/cortex/crush/hippocortex_25kb_combined_Eigenway.tsv"
UMAP = D / "analysis_old_method/hippo_cortex_oldcalls_25kb_chr1_umap_coords.csv"
GENES = "mm10genes.bed"
RES, TOP, NSPLIT = 25_000, 150, 20
TYPES = ["Neuron", "Mature Oligodendrocyte", "Microglia Etc.", "Adult Astrocyte",
         "Neonatal Astrocyte", "Oligodendrocyte Progenitor"]
COL = {"Neuron": "#d62728", "Neonatal Astrocyte": "#98df8a", "Adult Astrocyte": "#2ca02c",
       "Oligodendrocyte Progenitor": "#e377c2", "Mature Oligodendrocyte": "#9467bd",
       "Microglia Etc.": "#1f77b4"}
SHORT = {"Neuron": "Neuron", "Mature Oligodendrocyte": "Mature oligo.", "Microglia Etc.": "Microglia",
         "Adult Astrocyte": "Adult astro.", "Neonatal Astrocyte": "Neonatal astro.",
         "Oligodendrocyte Progenitor": "OPC"}


def gene_matrix():
    X = pd.read_csv(MATRIX, sep="\t", header=None, na_values=["NA"]).to_numpy(float)
    s = pd.read_csv(STARTS, sep="\t", header=None, usecols=[0])[0].to_numpy()
    assert len(s) == len(X)
    row = {v: i for i, v in enumerate(s)}
    g = pd.read_csv(GENES, sep="\t", header=None, usecols=[0, 1, 2, 3], dtype={0: str})
    g = g[g[0].isin(["1", "chr1"])].drop_duplicates(3)
    names, vals = [], []
    for _, (c, a, b, name) in g.iterrows():
        idx = [row[k] for k in range(a // RES * RES, b, RES) if k in row]
        if idx:
            with np.errstate(all="ignore"):
                v = np.nanmean(X[idx], axis=0)
            if np.isfinite(v).any():
                names.append(name)
                vals.append(np.where(np.isfinite(v), v, np.nanmean(v)))
    return pd.DataFrame(np.array(vals), index=names)


def bh(p):
    p = np.asarray(p)
    q = p * len(p) / rankdata(p)
    q = np.minimum.accumulate(q[np.argsort(p)[::-1]])[::-1]
    out = np.empty_like(q)
    out[np.argsort(p)] = q
    return np.minimum(out, 1)


def signature(G, y):
    """Top genes more A in the y==1 cells (spec = 2*AUC - 1, then FDR)."""
    rows = []
    for gene, v in zip(G.index, G.to_numpy()):
        if v.std() == 0:
            continue
        rows.append((gene, 2 * roc_auc_score(y, v) - 1, mannwhitneyu(v[y == 1], v[y == 0]).pvalue))
    r = pd.DataFrame(rows, columns=["gene", "spec", "p"])
    r["FDR"] = bh(r.p)
    return r[r.spec > 0].sort_values(["spec", "FDR"], ascending=[False, True]).head(TOP)


def main():
    U = pd.read_csv(UMAP)
    G = gene_matrix()
    assert G.shape[1] == len(U)
    print(f"{len(G)} chr1 genes x {G.shape[1]} cells")
    typed = (U.Group != "Unknown").to_numpy()
    Gt, grp = G.loc[:, typed], U.Group.to_numpy()[typed]
    rng = np.random.default_rng(1)

    score, roc, stats, sigs = {}, {}, [], []
    for t in TYPES:
        y = (grp == t).astype(int)
        sig = signature(Gt, y)
        sigs.append(sig.assign(cell_type=t))
        score[t] = G.loc[sig.gene].mean(axis=0).to_numpy()
        score[t] = (score[t] - score[t].mean()) / score[t].std()
        st = score[t][typed]
        roc[t] = roc_curve(y, st)[:2]
        held = []
        for _ in range(NSPLIT):
            tr = np.zeros(len(y), bool)
            for k in (0, 1):
                ix = np.flatnonzero(y == k)
                tr[rng.choice(ix, len(ix) // 2, replace=False)] = True
            sg = signature(Gt.loc[:, tr], y[tr]).gene
            held.append(roc_auc_score(y[~tr], Gt.loc[sg, ~tr].mean(axis=0)))
        stats.append({"cell_type": t, "n_cells": int(y.sum()), "n_other": int((1 - y).sum()),
                      "AUC_in_sample": roc_auc_score(y, st), "AUC_heldout_mean": np.mean(held),
                      "AUC_heldout_sd": np.std(held, ddof=1),
                      "min_spec_in_signature": sig.spec.min(), "max_FDR_in_signature": sig.FDR.max()})
        print(stats[-1], flush=True)
    pd.DataFrame(stats).to_csv(HERE / "fig3F_hc_stats.tsv", sep="\t", index=False, float_format="%.4g")
    pd.concat(sigs)[["cell_type", "gene", "spec", "p", "FDR"]].to_csv(
        HERE / "fig3F_hc_signatures.tsv", sep="\t", index=False, float_format="%.4g")
    auc = {s["cell_type"]: s["AUC_in_sample"] for s in stats}

    for notext in (False, True):
        fig = plt.figure(figsize=(10.5, 4.6))
        gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 0.35, 1.4], wspace=0.3, hspace=0.3)
        for k, t in enumerate(TYPES):
            ax = fig.add_subplot(gs[k // 3, k % 3])
            v = score[t]
            o = np.argsort(v)
            sc = ax.scatter(U.UMAP_1.to_numpy()[o], U.UMAP_2.to_numpy()[o], c=v[o], s=2, lw=0,
                            cmap="BrBG_r", vmin=-2, vmax=2)
            ax.set_xticks([])
            ax.set_yticks([])
            cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
            cb.ax.tick_params(labelsize=5)
            if not notext:
                ax.set_title(f"{SHORT[t]} (top {TOP} A genes)", fontsize=7)
        ax = fig.add_subplot(gs[:, 4])
        for t in TYPES:
            fpr, tpr = roc[t]
            ax.plot(fpr, tpr, color=COL[t], lw=1.2, label=f"{SHORT[t]} ({auc[t]:.2f})")
        ax.plot([0, 1], [0, 1], ls="--", color="grey", lw=0.8)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if not notext:
            ax.set_xlabel("FPR", fontsize=8)
            ax.set_ylabel("TPR", fontsize=8)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"fig3F_hc{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
