#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CA = Path("DipC/"
          "mouseOlfactory_NewCombinedEigenWay/compAnalysis")
UMAP = CA / "UMAP_projection_scores_Crush25kb_c1_noBulk.csv"
BED = CA / "mouseOLF_crush_25kb_nobulk_NAed_genesIntersected.bed"
SPEC = CA / "AUC_outputs_noBulk_re_rank/CRUSH_25kb_all_genes_specificity_lenient.csv"
GENE_COL, COMP_START, TOP = 3, 7, 150
GENES = ["Pdc", "Lgr6"]


def z(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def umap_plot(U, v, stem, title, notext):
    fig, ax = plt.subplots(figsize=(2.6, 2.4))
    o = np.argsort(v)
    sc = ax.scatter(U.UMAP_1.to_numpy()[o], U.UMAP_2.to_numpy()[o], c=v[o], s=4, lw=0,
                    cmap="BrBG_r", vmin=-2, vmax=2)
    ax.set_xticks([])
    ax.set_yticks([])
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=6)
    if not notext:
        ax.set_title(title, fontsize=8)
        cb.set_label("CRUSH (z-score)", fontsize=7)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"{stem}{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    U = pd.read_csv(UMAP)
    n = len(U)
    B = pd.read_csv(BED, sep="\t", header=None)
    genes = B.iloc[:, GENE_COL].astype(str)

    spec = pd.read_csv(SPEC)
    pos = spec.query("spec_score > 0").sort_values(["spec_score", "FDR"], ascending=[False, True]).head(TOP).gene
    neg = spec.query("spec_score < 0").sort_values(["spec_score", "FDR"], ascending=[True, True]).head(TOP).gene
    last = B.iloc[:, -n:].apply(pd.to_numeric, errors="coerce").to_numpy()
    sets = {"retina": pos, "MOE": neg}
    for notext in (False, True):
        for name, gs in sets.items():
            m = genes.isin(gs.astype(str)).to_numpy()
            v = z(np.nanmean(last[m], axis=0))
            if not notext:
                print(f"3C {name}-like: {m.sum()} rows from {gs.size} genes")
            umap_plot(U, v, f"fig3C_z_{name}", f"Top {TOP} {name} A genes", notext)

    vals = B.iloc[:, COMP_START:COMP_START + n].apply(pd.to_numeric, errors="coerce")
    vals["gene"] = genes
    G = vals.groupby("gene").mean()
    for g in GENES:
        x = G.loc[g].to_numpy(float)
        x = np.where(np.isnan(x), np.nanmean(x), x)
        for notext in (False, True):
            umap_plot(U, z(x), f"fig3D_z_{g}", g, notext)


if __name__ == "__main__":
    main()
