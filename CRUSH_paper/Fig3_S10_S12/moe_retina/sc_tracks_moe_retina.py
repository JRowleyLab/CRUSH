#!/usr/bin/env python3
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, optimal_leaf_ordering
from scipy.spatial.distance import pdist

HERE = Path(__file__).resolve().parent
CA = Path("DipC/"
          "mouseOlfactory_NewCombinedEigenWay/compAnalysis")
MATRIX = CA / "mouseOLF_crush_25kb_nobulk_NAed.bed"
UMAP = CA / "UMAP_projection_scores_Crush25kb_c1_noBulk.csv"
GENES = "mm10genes.bed"
REGION = (120.0, 125.0)
N_PER, CLIP, NPC = 25, 300, 10
A_COL, B_COL = "#8c510a", "#01665e"
TISSUE_COL = {"MOE": "#e41a1c", "retina": "#377eb8"}
SHOW_GENES = ["Steap3", "En1", "Insig2", "Dpp10", "Pdc", "Lgr6"]


def load():
    B = pd.read_csv(MATRIX, sep="\t", header=None, na_values=["NA"])
    U = pd.read_csv(UMAP)
    X = B.iloc[:, 3:].to_numpy(float)
    assert X.shape[1] == len(U)
    return B[1].to_numpy(), X, U.TissueType.to_numpy()


def impute(X):
    return np.where(np.isnan(X), np.nanmean(X, axis=1, keepdims=True), X)


def order_score(X, tissue):
    Xi = impute(X)
    r = lambda m: np.array([np.corrcoef(Xi[:, j], m)[0, 1] for j in range(Xi.shape[1])])
    return r(Xi[:, tissue == "retina"].mean(1)) - r(Xi[:, tissue == "MOE"].mean(1))


def pick(score, tissue):
    idx = []
    for t in ("MOE", "retina"):
        ix = np.flatnonzero(tissue == t)
        ix = ix[np.argsort(score[ix])]
        idx += list(ix[np.round(np.linspace(0, len(ix) - 1, N_PER)).astype(int)])
    return np.array(sorted(idx, key=lambda j: score[j]))


def cluster(X, sel, score):
    """Linkage of the shown cells and their leaf order, MOE-like end first."""
    P = PCA(NPC, random_state=0).fit_transform(
        StandardScaler().fit_transform(KNNImputer(n_neighbors=5).fit_transform(X.T)))[sel]
    Z = optimal_leaf_ordering(linkage(P, "ward"), pdist(P))
    leaves = dendrogram(Z, no_plot=True)["leaves"]
    if score[sel[leaves[: len(sel) // 2]]].mean() > score[sel[leaves[len(sel) // 2:]]].mean():
        return Z, sel[leaves[::-1]], True
    return Z, sel[leaves], False


def fill(ax, x, v, clip=CLIP):
    v = np.clip(v, -clip, clip)
    ax.fill_between(x, 0, v, where=v > 0, color=A_COL, lw=0, step="post")
    ax.fill_between(x, 0, v, where=v < 0, color=B_COL, lw=0, step="post")
    ax.set_ylim(-clip, clip)
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def genes_in(a, b):
    g = pd.read_csv(GENES, sep="\t", header=None, usecols=[0, 1, 2, 3, 5], dtype={0: str})
    g = g[g[0].isin(["1", "chr1"]) & (g[2] > a) & (g[1] < b)].drop_duplicates(3)
    return g.sort_values(1)


def gene_track(ax, g, a, b, notext):
    ends = []
    for _, (c, s, e, name, strand) in g.iterrows():
        row = next((k for k, x in enumerate(ends) if x < s - 0.02 * (b - a)), len(ends))
        if row == len(ends):
            ends.append(e)
        else:
            ends[row] = e
        lab = name in SHOW_GENES
        ax.plot([s / 1e6, e / 1e6], [-row, -row], color="black" if lab else "grey", lw=2.5 if lab else 1.5,
                solid_capstyle="butt")
        if lab and not notext:
            ax.text((s + e) / 2e6, -row - 0.45, name, ha="center", va="top", fontsize=6, style="italic",
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.5), zorder=3)
    ax.set_ylim(-max(len(ends), 1) - 0.6, 0.6)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


def main():
    args = [v for v in sys.argv[1:] if v != "--cluster"]
    clustered = "--cluster" in sys.argv
    a, b = (float(args[0]), float(args[1])) if len(args) == 2 else REGION
    tag = f"chr1_{a:g}-{b:g}Mb"
    pos, X, tissue = load()
    score = order_score(X, tissue)
    sel = pick(score, tissue)
    suffix = f"_{len(sel)}cells"
    if clustered:
        Z, sel, flipped = cluster(X, sel, score)
        suffix += "_dendro"
    w = (pos >= a * 1e6) & (pos < b * 1e6)
    x = np.r_[pos[w], pos[w][-1] + 25_000] / 1e6
    V = X[w]
    pd.DataFrame({"rank": range(1, len(sel) + 1), "cell": [f"cell_{j + 1:03d}" for j in sel],
                  "tissue": tissue[sel], "retina_minus_MOE_corr": score[sel]}).to_csv(
        HERE / f"sc_tracks_{tag}{suffix}.tsv", sep="\t", index=False, float_format="%.4f")
    g = genes_in(a * 1e6, b * 1e6)
    pb = {t: np.nanmean(V[:, tissue == t], axis=1) for t in ("MOE", "retina")}
    pb_lim = np.ceil(max(np.abs(v).max() for v in pb.values()) / 10) * 10

    for notext in (False, True):
        n = len(sel)
        fig = plt.figure(figsize=(6.5, 0.22 * n + 2.2))
        gs = fig.add_gridspec(n + 4, 2, width_ratios=[0.015, 1], wspace=0.01, hspace=0.15,
                              height_ratios=[1.6, 1.6, 0.6] + [1] * n + [2.5])
        cell_axes = []
        rows = [("MOE", pb["MOE"]), ("retina", pb["retina"]), None] + \
               [(tissue[j], V[:, j]) for j in sel] + [None]
        for k, item in enumerate(rows):
            if item is None:
                continue
            t, v = item
            ax = fig.add_subplot(gs[k, 1])
            fill(ax, x, np.r_[v, v[-1]], pb_lim if k < 2 else CLIP)
            ax.set_xlim(a, b)
            ax.set_xticks([])
            ax.axhline(0, color="black", lw=0.3)
            if 2 < k < n + 3:
                cell_axes.append(ax)
            s = fig.add_subplot(gs[k, 0])
            s.axis("off")
            s.add_patch(plt.Rectangle((0, 0), 1, 1, color=TISSUE_COL[t], transform=s.transAxes))
            if not notext and k < 2:
                ax.text(1.01, 0.5, f"\u00b1{pb_lim:g}", transform=ax.transAxes, ha="left", va="center",
                        fontsize=5)
                ax.text(-0.04, 0.5, f"{t} (all cells)", transform=ax.transAxes, ha="right", va="center",
                        fontsize=6)
        if clustered:
            top, bot = cell_axes[0].get_position(), cell_axes[-1].get_position()
            yc0, yc1 = (top.y0 + top.y1) / 2, (bot.y0 + bot.y1) / 2
            dax = fig.add_axes([top.x0 - 0.13, yc1, 0.11, yc0 - yc1])
            dendrogram(Z, orientation="left", ax=dax, color_threshold=0, above_threshold_color="black",
                       no_labels=True)
            for ln in dax.collections:
                ln.set_linewidth(0.5)
            dax.set_ylim((10 * n - 5, 5) if not flipped else (5, 10 * n - 5))
            dax.axis("off")
        ax = fig.add_subplot(gs[-1, 1])
        gene_track(ax, g, a * 1e6, b * 1e6, notext)
        ax.set_xlim(a, b)
        ax.tick_params(labelsize=6)
        if not notext:
            ax.set_xlabel("chr1 (Mb)", fontsize=7)
            fig.text(0.92, 0.45, f"\u00b1{CLIP}", fontsize=5, ha="left")
            fig.text(0.08 if not clustered else 0.0, 0.5,
                     f"{n} single cells ({'clustered' if clustered else 'MOE-like → retina-like'})", rotation=90,
                     ha="center", va="center", fontsize=7)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"sc_tracks_{tag}{suffix}{'_notext' if notext else ''}.{ext}", dpi=300,
                        bbox_inches="tight")
        plt.close(fig)

        o = np.argsort(score)
        fig, (s, ax) = plt.subplots(1, 2, figsize=(6.5, 4), gridspec_kw={"width_ratios": [0.02, 1],
                                                                       "wspace": 0.01})
        im = ax.imshow(np.clip(V[:, o].T, -CLIP, CLIP), aspect="auto", cmap="BrBG_r", vmin=-CLIP, vmax=CLIP,
                       extent=(x[0], x[-1], len(o), 0), interpolation="nearest")
        ax.set_yticks([])
        ax.tick_params(labelsize=6)
        s.imshow(np.array([[matplotlib.colors.to_rgb(TISSUE_COL[t])] for t in tissue[o]]), aspect="auto",
                 interpolation="nearest")
        s.axis("off")
        cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        cb.ax.tick_params(labelsize=6)
        if not notext:
            ax.set_xlabel("chr1 (Mb)", fontsize=7)
            cb.set_label("CRUSH", fontsize=7)
            s.set_title(f"{len(o)} cells", fontsize=6)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"sc_heatmap_{tag}{'_notext' if notext else ''}.{ext}", dpi=300,
                        bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
