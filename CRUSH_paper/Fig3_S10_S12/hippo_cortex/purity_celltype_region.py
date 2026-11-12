#!/usr/bin/env python3
import re
import unicodedata
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = Path("hippoCortex_NewCombinedEigenWay")
MATRICES = {"25kb": "hippocortex_25kb_combined_Eigenway_NAed_startRemoved.txt",
            "50kb": "hippocortex_50kb_combined_Eigenway_NAed_startRemoved.txt",
            "100kb": "hippocortex_100kb_combined_Eigenway_NAed_startRemoved.txt"}
SEEDS, K = range(10), 15
BAR_COLORS = {"25kb": "#1f77b4", "50kb": "#2ca02c", "100kb": "#17becf"}


def clean(v):
    if pd.isna(v):
        return "Unknown"
    v = re.sub(r"\s+", " ", unicodedata.normalize("NFKC", str(v))).strip()
    return "Unknown" if v.casefold() in {"", "na", "n/a", "nan", "none", "null", "unknown"} else v


def neurons_merged(t):
    if t == "Unknown" or any(k in t for k in ("Astrocyte", "Oligodendrocyte", "Microglia")):
        return t
    return "Neuron"


def labels():
    reg = pd.read_csv(DATA / "updated_metadata_cells_contacts_cortex_hippo_names.txt", sep="\t",
                      header=None, usecols=[0, 1], dtype=str, names=["UniqueID", "Region"])
    reg["UniqueID"] = reg.UniqueID.str.strip()
    reg["Region"] = reg.Region.str.strip().str.lower().map({"cortex": "Cortex", "hippo": "Hippocampus",
                                                           "hippocampus": "Hippocampus"})
    ct = pd.read_csv(DATA / "metadata.cells_contacts_100k.txt", sep="\t", header=None, usecols=[0, 1],
                     dtype=str, names=["UniqueID", "CellType"])
    ct["UniqueID"] = ct.UniqueID.str.strip()
    ct["CellType"] = ct.CellType.apply(clean)
    lab = reg.merge(ct, on="UniqueID", how="left")
    lab["Group"] = lab.CellType.fillna("Unknown").map(neurons_merged)
    return lab


def purity(e, y):
    from scipy.spatial.distance import cdist
    D = cdist(e, e)
    np.fill_diagonal(D, np.inf)
    return (y[np.argsort(D, axis=1)[:, :K]] == y[:, None]).mean(axis=1).mean()


def run(res):
    import umap.umap_ as umap
    from sklearn.impute import KNNImputer
    from sklearn.preprocessing import StandardScaler
    lab = labels()
    X = pd.read_csv(DATA / MATRICES[res], sep="\t", header=None, na_values=["NA"]).T.values.astype(float)
    assert len(X) == len(lab)
    Xs = StandardScaler().fit_transform(KNNImputer(n_neighbors=5).fit_transform(X))
    keep = (lab.Group != "Unknown").to_numpy()
    g, r = lab.Group.to_numpy(str)[keep], lab.Region.to_numpy(str)[keep]
    rows = []
    for seed in SEEDS:
        e = umap.UMAP(n_neighbors=15, min_dist=0.01, metric="cosine", random_state=seed).fit_transform(Xs)[keep]
        rows.append({"resolution": res, "seed": seed, "celltype": purity(e, g), "region": purity(e, r)})
        print(rows[-1], flush=True)
    return rows


def main():
    out = HERE / "purity_celltype_region.csv"
    if out.exists():
        t = pd.read_csv(out)
    else:
        with ProcessPoolExecutor(3) as ex:
            rows = [r for part in ex.map(run, MATRICES) for r in part]
        t = pd.DataFrame(rows)
        t.to_csv(out, index=False)
    lab = labels()
    lab = lab[lab.Group != "Unknown"]
    chance = {k: float((lab[c].value_counts(normalize=True) ** 2).sum())
              for k, c in (("celltype", "Group"), ("region", "Region"))}
    s = t.groupby("resolution", sort=False)[["celltype", "region"]].agg(["mean", "std"])
    print(s.round(3).to_string(), "\nchance", {k: round(v, 3) for k, v in chance.items()})

    res = list(MATRICES)
    x = np.arange(len(res))
    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(4.5, 3.2))
        for off, key, hatch in ((-0.2, "celltype", None), (0.2, "region", "////")):
            m, sd = s[(key, "mean")].to_numpy(), s[(key, "std")].to_numpy()
            ax.bar(x + off, m, 0.38, yerr=sd, capsize=3, color=[BAR_COLORS[r] for r in res],
                   edgecolor="black", linewidth=0.6, hatch=hatch, alpha=1 if hatch is None else 0.45)
            ax.hlines([chance[key]] * len(x), x + off - 0.19, x + off + 0.19, ls="--", color="black", lw=0.8)
        ax.set_xticks(x)
        ax.set_ylim(0.4, 1)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if notext:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(res)
            ax.set_ylabel(f"Neighborhood purity (k={K})")
            ax.set_xlabel("CRUSH")
            from matplotlib.patches import Patch
            ax.legend([Patch(facecolor="white", edgecolor="black"),
                       Patch(facecolor="white", edgecolor="black", hatch="////")],
                      ["Cell type", "Region (cortex vs hippocampus)"], frameon=False, fontsize=7,
                      loc="upper right")
            ax.text(x[-1] + 0.42, chance["region"], "chance", fontsize=6, color="black", va="center")
        fig.tight_layout()
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"figS11B_purity{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
