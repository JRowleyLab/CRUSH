#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import gaussian_kde, spearmanr

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
REV = FIG.parent
OUT = FIG / "results"

TIM = Path("Chicken_Study")
GENES_TSV = TIM / "Chicken_Gal6_genes_nochr_gThan5kb_simple.bed"
RNA_TSV = TIM / "RNA_seq" / "Counts" / "AllSample_Counts.txt"
PC1 = Path("hicPCA")
TRACKS = REV / "tracks" / "chicken"
DENS_BED = Path("chicken_nochr_genes.bed")

WINDOW_BP = 10_000
DENS_RES = 50_000
ORDER = ["Erythro", "Fibro", "eRBC", "aRBC"]
NAME = {"Erythro": "Erythroblast", "Fibro": "Fibroblast",
        "eRBC": "eRBC", "aRBC": "aRBC"}
ARMS = ["PC1 50 kb", "CRUSH 50 kb", "CRUSH 1 kb"]
ORIENT = {"PC1 50 kb": False, "CRUSH 50 kb": True, "CRUSH 1 kb": True}


def rank_transform(s):
    return (s.rank(method="average") - 1) / (len(s) - 1) if len(s) > 1 else 0.0


def arm_bw(arm, cell):
    if arm == "PC1 50 kb":
        return PC1 / f"{cell}1.bw"
    res = 50000 if "50" in arm else 1000
    return TRACKS / f"{NAME[cell]}_r{res}" / f"{NAME[cell]}_mergedCrush_{res}.bw"


def load_rna():
    rna = pd.read_csv(RNA_TSV, sep="\t")
    for c in ("Genename", "CellType"):
        rna[c] = rna[c].astype(str).str.strip()
    p = (rna.groupby(["Genename", "CellType"], as_index=False)["pseudo_TPM"]
            .median().rename(columns={"pseudo_TPM": "TPM_median"}))
    p["Expr"] = p.groupby("CellType")["TPM_median"].transform(rank_transform)
    return p


def gene_density():
    g = pd.read_csv(DENS_BED, sep="\t", header=None, usecols=[0, 1, 3],
                    names=["c", "s", "kind"], dtype={"c": str})
    g = g[g["kind"] == "transcript"]
    g["b"] = (g["s"] // DENS_RES).astype(np.int64)
    return g.groupby(["c", "b"]).size().rename("n").reset_index()


def flip_map(bw_path, dens):
    """Chromosomes whose track anti-correlates with gene density."""
    bw = pyBigWig.open(str(bw_path))
    flip = set()
    for c, n in bw.chroms().items():
        key = c.removeprefix("chr")
        sub = dens[dens["c"] == key]
        if len(sub) < 50:
            continue
        v = []
        for b in sub["b"]:
            s = int(b) * DENS_RES
            if s + DENS_RES > n:
                v.append(np.nan); continue
            q = bw.stats(c, s, s + DENS_RES, type="mean")
            v.append(q[0] if q and q[0] is not None else np.nan)
        v = np.asarray(v, float); d = sub["n"].to_numpy(float)
        ok = np.isfinite(v)
        if ok.sum() < 50 or np.std(v[ok]) == 0:
            continue
        r = np.corrcoef(v[ok], d[ok])[0, 1]
        if np.isfinite(r) and r < 0:
            flip.add(key)
    bw.close()
    return flip


def regions(bw_path):
    g = pd.read_csv(GENES_TSV, sep="\t", dtype={"Chr": str, "Strand": str})
    g["Genename"] = g["Genename"].astype(str).str.strip()
    bw = pyBigWig.open(str(bw_path)); sizes = bw.chroms(); bw.close()
    if any(k.startswith("chr") for k in sizes):
        g["Chr"] = "chr" + g["Chr"].str.removeprefix("chr")
    tss = np.where(g["Strand"].to_numpy() == "+",
                   g["Start"].to_numpy(), g["Stop"].to_numpy())
    lo = np.maximum(0, tss - WINDOW_BP)
    hi = tss + WINDOW_BP
    cap = g["Chr"].map(sizes)
    keep = cap.notna()
    out = g.loc[keep, ["Genename", "Chr"]].copy()
    out["WinStart"] = lo[keep.to_numpy()]
    out["WinEnd"] = np.minimum(hi[keep.to_numpy()], cap[keep].to_numpy()).astype(int)
    return out[out["WinEnd"] > out["WinStart"]].reset_index(drop=True)


def score(bw_path, flip):
    reg = regions(bw_path)
    bw = pyBigWig.open(str(bw_path)); sizes = bw.chroms()
    vals = np.full(len(reg), np.nan)
    for i, (c, s, e) in enumerate(zip(reg["Chr"], reg["WinStart"], reg["WinEnd"])):
        if c not in sizes or e > sizes[c]:
            continue
        try:
            q = bw.stats(c, int(s), int(e), type="mean")
        except RuntimeError:
            continue
        if q and q[0] is not None and np.isfinite(q[0]):
            vals[i] = q[0]
    bw.close()
    reg["CompartmentScore"] = vals
    if flip:
        key = reg["Chr"].astype(str).str.removeprefix("chr")
        reg.loc[key.isin(flip), "CompartmentScore"] *= -1
    return reg


LABELS = True


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rna = load_rna()
    dens = gene_density()
    panels, rows = {}, []
    for arm in ARMS:
        for cell in ORDER:
            bw = arm_bw(arm, cell)
            if not bw.exists():
                raise SystemExit(f"missing {bw}")
            fl = flip_map(bw, dens)
            plotted = "oriented" if ORIENT[arm] else "raw"
            for tag, flip in (("oriented", fl), ("raw", None)):
                d = score(bw, flip)
                d["CellType"] = cell
                mm = d.merge(rna, on=["Genename", "CellType"], how="inner") \
                      .dropna(subset=["CompartmentScore", "Expr"])
                mm["CompPct"] = rank_transform(mm["CompartmentScore"])
                rho, p = spearmanr(mm["CompPct"], mm["Expr"])
                rows.append({"arm": arm, "cell_type": cell, "orientation": tag,
                             "n_genes": len(mm), "spearman_rho": rho, "p": p,
                             "chroms_flipped": len(fl) if tag == "oriented" else 0})
                rows[-1]["plotted"] = (tag == plotted)
                if tag == plotted:
                    panels[(arm, cell)] = (mm, rho)
                    note = (f"{len(fl)} chroms flipped" if tag == "oriented"
                            else f"as produced; orienting would flip {len(fl)}")
                    print(f"  {arm:<12}{cell:<9} n={len(mm):>7,} rho={rho:+.3f} "
                          f"[{tag}, {note}]", flush=True)

    q = OUT / "figS4C_rankrank_sourcedata.csv"
    with q.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {q}")

    for labels in (True, False):
        draw(panels, labels)

    df = pd.DataFrame(rows)
    print("\n" + df.pivot_table(index=["arm", "orientation"], columns="cell_type",
                                values="spearman_rho")[ORDER].round(3).to_string())


def draw(panels, labels):
    fig, axes = plt.subplots(len(ARMS), len(ORDER),
                             figsize=(3.5 * len(ORDER), 3.0 * len(ARMS)),
                             squeeze=False, facecolor="white")
    qs = [0, 0.05, 0.10, 0.20, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.99, 1.0]
    cmap = plt.get_cmap("magma")
    for r, arm in enumerate(ARMS):
        for c, cell in enumerate(ORDER):
            ax = axes[r][c]
            mm, rho = panels[(arm, cell)]
            x = mm["CompPct"].to_numpy(); y = mm["Expr"].to_numpy()
            kde = gaussian_kde(np.vstack([x, y]))
            g = np.linspace(0, 1, 100)
            X, Y = np.meshgrid(g, g)
            Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
            Z = (Z - Z.min()) / (Z.max() - Z.min() + 1e-12)
            levs = np.unique(np.quantile(Z, qs))
            ax.contourf(X, Y, Z, levels=levs, cmap=cmap)
            ax.contour(X, Y, Z, levels=levs, colors="k", linewidths=0.5, alpha=0.4)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            if labels:
                ax.set_title((f"{cell}\n" if r == 0 else "")
                             + f"ρ={rho:.3f}, n={len(mm):,}", fontsize=10)
                if r == len(ARMS) - 1:
                    ax.set_xlabel("Compartment percentile")
                if c == 0:
                    tag = "" if ORIENT[arm] else " (raw)"
                    ax.set_ylabel(f"{arm}{tag}\nExpression percentile", fontsize=10)
            else:
                ax.set_xticklabels([]); ax.set_yticklabels([])
    fig.subplots_adjust(right=0.90, wspace=0.28, hspace=0.38)
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.70])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=mpl.colors.Normalize(0, 1))
    if labels:
        cb.set_label("Relative density")
    else:
        cb.set_ticks([])
    stem = "figS4C_rankrank" + ("_labeled" if labels else "")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT / (stem + '.png')} (+ .pdf)")


if __name__ == "__main__":
    main()
