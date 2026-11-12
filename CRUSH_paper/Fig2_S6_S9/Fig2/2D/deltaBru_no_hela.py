#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

CELL_LINES = ["hct", "hepg2", "imr90", "k562"]
CELL_COLORS = {
    "hct": "tab:red",
    "hepg2": "tab:green",
    "imr90": "tab:purple",
    "k562": "tab:orange",
}

def row_center(x):
    return x - np.nanmean(x, axis=1, keepdims=True)

def col_zscore(x):
    mu = np.nanmean(x, axis=0, keepdims=True)
    sd = np.nanstd(x, axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    return (x - mu) / sd

def finite_median(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return np.nan if x.size == 0 else float(np.median(x))

def plot_summary(unchanged, changed, med_unch, med_ch, ylabel, title, outpath, showfliers=False, alpha=0.9):
    plt.figure(figsize=(5, 4), dpi=200)
    plt.boxplot(
        [unchanged, changed],
        labels=["Unchanged", "Changed"],
        showfliers=showfliers,
        medianprops=dict(color="black", linewidth=1.5),
        boxprops=dict(color="black"),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
    )
    plt.axhline(0, linestyle="--", linewidth=1, color="black")

    offsets = np.linspace(-0.1, 0.1, len(CELL_LINES))
    for off, cl in zip(offsets, CELL_LINES):
        if np.isfinite(med_unch[cl]):
            plt.scatter(
                1 + off, med_unch[cl], color=CELL_COLORS[cl],
                s=50, alpha=alpha, edgecolor="black", linewidth=0.4,
                label=cl.upper()
            )
        if np.isfinite(med_ch[cl]):
            plt.scatter(
                2 + off, med_ch[cl], color=CELL_COLORS[cl],
                s=50, alpha=alpha, edgecolor="black", linewidth=0.4
            )

    handles, labels = plt.gca().get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    plt.legend(uniq.values(), uniq.keys(), frameon=False, title="Cell line")

    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()
    print(f"Saved: {outpath}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--out_prefix", default="summary_changed_vs_unchanged")
    ap.add_argument("--bru_pct", type=float, default=95.0)
    ap.add_argument("--B_max", type=float, default=0.0)
    ap.add_argument("--log1p", action="store_true")
    ap.add_argument("--col_zscore", action="store_true")
    ap.add_argument("--center_rows", action="store_true")
    ap.add_argument("--showfliers", action="store_true")
    ap.add_argument("--alpha_points", type=float, default=0.9)
    args = ap.parse_args()

    cols = [
        "chr","start","end",
        "hct_comp","hela_comp","hepg2_comp","imr90_comp","k562_comp",
        "chr2","start2","end2",
        "hct_bru","hela_bru","hepg2_bru","imr90_bru","k562_bru"
    ]
    df = pd.read_csv(args.infile, sep=args.sep, header=None, names=cols)

    comp_cols = ["hct_comp", "hepg2_comp", "imr90_comp", "k562_comp"]
    bru_cols = ["hct_bru", "hepg2_bru", "imr90_bru", "k562_bru"]

    for c in comp_cols + bru_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    comp = df[comp_cols].to_numpy(dtype=float)
    bru = df[bru_cols].to_numpy(dtype=float)

    if args.log1p:
        bru = np.log1p(bru)
    if args.col_zscore:
        bru = col_zscore(bru)
    if args.center_rows:
        bru = row_center(bru)

    bru_thresh = {}
    for j, cl in enumerate(CELL_LINES):
        x = bru[:, j]
        x = x[np.isfinite(x)]
        bru_thresh[cl] = float(np.percentile(x, args.bru_pct))

    print("\nPer-cell Bru thresholds:")
    for cl in CELL_LINES:
        print(f"{cl.upper():5s}: p{args.bru_pct:g} = {bru_thresh[cl]:.4f}")

    all_bru_changed, all_bru_unchanged = [], []
    all_comp_changed, all_comp_unchanged = [], []
    med_bru_changed, med_bru_unchanged = {}, {}
    med_comp_changed, med_comp_unchanged = {}, {}

    for j, cl in enumerate(CELL_LINES):
        changed = (bru[:, j] >= bru_thresh[cl]) & (np.all(np.delete(bru, j, axis=1) <= args.B_max, axis=1))
        unchanged = ~changed

        delta_bru = bru[:, j] - np.mean(np.delete(bru, j, axis=1), axis=1)
        delta_comp = comp[:, j] - np.mean(np.delete(comp, j, axis=1), axis=1)

        db_c, db_u = delta_bru[changed], delta_bru[unchanged]
        dc_c, dc_u = delta_comp[changed], delta_comp[unchanged]

        all_bru_changed.append(db_c)
        all_bru_unchanged.append(db_u)
        all_comp_changed.append(dc_c)
        all_comp_unchanged.append(dc_u)

        med_bru_changed[cl] = finite_median(db_c)
        med_bru_unchanged[cl] = finite_median(db_u)
        med_comp_changed[cl] = finite_median(dc_c)
        med_comp_unchanged[cl] = finite_median(dc_u)

        print(
            f"{cl.upper():5s}: "
            f"n_changed={changed.sum():7d}, n_unchanged={unchanged.sum():7d}, "
            f"ΔBru med={med_bru_unchanged[cl]:.4f}/{med_bru_changed[cl]:.4f}, "
            f"ΔComp med={med_comp_unchanged[cl]:.4f}/{med_comp_changed[cl]:.4f}"
        )

    bru_c = np.concatenate(all_bru_changed)
    bru_u = np.concatenate(all_bru_unchanged)
    comp_c = np.concatenate(all_comp_changed)
    comp_u = np.concatenate(all_comp_unchanged)

    bru_c = bru_c[np.isfinite(bru_c)]
    bru_u = bru_u[np.isfinite(bru_u)]
    comp_c = comp_c[np.isfinite(comp_c)]
    comp_u = comp_u[np.isfinite(comp_u)]

    _, p_bru = mannwhitneyu(bru_c, bru_u, alternative="greater")
    _, p_comp = mannwhitneyu(comp_c, comp_u, alternative="greater")

    print("\nPooled tests:")
    print(f"ΔBru  p={p_bru:.3e}, median(changed)={np.median(bru_c):.4f}, median(unchanged)={np.median(bru_u):.4f}")
    print(f"ΔCRUSH p={p_comp:.3e}, median(changed)={np.median(comp_c):.4f}, median(unchanged)={np.median(comp_u):.4f}")

    plot_summary(
        bru_u, bru_c,
        med_bru_unchanged, med_bru_changed,
        ylabel="ΔBru (focal − mean others)",
        title=f"ΔBru: Changed vs Unchanged (p={p_bru:.2e})",
        outpath=f"{args.out_prefix}.deltaBru.summary.svg",
        showfliers=args.showfliers,
        alpha=args.alpha_points,
    )

    plot_summary(
        comp_u, comp_c,
        med_comp_unchanged, med_comp_changed,
        ylabel="ΔCRUSH (focal − mean others)",
        title=f"ΔCRUSH: Changed vs Unchanged (p={p_comp:.2e})",
        outpath=f"{args.out_prefix}.deltaCRUSH.summary.svg",
        showfliers=args.showfliers,
        alpha=args.alpha_points,
    )

if __name__ == "__main__":
    main()
