#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
import warnings
warnings.filterwarnings('ignore')


COMP_DIR = "compBigWigs/"

COMPARTMENTS = {
    "HCT116": f"{COMP_DIR}HCT116_compartment.bw",
    "HepG2":  f"{COMP_DIR}HepG2_compartment.bw",
    "IMR90":  f"{COMP_DIR}IMR90_compartment.bw",
    "K562":   f"{COMP_DIR}K562_compartment.bw",
}

BRUSEQ = {
    "HCT116": "HCT116_bruseq.bw",
    "HepG2":  "hepg2_bruseq.bw",
    "IMR90":  "imr90_bruseq.bw",
    "K562":   "k562_bruseq.bw",
}

HISTONES = {
    "HCT116": {
        "H3K4me1":  "ENCFF205ZML_k4me1.bigWig",
        "H3K27me3": "ENCFF977LFP_k27me3.bigWig",
        "H3K27ac":  "ENCFF445BLD_k27ac.bigWig",
    },
    "HepG2": {
        "H3K4me1":  "ENCFF436IHQ_k4me1.bigWig",
        "H3K27me3": "ENCFF695IBM_k27me3.bigWig",
        "H3K27ac":  "ENCFF279IEE_hg19_k27ac.bigWig",
    },
    "IMR90": {
        "H3K4me1":  "ENCFF159UXL_k4me1.bigWig",
        "H3K27me3": "ENCFF328UMQ_k27me3.bigWig",
        "H3K27ac":  "ENCFF801HPN_k27ac.bigWig",
    },
    "K562": {
        "H3K4me1":  "ENCFF900IMQ_k4me1.bigWig",
        "H3K27me3": "ENCFF552AKP_k27me3.bigWig",
        "H3K27ac":  "ENCFF840LLW_k27ac.bigWig",
    },
}


BIN_SIZE   = 1000
CHROMS_CHR = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CHROMS_NUM = [str(i) for i in range(1, 23)] + ["X"]
N_CV_FOLDS = 5
OUTDIR     = "compartment_prediction_results_nok9_nohela_1kb"
os.makedirs(OUTDIR, exist_ok=True)

FEATURES = ["BRU-seq", "H3K4me1", "H3K27me3", "H3K27ac"]


def detect_chroms(bw_path):
    bw = pyBigWig.open(bw_path)
    keys = set(bw.chroms().keys())
    bw.close()
    if "chr1" in keys:
        print(f"    [chrom style: UCSC chr-prefix]")
        return CHROMS_CHR
    elif "1" in keys:
        print(f"    [chrom style: Ensembl numeric]")
        return CHROMS_NUM
    else:
        raise ValueError(f"Cannot detect chrom naming in {bw_path}. Keys: {list(keys)[:5]}")


def normalize_chrom(chrom, target_style):
    if target_style == "chr":
        return chrom if chrom.startswith("chr") else f"chr{chrom}"
    else:
        return chrom.replace("chr", "") if chrom.startswith("chr") else chrom


def score_bw_bins(bw_path, chroms, bin_size):
    records = []
    try:
        bw = pyBigWig.open(bw_path)
    except Exception as e:
        raise FileNotFoundError(f"Cannot open {bw_path}: {e}")

    bw_keys = set(bw.chroms().keys())
    if "chr1" in bw_keys:
        bw_chroms = {c: normalize_chrom(c, "chr") for c in chroms}
    elif "1" in bw_keys:
        bw_chroms = {c: normalize_chrom(c, "num") for c in chroms}
    else:
        bw.close()
        return records

    chrom_sizes = bw.chroms()
    for our_chrom, bw_chrom in bw_chroms.items():
        if bw_chrom not in chrom_sizes:
            continue
        clen = chrom_sizes[bw_chrom]
        for s in range(0, clen - bin_size, bin_size):
            e = min(s + bin_size, clen)
            try:
                val = bw.stats(bw_chrom, s, e, type="mean")[0]
            except Exception:
                val = None
            ucsc_chrom = our_chrom if our_chrom.startswith("chr") else f"chr{our_chrom}"
            records.append((ucsc_chrom, s, e, val if val is not None else np.nan))
    bw.close()
    return records


def build_matrix(cell_line):
    print(f"\n{'='*60}")
    print(f"  Building matrix for {cell_line}")
    print(f"{'='*60}")

    chroms = detect_chroms(COMPARTMENTS[cell_line])

    print("  Scoring compartments ...", flush=True)
    comp_recs = score_bw_bins(COMPARTMENTS[cell_line], chroms, BIN_SIZE)
    df = pd.DataFrame(comp_recs, columns=["chrom", "start", "end", "PC1"])
    print(f"    → {len(df)} raw bins, non-null: {df['PC1'].notna().sum()}")

    print("  Scoring BRU-seq ...", flush=True)
    bru_recs = score_bw_bins(BRUSEQ[cell_line], chroms, BIN_SIZE)
    tmp = pd.DataFrame(bru_recs, columns=["chrom", "start", "end", "BRU-seq"])
    df  = df.merge(tmp, on=["chrom", "start", "end"], how="left")

    for feat in [f for f in FEATURES if f != "BRU-seq"]:
        bw_path = HISTONES[cell_line].get(feat)
        if not bw_path or not os.path.exists(bw_path):
            print(f"  WARNING: Missing {feat} for {cell_line}, filling NaN")
            df[feat] = np.nan
            continue
        print(f"  Scoring {feat} ...", flush=True)
        recs = score_bw_bins(bw_path, chroms, BIN_SIZE)
        tmp  = pd.DataFrame(recs, columns=["chrom", "start", "end", feat])
        df   = df.merge(tmp, on=["chrom", "start", "end"], how="left")

    df = df[df["PC1"].notna() & (df["PC1"] != 0)].copy()

    active_cols = ["PC1"] + FEATURES
    df = df.dropna(subset=[c for c in active_cols if c in df.columns]).copy()

    df["compartment"] = (df["PC1"] > 0).astype(int)

    n_A = df["compartment"].sum()
    n_B = (df["compartment"] == 0).sum()
    print(f"  A compartment bins: {n_A:,}  |  B compartment bins: {n_B:,}")

    return df


def run_models(df, cell_line):
    y = df["compartment"].values
    results = []

    pipe = lambda: Pipeline([
        ("scaler", StandardScaler()),
        ("lr", LogisticRegression(max_iter=1000, solver="lbfgs",
                                   class_weight="balanced")),
    ])
    cv = StratifiedKFold(n_splits=N_CV_FOLDS, shuffle=True, random_state=42)

    for feat in FEATURES:
        if feat not in df.columns or df[feat].isna().all():
            continue
        X    = df[[feat]].values
        aucs = cross_val_score(pipe(), X, y, cv=cv, scoring="roc_auc")
        print(f"  {feat:12s}  AUC = {aucs.mean():.3f} ± {aucs.std():.3f}")
        results.append({"cell_line": cell_line, "model": feat,
                        "auc_mean": aucs.mean(), "auc_std": aucs.std()})

    available = [f for f in FEATURES if f in df.columns and not df[f].isna().all()]
    X_all    = df[available].values
    aucs_all = cross_val_score(pipe(), X_all, y, cv=cv, scoring="roc_auc")
    print(f"  {'ALL (multi)':12s}  AUC = {aucs_all.mean():.3f} ± {aucs_all.std():.3f}")
    results.append({"cell_line": cell_line, "model": "All (multi)",
                    "auc_mean": aucs_all.mean(), "auc_std": aucs_all.std()})

    return results


all_results = []

for cl in COMPARTMENTS.keys():
    try:
        df  = build_matrix(cl)
        res = run_models(df, cl)
        all_results.extend(res)
    except FileNotFoundError as e:
        print(f"  SKIPPING {cl}: {e}")

summary = pd.DataFrame(all_results)
summary.to_csv(f"{OUTDIR}/compartment_prediction_summary.csv", index=False)
print(f"\nSummary saved → {OUTDIR}/compartment_prediction_summary.csv")
print(summary.to_string(index=False))


cell_lines = summary["cell_line"].unique()
models     = summary["model"].unique()

fig, ax = plt.subplots(figsize=(12, 5))
n_models = len(models)
x        = np.arange(len(cell_lines))
width    = 0.75 / n_models
colors   = sns.color_palette("Set2", n_models)

for i, model in enumerate(models):
    sub    = summary[summary["model"] == model].set_index("cell_line")
    means  = [sub.loc[cl, "auc_mean"] if cl in sub.index else np.nan for cl in cell_lines]
    stds   = [sub.loc[cl, "auc_std"]  if cl in sub.index else 0      for cl in cell_lines]
    offset = (i - n_models / 2 + 0.5) * width
    ax.bar(x + offset, means, width=width * 0.9, yerr=stds,
           capsize=2, label=model, color=colors[i])

ax.axhline(0.5, color="gray", linestyle="--", lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels(cell_lines)
ax.set_ylabel("CV AUC-ROC")
ax.set_ylim(0, 1)
ax.set_title("Predictiveness of Each Feature for A Compartment Identity\n(univariate logistic regression per feature)")
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
plt.tight_layout()
fig.savefig(f"{OUTDIR}/grouped_bar_auc.pdf", dpi=150)
fig.savefig(f"{OUTDIR}/grouped_bar_auc.png", dpi=150)
plt.close(fig)

pivot = summary[summary["model"] != "All (multi)"].pivot(
    index="model", columns="cell_line", values="auc_mean")
pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(7, 4))
sns.heatmap(pivot, annot=True, fmt=".3f", cmap="YlOrRd",
            vmin=0.5, vmax=1.0, linewidths=0.5, ax=ax)
ax.set_title("AUC-ROC: Predicting A Compartment Identity\n(rows sorted by mean AUC, univariate)")
ax.set_ylabel("")
plt.tight_layout()
fig.savefig(f"{OUTDIR}/auc_heatmap.pdf", dpi=150)
fig.savefig(f"{OUTDIR}/auc_heatmap.png", dpi=150)
plt.close(fig)

rank = summary[summary["model"] != "All (multi)"].groupby("model")["auc_mean"].mean().sort_values(ascending=False)

fig, ax = plt.subplots(figsize=(6, 4))
colors_rank = sns.color_palette("Set2", len(rank))
ax.barh(rank.index[::-1], rank.values[::-1], color=colors_rank[::-1])
ax.axvline(0.5, color="gray", linestyle="--", lw=0.8)
ax.set_xlabel("Mean CV AUC-ROC (across cell lines)")
ax.set_title("Feature Ranking: A Compartment Predictiveness")
ax.set_xlim(0, 1)
plt.tight_layout()
fig.savefig(f"{OUTDIR}/feature_ranking.pdf", dpi=150)
fig.savefig(f"{OUTDIR}/feature_ranking.png", dpi=150)
plt.close(fig)

print(f"\nAll plots saved to {OUTDIR}/")
print("\nFeature ranking (mean AUC across cell lines):")
print(rank.to_string())
print("\nDone!")
