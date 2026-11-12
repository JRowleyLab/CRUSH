import argparse
import os
import subprocess
import tempfile

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from itertools import combinations
import numpy as np


def detect_chr_prefix(filepath):
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            return line.split()[0].startswith("chr")
    return False


def harmonize_bed_file(input_file, target_has_chr):
    source_has_chr = detect_chr_prefix(input_file)
    out = tempfile.NamedTemporaryFile(delete=False, suffix=".bed").name
    if source_has_chr == target_has_chr:
        subprocess.run(f"grep -v '^track\\|^#' {input_file} > {out}", shell=True, check=True)
        return out
    with open(input_file) as fin, open(out, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom = fields[0]
            if target_has_chr and not chrom.startswith("chr"):
                fields[0] = "chr" + chrom
            elif not target_has_chr and chrom.startswith("chr"):
                fields[0] = chrom.replace("chr", "", 1)
            fout.write("\t".join(fields) + "\n")
    return out


def sort_and_merge(peak_file):
    merged = tempfile.NamedTemporaryFile(delete=False, suffix=".bed").name
    subprocess.run(f"sort -k1,1 -k2,2n {peak_file} | bedtools merge > {merged}",
                   shell=True, check=True)
    return merged


def bedtools_coverage(bins_bed, peaks_bed):
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv").name
    subprocess.run(f"bedtools coverage -a {bins_bed} -b {peaks_bed} > {tmp}",
                   shell=True, check=True)
    df = pd.read_csv(tmp, sep="\t", header=None)
    os.remove(tmp)
    df.columns = ["chr", "start", "end", "n_features", "n_bases", "bin_len", "frac"]
    return df


def get_crush_per_bin(bins_bed, crush_harmonized):
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv").name
    subprocess.run(f"bedtools intersect -a {bins_bed} -b {crush_harmonized} -wa -wb > {tmp}",
                   shell=True, check=True)
    if os.path.getsize(tmp) == 0:
        os.remove(tmp)
        return pd.Series(dtype=float)
    df = pd.read_csv(tmp, sep="\t", header=None)
    os.remove(tmp)
    df.columns = ["chr", "start", "end", "c_chr", "c_start", "c_end", "crush_score"]
    return df.groupby(["chr", "start", "end"])["crush_score"].mean()


def significance_stars(pval):
    if pval < 1e-4:   return "****"
    elif pval < 1e-3:  return "***"
    elif pval < 1e-2:  return "**"
    elif pval < 0.05:  return "*"
    return "ns"


def run_pairwise_tests(groups_dict):
    rows = []
    for g1, g2 in combinations(groups_dict.keys(), 2):
        s1, s2 = groups_dict[g1], groups_dict[g2]
        if len(s1) < 2 or len(s2) < 2:
            rows.append({"Group1": g1, "Group2": g2, "N1": len(s1), "N2": len(s2),
                         "Median1": np.nan, "Median2": np.nan, "U_stat": np.nan,
                         "P_value": np.nan, "EffectSize_RBC": np.nan, "Significance": "NA"})
            continue
        stat, pval = mannwhitneyu(s1, s2, alternative="two-sided")
        n1, n2 = len(s1), len(s2)
        rbc = 1 - (2 * stat) / (n1 * n2)
        rows.append({"Group1": g1, "Group2": g2, "N1": n1, "N2": n2,
                      "Median1": round(np.median(s1), 4), "Median2": round(np.median(s2), 4),
                      "U_stat": round(stat, 4), "P_value": pval,
                      "EffectSize_RBC": round(rbc, 4), "Significance": significance_stars(pval)})
    return pd.DataFrame(rows)


def remove_outliers_iqr_per_group(df, col="Score", group_col="Group"):
    cleaned = []
    for grp, sub in df.groupby(group_col):
        q1 = sub[col].quantile(0.25)
        q3 = sub[col].quantile(0.75)
        iqr = q3 - q1
        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr
        cleaned.append(sub[(sub[col] >= lower) & (sub[col] <= upper)])
    return pd.concat(cleaned, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(
        description="Active-only dose-response: CRUSH by number of active marks (no K27me3)."
    )
    ap.add_argument("--bins", required=True, help="1kb bin BED file")
    ap.add_argument("--crush", required=True, help="CRUSH bedGraph")
    ap.add_argument("--k27", required=True, help="K27me3 peak file")
    ap.add_argument("--active", nargs="+", required=True,
                    help="Active marks as name:path (e.g. K27ac:k27ac.bed)")
    ap.add_argument("--k27_thresh", type=float, default=0.2,
                    help="K27me3 threshold (bins BELOW this are 'no K27me3')")
    ap.add_argument("--active_thresh", type=float, default=0.2,
                    help="Active mark threshold")
    ap.add_argument("--label", required=True, help="Cell line label")
    ap.add_argument("--outprefix", required=True, help="Output prefix")

    args = ap.parse_args()

    active_names, active_files = [], []
    for item in args.active:
        if ":" not in item:
            raise ValueError(f"Expected name:path, got: {item}")
        name, path = item.split(":", 1)
        active_names.append(name)
        active_files.append(path)
    n_active = len(active_files)

    print(f"Cell line: {args.label}")
    print(f"Active marks ({n_active}): {', '.join(active_names)}")
    print(f"K27 thresh: {args.k27_thresh}, Active thresh: {args.active_thresh}")

    bins_has_chr = detect_chr_prefix(args.bins)
    print(f"Bins use {'chr' if bins_has_chr else 'no-chr'} prefix.")

    temp_files = []

    crush_h = harmonize_bed_file(args.crush, bins_has_chr)
    k27_h = harmonize_bed_file(args.k27, bins_has_chr)
    temp_files.extend([crush_h, k27_h])

    active_harmonized = []
    for i, f in enumerate(active_files):
        h = harmonize_bed_file(f, bins_has_chr)
        active_harmonized.append(h)
        temp_files.append(h)

    print("Computing K27me3 coverage...")
    k27_merged = sort_and_merge(k27_h)
    temp_files.append(k27_merged)
    k27_cov = bedtools_coverage(args.bins, k27_merged)

    print("Computing active mark coverage...")
    active_covs = []
    for i, h in enumerate(active_harmonized):
        print(f"  {active_names[i]}...")
        merged = sort_and_merge(h)
        temp_files.append(merged)
        cov = bedtools_coverage(args.bins, merged)
        active_covs.append(cov)

    print("Getting CRUSH scores...")
    crush_scores = get_crush_per_bin(args.bins, crush_h)

    bins_df = k27_cov[["chr", "start", "end"]].copy()
    bins_df["k27_frac"] = k27_cov["frac"].values
    for i, name in enumerate(active_names):
        bins_df[f"{name}_frac"] = active_covs[i]["frac"].values

    bins_df = bins_df.set_index(["chr", "start", "end"])
    bins_df["crush"] = crush_scores
    bins_df = bins_df.dropna(subset=["crush"]).reset_index()

    print(f"Total bins with CRUSH: {len(bins_df):,}")

    for tf in temp_files:
        try:
            os.remove(tf)
        except Exception:
            pass

    bins_df["k27_bin"] = (bins_df["k27_frac"] >= args.k27_thresh).astype(int)

    active_bin_cols = []
    for name in active_names:
        col = f"{name}_bin"
        bins_df[col] = (bins_df[f"{name}_frac"] >= args.active_thresh).astype(int)
        active_bin_cols.append(col)

    bins_df["active_sum"] = bins_df[active_bin_cols].sum(axis=1)

    active_only = bins_df[(bins_df["k27_bin"] == 0) & (bins_df["active_sum"] >= 1)].copy()

    print(f"Active-only bins (no K27me3): {len(active_only):,}")

    active_only["group"] = active_only["active_sum"].apply(lambda x: f"{int(x)} active mark{'s' if x > 1 else ''}")

    group_order = [f"{n} active mark{'s' if n > 1 else ''}" for n in range(1, n_active + 1)]
    group_order = [g for g in group_order if g in active_only["group"].unique()]

    scored = active_only[["group", "crush"]].copy()
    scored.columns = ["Group", "Score"]
    scored["CellLine"] = args.label

    scores_out = f"{args.outprefix}_{args.label}_active_doseresponse_scores.tsv"
    scored.to_csv(scores_out, sep="\t", index=False)
    print(f"Scores TSV       -> {scores_out}")

    summary_rows = []
    for grp in group_order:
        vals = active_only.loc[active_only["group"] == grp, "crush"]
        if len(vals) == 0:
            continue
        summary_rows.append({
            "Group": grp, "N": len(vals),
            "Median": round(vals.median(), 4),
            "Mean": round(vals.mean(), 4),
            "Std": round(vals.std(), 4),
            "Q1": round(vals.quantile(0.25), 4),
            "Q3": round(vals.quantile(0.75), 4),
            "Pct_positive": round((vals > 0).mean() * 100, 2),
        })
    pd.DataFrame(summary_rows).to_csv(
        f"{args.outprefix}_{args.label}_active_doseresponse_summary.tsv",
        sep="\t", index=False)
    print(f"Summary          -> {args.outprefix}_{args.label}_active_doseresponse_summary.tsv")

    groups_dict = {}
    for grp in group_order:
        vals = active_only.loc[active_only["group"] == grp, "crush"].values
        if len(vals) >= 2:
            groups_dict[grp] = vals

    if len(groups_dict) >= 2:
        stats_df = run_pairwise_tests(groups_dict)
        stats_out = f"{args.outprefix}_{args.label}_active_doseresponse_pairwise_stats.tsv"
        stats_df.to_csv(stats_out, sep="\t", index=False)
        print(f"Pairwise stats   -> {stats_out}")

    scored_plot = remove_outliers_iqr_per_group(scored)

    palette_colors = sns.color_palette("YlOrRd", n_colors=max(len(group_order), 3))
    palette = {g: palette_colors[i] for i, g in enumerate(group_order)}

    plt.figure(figsize=(max(6, len(group_order) * 1.5), 5))
    sns.boxplot(data=scored_plot, x="Group", y="Score",
                order=group_order, palette=palette,
                width=0.5, linewidth=1.2, fliersize=0,
                medianprops=dict(color="black", linewidth=1.8),
                boxprops=dict(edgecolor="black"),
                whiskerprops=dict(color="black"),
                capprops=dict(color="black"))

    plt.axhline(0, linestyle="--", color="gray", alpha=0.6)
    plt.title(f"{args.label}: CRUSH by active mark count (no K27me3)\n"
              f"(K27<{args.k27_thresh}, Active≥{args.active_thresh})")
    plt.ylabel("CRUSH Compartment Score")
    plt.xlabel("")
    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()

    plt.savefig(f"{args.outprefix}_{args.label}_active_doseresponse_boxplot.svg", format="svg")
    plt.close()
    print(f"Boxplot          -> {args.outprefix}_{args.label}_active_doseresponse_boxplot.svg")

    print("\nDone.")


if __name__ == "__main__":
    main()
