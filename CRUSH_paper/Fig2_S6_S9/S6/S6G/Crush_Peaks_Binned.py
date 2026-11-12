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
    """Check whether the first data line of a file uses 'chr' prefix."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            chrom = line.split()[0]
            return chrom.startswith("chr")
    return False


def harmonize_bed_file(input_file, target_has_chr):
    """
    Read a BED/peak/bedGraph file, harmonize chr naming to match target,
    strip track/comment lines, write to a new temp file. Returns path.
    """
    source_has_chr = detect_chr_prefix(input_file)

    if source_has_chr == target_has_chr:
        out = tempfile.NamedTemporaryFile(delete=False, suffix=".bed").name
        cmd = f"grep -v '^track\\|^#' {input_file} > {out}"
        subprocess.run(cmd, shell=True, check=True)
        return out

    out = tempfile.NamedTemporaryFile(delete=False, suffix=".bed").name
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
    """Sort and merge a (already harmonized) peak/BED file. Returns path."""
    merged = tempfile.NamedTemporaryFile(delete=False, suffix=".bed").name
    cmd = f"sort -k1,1 -k2,2n {peak_file} | bedtools merge > {merged}"
    subprocess.run(cmd, shell=True, check=True)
    return merged


def bedtools_coverage(bins_bed, peaks_bed):
    """
    Run bedtools coverage -a bins -b peaks.
    Returns DataFrame: chr, start, end, n_features, n_bases, bin_len, frac.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv").name
    cmd = f"bedtools coverage -a {bins_bed} -b {peaks_bed} > {tmp}"
    subprocess.run(cmd, shell=True, check=True)

    df = pd.read_csv(tmp, sep="\t", header=None)
    os.remove(tmp)
    df.columns = ["chr", "start", "end", "n_features", "n_bases", "bin_len", "frac"]
    return df


def get_crush_per_bin(bins_bed, crush_bed_harmonized):
    """
    Intersect bins with CRUSH bedGraph.
    Returns Series indexed by (chr, start, end) with CRUSH score.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv").name
    cmd = f"bedtools intersect -a {bins_bed} -b {crush_bed_harmonized} -wa -wb > {tmp}"
    subprocess.run(cmd, shell=True, check=True)

    if os.path.getsize(tmp) == 0:
        os.remove(tmp)
        return pd.Series(dtype=float)

    df = pd.read_csv(tmp, sep="\t", header=None)
    os.remove(tmp)
    df.columns = ["chr", "start", "end", "c_chr", "c_start", "c_end", "crush_score"]
    grouped = df.groupby(["chr", "start", "end"])["crush_score"].mean()
    return grouped


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
        description="Binned 1kb combinatorial CRUSH analysis with threshold sweep."
    )
    ap.add_argument("--bins", required=True,
                    help="1kb bin BED file (chr start end) — used as reference for chr naming")
    ap.add_argument("--crush", required=True,
                    help="CRUSH bedGraph (chr start end score)")
    ap.add_argument("--k27", required=True,
                    help="K27me3 peak file")
    ap.add_argument("--active", nargs="+", required=True,
                    help="Active mark peak files as name:path pairs "
                         "(e.g. K27ac:k27ac.bed K36me3:k36me3.bed K4me1:k4me1.bed K4me3:k4me3.bed)")
    ap.add_argument("--k27_thresh", nargs="+", type=float, default=[0.01],
                    help="K27me3 coverage fraction threshold(s)")
    ap.add_argument("--active_thresh", nargs="+", type=float, default=[0.01],
                    help="Active mark coverage fraction threshold(s)")
    ap.add_argument("--label", required=True, help="Cell line label")
    ap.add_argument("--outprefix", required=True, help="Output prefix")

    args = ap.parse_args()

    active_names = []
    active_files = []
    for item in args.active:
        if ":" not in item:
            raise ValueError(f"Expected name:path for --active but got: {item}\n"
                             f"Format: K27ac:k27ac.bed K36me3:k36me3.bed ...")
        name, path = item.split(":", 1)
        active_names.append(name)
        active_files.append(path)

    n_active = len(active_files)

    print(f"Cell line: {args.label}")
    print(f"Active marks ({n_active}): {', '.join(active_names)}")
    print(f"K27 thresholds: {args.k27_thresh}")
    print(f"Active thresholds: {args.active_thresh}")

    bins_has_chr = detect_chr_prefix(args.bins)
    print(f"Bins file uses {'chr' if bins_has_chr else 'no-chr'} prefix — harmonizing all inputs to match.")

    temp_files = []

    print("Harmonizing CRUSH bedGraph...")
    crush_harmonized = harmonize_bed_file(args.crush, bins_has_chr)
    temp_files.append(crush_harmonized)

    print("Harmonizing K27me3 peaks...")
    k27_harmonized = harmonize_bed_file(args.k27, bins_has_chr)
    temp_files.append(k27_harmonized)

    active_harmonized = []
    for i, peak_file in enumerate(active_files):
        print(f"Harmonizing {active_names[i]}...")
        h = harmonize_bed_file(peak_file, bins_has_chr)
        active_harmonized.append(h)
        temp_files.append(h)

    print("\nComputing K27me3 coverage per bin...")
    k27_merged = sort_and_merge(k27_harmonized)
    temp_files.append(k27_merged)
    k27_cov = bedtools_coverage(args.bins, k27_merged)

    print("Computing active mark coverage per bin...")
    active_covs = []
    for i, h_file in enumerate(active_harmonized):
        print(f"  {active_names[i]}...")
        merged = sort_and_merge(h_file)
        temp_files.append(merged)
        cov = bedtools_coverage(args.bins, merged)
        active_covs.append(cov)

    print("Getting CRUSH scores per bin...")
    crush_scores = get_crush_per_bin(args.bins, crush_harmonized)

    bins_df = k27_cov[["chr", "start", "end"]].copy()
    bins_df["k27_frac"] = k27_cov["frac"].values

    for i, name in enumerate(active_names):
        bins_df[f"{name}_frac"] = active_covs[i]["frac"].values

    bins_df = bins_df.set_index(["chr", "start", "end"])
    bins_df["crush"] = crush_scores
    bins_df = bins_df.dropna(subset=["crush"]).reset_index()

    print(f"Total bins with CRUSH scores: {len(bins_df):,}")

    for tf in temp_files:
        try:
            os.remove(tf)
        except Exception:
            pass

    for k27_t in args.k27_thresh:
        for act_t in args.active_thresh:
            tag = f"k27t{k27_t}_actt{act_t}"
            print(f"\n--- Threshold: K27≥{k27_t}, Active≥{act_t} ---")

            df = bins_df.copy()
            df["k27_bin"] = (df["k27_frac"] >= k27_t).astype(int)

            active_bin_cols = []
            for name in active_names:
                col = f"{name}_bin"
                df[col] = (df[f"{name}_frac"] >= act_t).astype(int)
                active_bin_cols.append(col)

            df["active_sum"] = df[active_bin_cols].sum(axis=1)

            for n in range(n_active + 1):
                mask = (df["k27_bin"] == 1) & (df["active_sum"] == n)
                if n == 0:
                    lbl = "K27me3 only"
                else:
                    lbl = f"K27me3 + {n} active"
                df.loc[mask, "group"] = lbl

            mask_active_only = (df["k27_bin"] == 0) & (df["active_sum"] >= 1)
            df.loc[mask_active_only, "group"] = "Active only"

            df = df.dropna(subset=["group"])

            group_order = ["K27me3 only"]
            for n in range(1, n_active + 1):
                group_order.append(f"K27me3 + {n} active")
            group_order.append("Active only")
            group_order = [g for g in group_order if g in df["group"].unique()]

            matrix_out = f"{args.outprefix}_{args.label}_{tag}_bin_matrix.tsv.gz"
            save_cols = ["chr", "start", "end", "k27_bin"] + active_bin_cols + ["active_sum", "crush", "group"]
            df[save_cols].to_csv(matrix_out, sep="\t", index=False, compression="gzip")
            print(f"  Bin matrix     -> {matrix_out}")

            scored = df[["group", "crush"]].copy()
            scored.columns = ["Group", "Score"]
            scored["CellLine"] = args.label
            scored["K27_thresh"] = k27_t
            scored["Active_thresh"] = act_t

            scores_out = f"{args.outprefix}_{args.label}_{tag}_binned_scores.tsv"
            scored.to_csv(scores_out, sep="\t", index=False)
            print(f"  Scores TSV     -> {scores_out}")

            summary_rows = []
            for grp in group_order:
                vals = df.loc[df["group"] == grp, "crush"]
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
            summary_df = pd.DataFrame(summary_rows)
            summary_out = f"{args.outprefix}_{args.label}_{tag}_distribution_summary.tsv"
            summary_df.to_csv(summary_out, sep="\t", index=False)
            print(f"  Summary        -> {summary_out}")

            k27_bins = df[df["k27_bin"] == 1]
            total_k27 = len(k27_bins)
            overlap_rows = []
            for n in range(n_active + 1):
                cnt = (k27_bins["active_sum"] == n).sum()
                pct = 100.0 * cnt / total_k27 if total_k27 > 0 else 0
                overlap_rows.append({"Active_count": n, "N": cnt, "Percent": round(pct, 2)})

            for ge in range(1, n_active + 1):
                cnt = (k27_bins["active_sum"] >= ge).sum()
                pct = 100.0 * cnt / total_k27 if total_k27 > 0 else 0
                overlap_rows.append({"Active_count": f"≥{ge}", "N": cnt, "Percent": round(pct, 2)})

            overlap_df = pd.DataFrame(overlap_rows)
            overlap_out = f"{args.outprefix}_{args.label}_{tag}_overlap_counts.tsv"
            overlap_df.to_csv(overlap_out, sep="\t", index=False)
            print(f"  Overlap counts -> {overlap_out}")

            n_active_overlap_k27 = len(df[(df["k27_bin"] == 1) & (df["active_sum"] >= 1)])
            n_k27_only_bins = len(df[(df["k27_bin"] == 1) & (df["active_sum"] == 0)])
            n_active_with_k27 = n_active_overlap_k27
            n_active_only_bins = len(df[(df["k27_bin"] == 0) & (df["active_sum"] >= 1)])
            total_active_bins = n_active_with_k27 + n_active_only_bins

            fig_ov, (ax_k, ax_a) = plt.subplots(1, 2, figsize=(10, 4.5))

            k27_cats = [f"{n} active" for n in range(n_active + 1)]
            k27_cnts = [(k27_bins["active_sum"] == n).sum() for n in range(n_active + 1)]
            k27_pcts = [100.0 * c / total_k27 if total_k27 > 0 else 0 for c in k27_cnts]
            k27_colors = list(reversed(sns.color_palette("RdYlGn", n_colors=max(n_active + 1, 3))))

            bars_k = ax_k.bar(k27_cats, k27_pcts, color=k27_colors, edgecolor="black", linewidth=0.6)
            for bar, pct, cnt in zip(bars_k, k27_pcts, k27_cnts):
                ax_k.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                         f"{pct:.1f}%\n(n={cnt:,})", ha="center", va="bottom", fontsize=7)
            ax_k.set_ylabel("% of K27me3 bins")
            ax_k.set_title("K27me3 perspective", fontsize=11, fontweight="bold")
            ax_k.set_ylim(0, max(k27_pcts + [1]) + 15)
            ax_k.tick_params(axis="x", rotation=15, labelsize=9)

            if total_active_bins > 0:
                act_pct_k27 = 100.0 * n_active_with_k27 / total_active_bins
                act_pct_only = 100.0 * n_active_only_bins / total_active_bins
            else:
                act_pct_k27 = act_pct_only = 0

            bars_a = ax_a.bar(["∩ K27me3", "Active\nonly"],
                              [act_pct_k27, act_pct_only],
                              color=["#66c2a5", "#fc8d62"], edgecolor="black", linewidth=0.6)
            for bar, pct, cnt in zip(bars_a, [act_pct_k27, act_pct_only],
                                      [n_active_with_k27, n_active_only_bins]):
                ax_a.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                         f"{pct:.1f}%\n(n={cnt:,})", ha="center", va="bottom", fontsize=9)
            ax_a.set_ylabel("% of Active bins")
            ax_a.set_title("Active-mark perspective", fontsize=11, fontweight="bold")
            ax_a.set_ylim(0, max(act_pct_k27, act_pct_only, 1) + 15)

            fig_ov.suptitle(f"{args.label}: Bin overlap — K27me3 & Active marks\n"
                            f"(K27≥{k27_t}, Active≥{act_t})", fontsize=12, y=1.03)
            plt.tight_layout()
            ov_svg = f"{args.outprefix}_{args.label}_{tag}_overlap_pct.svg"
            fig_ov.savefig(ov_svg, format="svg", bbox_inches="tight")
            plt.close(fig_ov)
            print(f"  Overlap plot   -> {ov_svg}")

            ov_pct_rows = [
                {"Perspective": "K27me3", "Category": f"{n} active",
                 "Count": (k27_bins["active_sum"] == n).sum(),
                 "Total": total_k27,
                 "Percent": round(100.0 * (k27_bins["active_sum"] == n).sum() / total_k27, 2) if total_k27 else 0}
                for n in range(n_active + 1)
            ] + [
                {"Perspective": "Active", "Category": "∩ K27me3",
                 "Count": n_active_with_k27, "Total": total_active_bins,
                 "Percent": round(act_pct_k27, 2)},
                {"Perspective": "Active", "Category": "Active only",
                 "Count": n_active_only_bins, "Total": total_active_bins,
                 "Percent": round(act_pct_only, 2)},
            ]
            ov_pct_out = f"{args.outprefix}_{args.label}_{tag}_overlap_pct.tsv"
            pd.DataFrame(ov_pct_rows).to_csv(ov_pct_out, sep="\t", index=False)
            print(f"  Overlap % TSV  -> {ov_pct_out}")

            groups_dict = {}
            for grp in group_order:
                vals = df.loc[df["group"] == grp, "crush"].values
                if len(vals) >= 2:
                    groups_dict[grp] = vals

            if len(groups_dict) >= 2:
                stats_df = run_pairwise_tests(groups_dict)
                stats_out = f"{args.outprefix}_{args.label}_{tag}_pairwise_stats.tsv"
                stats_df.to_csv(stats_out, sep="\t", index=False)
                print(f"  Pairwise stats -> {stats_out}")

            scored_plot = remove_outliers_iqr_per_group(scored, col="Score", group_col="Group")

            n_groups = len(group_order)
            palette_colors = list(reversed(sns.color_palette("RdYlGn", n_colors=max(n_groups, 3))))
            palette = {g: palette_colors[i] for i, g in enumerate(group_order)}

            plt.figure(figsize=(max(8, n_groups * 1.2), 5))
            sns.boxplot(data=scored_plot, x="Group", y="Score",
                        order=group_order, palette=palette,
                        width=0.6, linewidth=1.2, fliersize=0,
                        medianprops=dict(color="black", linewidth=1.8),
                        boxprops=dict(edgecolor="black"),
                        whiskerprops=dict(color="black"),
                        capprops=dict(color="black"))

            plt.axhline(0, linestyle="--", color="gray", alpha=0.6)
            plt.title(f"{args.label}: CRUSH by K27me3 + active mark count\n"
                       f"(K27≥{k27_t}, Active≥{act_t})")
            plt.ylabel("CRUSH Compartment Score")
            plt.xlabel("")
            plt.xticks(rotation=20, ha="right")
            plt.tight_layout()

            box_svg = f"{args.outprefix}_{args.label}_{tag}_boxplot.svg"
            plt.savefig(box_svg, format="svg")
            plt.close()
            print(f"  Boxplot        -> {box_svg}")

    print("\nDone.")


if __name__ == "__main__":
    main()
