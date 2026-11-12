import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu
from itertools import combinations


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
    ap = argparse.ArgumentParser(description="Combined active-only dose-response plots.")
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--outprefix", required=True)
    ap.add_argument("--remove_outliers", action="store_true")
    args = ap.parse_args()

    frames = []
    for path in args.inputs:
        tmp = pd.read_csv(path, sep="\t")
        if "CellLine" not in tmp.columns:
            cl = path.rsplit("/", 1)[-1].split("_")[0]
            tmp["CellLine"] = cl
        frames.append(tmp)

    df_raw = pd.concat(frames, ignore_index=True)
    for c in ["Group", "Score", "CellLine"]:
        if c not in df_raw.columns:
            raise ValueError(f"Missing column '{c}'")

    if args.remove_outliers:
        df_plot = remove_outliers_iqr_per_group(df_raw)
    else:
        df_plot = df_raw.copy()

    all_groups = df_raw["Group"].unique()
    order = sorted([g for g in all_groups],
                   key=lambda x: int(x.split()[0]))

    cell_lines = sorted(df_raw["CellLine"].unique())

    palette_colors = sns.color_palette("YlOrRd", n_colors=max(len(order), 3))
    palette = {g: palette_colors[i] for i, g in enumerate(order)}

    cl_colors = sns.color_palette("tab10", n_colors=max(len(cell_lines), 3))
    cl_palette = {cl: cl_colors[i] for i, cl in enumerate(cell_lines)}

    medians = df_raw.groupby(["Group", "CellLine"])["Score"].median().reset_index()
    medians.columns = ["Group", "CellLine", "Median"]

    fig1, ax1 = plt.subplots(figsize=(max(8, len(order) * 1.5), 5.5))

    sns.boxplot(data=df_plot, x="Group", y="Score", order=order,
                palette=palette, width=0.5, fliersize=0,
                linewidth=1.2,
                medianprops=dict(color="black", linewidth=2.0),
                boxprops=dict(edgecolor="black"),
                whiskerprops=dict(color="black"),
                capprops=dict(color="black"),
                ax=ax1, zorder=1)

    np.random.seed(42)
    group_to_x = {g: i for i, g in enumerate(order)}
    for _, row in medians.iterrows():
        grp, cl, med_val = row["Group"], row["CellLine"], row["Median"]
        if grp not in group_to_x:
            continue
        x = group_to_x[grp] + np.random.uniform(-0.15, 0.15)
        ax1.scatter(x, med_val, color=cl_palette[cl], s=60,
                    edgecolor="black", linewidth=0.6, zorder=3,
                    label=cl if grp == order[0] else None)

    ax1.axhline(0, linestyle=":", color="gray", alpha=0.5)
    ax1.set_ylabel("CRUSH Compartment Score", fontsize=12)
    ax1.set_xlabel("")
    ax1.set_title("Combined CRUSH by Active Mark Count (no K27me3)", fontsize=13)
    ax1.tick_params(axis="x", rotation=15, labelsize=10)

    handles, labels = ax1.get_legend_handles_labels()
    seen = {}
    uh, ul = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen[l] = True
            uh.append(h)
            ul.append(l)
    ax1.legend(uh, ul, title="Cell Line (median)", bbox_to_anchor=(1.02, 1),
               loc="upper left", frameon=False, fontsize=9, title_fontsize=10)

    plt.tight_layout()
    fig1.savefig(f"{args.outprefix}_combined_boxplot.svg", format="svg", bbox_inches="tight")
    fig1.savefig(f"{args.outprefix}_combined_boxplot.png", format="png", dpi=300, bbox_inches="tight")
    plt.close(fig1)
    print(f"Combined boxplot -> {args.outprefix}_combined_boxplot.svg")

    n_cl = len(cell_lines)
    ncols = min(4, n_cl)
    nrows = int(np.ceil(n_cl / ncols))

    fig2, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.5 * nrows),
                               sharey=True, squeeze=False)

    for idx, cl in enumerate(cell_lines):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        cl_raw = df_raw[df_raw["CellLine"] == cl]
        cl_plot = remove_outliers_iqr_per_group(cl_raw) if args.remove_outliers else cl_raw
        cl_order = [g for g in order if g in cl_plot["Group"].unique()]

        sns.boxplot(data=cl_plot, x="Group", y="Score", order=cl_order,
                    palette=palette, width=0.5,
                    fliersize=1.5, flierprops=dict(marker=".", alpha=0.3),
                    linewidth=1.0,
                    medianprops=dict(color="black", linewidth=1.5),
                    boxprops=dict(edgecolor="black"),
                    whiskerprops=dict(color="black"),
                    capprops=dict(color="black"),
                    ax=ax)

        ax.axhline(0, linestyle=":", color="gray", alpha=0.4)
        ax.set_title(cl, fontsize=11, fontweight="bold")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=25, labelsize=8)
        ax.set_ylabel("CRUSH Score" if col == 0 else "", fontsize=10)

    for idx in range(n_cl, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    fig2.suptitle("CRUSH by Active Mark Count (no K27me3) — per Cell Line", fontsize=14, y=1.01)
    plt.tight_layout()
    fig2.savefig(f"{args.outprefix}_faceted_boxplot.svg", format="svg", bbox_inches="tight")
    fig2.savefig(f"{args.outprefix}_faceted_boxplot.png", format="png", dpi=300, bbox_inches="tight")
    plt.close(fig2)
    print(f"Faceted boxplot  -> {args.outprefix}_faceted_boxplot.svg")

    pooled = df_raw.groupby("Group")["Score"].agg(
        N="count", Median="median", Mean="mean", Std="std",
        Q1=lambda x: x.quantile(0.25), Q3=lambda x: x.quantile(0.75),
        Pct_positive=lambda x: round((x > 0).mean() * 100, 2)
    ).reindex(order)
    pooled.to_csv(f"{args.outprefix}_pooled_stats.tsv", sep="\t")

    detailed = df_raw.groupby(["Group", "CellLine"])["Score"].agg(
        N="count", Median="median", Mean="mean", Std="std",
        Q1=lambda x: x.quantile(0.25), Q3=lambda x: x.quantile(0.75)
    )
    detailed.to_csv(f"{args.outprefix}_cellline_stats.tsv", sep="\t")

    groups_dict = {}
    for grp in order:
        vals = df_raw.loc[df_raw["Group"] == grp, "Score"].values
        if len(vals) >= 2:
            groups_dict[grp] = vals

    if len(groups_dict) >= 2:
        rows = []
        for g1, g2 in combinations(groups_dict.keys(), 2):
            s1, s2 = groups_dict[g1], groups_dict[g2]
            stat, pval = mannwhitneyu(s1, s2, alternative="two-sided")
            n1, n2 = len(s1), len(s2)
            rbc = 1 - (2 * stat) / (n1 * n2)
            if pval < 1e-4:     stars = "****"
            elif pval < 1e-3:   stars = "***"
            elif pval < 1e-2:   stars = "**"
            elif pval < 0.05:   stars = "*"
            else:               stars = "ns"
            rows.append({
                "Group1": g1, "Group2": g2, "N1": n1, "N2": n2,
                "Median1": round(np.median(s1), 4), "Median2": round(np.median(s2), 4),
                "U_stat": round(stat, 4), "P_value": pval,
                "EffectSize_RBC": round(rbc, 4), "Significance": stars
            })
        pd.DataFrame(rows).to_csv(f"{args.outprefix}_pairwise_stats.tsv", sep="\t", index=False)
        print(f"Pairwise stats   -> {args.outprefix}_pairwise_stats.tsv")

    print(f"Pooled stats     -> {args.outprefix}_pooled_stats.tsv")
    print(f"Cell-line stats  -> {args.outprefix}_cellline_stats.tsv")
    print("\nDone.")


if __name__ == "__main__":
    main()
