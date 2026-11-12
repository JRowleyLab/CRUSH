import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu
from matplotlib.patches import FancyBboxPatch


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


def sig_label(pval):
    if pval < 1e-4:   return "****"
    elif pval < 1e-3: return "***"
    elif pval < 1e-2: return "**"
    elif pval < 0.05: return "*"
    else:             return "ns"


def draw_significance_bracket(ax, x1, x2, y, h, label, color="black", fontsize=9):
    """Draw a bracket with significance label above two boxes."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, color=color)
    ax.text((x1 + x2) / 2, y + h * 1.1, label,
            ha="center", va="bottom", fontsize=fontsize,
            color=color, fontweight="bold")


def draw_pair_box(ax, x_left, x_right, y_min, y_max, color, label=None, alpha=0.08):
    """Draw a rounded rectangle behind a pair of boxes."""
    pad_x = 0.42
    pad_y = (y_max - y_min) * 0.02
    width  = (x_right - x_left) + 2 * pad_x
    height = (y_max - y_min) + 2 * pad_y
    rect = FancyBboxPatch(
        (x_left - pad_x, y_min - pad_y),
        width, height,
        boxstyle="round,pad=0.05",
        linewidth=1.4,
        edgecolor=color,
        facecolor=color,
        alpha=alpha,
        zorder=0,
        clip_on=False
    )
    ax.add_patch(rect)
    rect_border = FancyBboxPatch(
        (x_left - pad_x, y_min - pad_y),
        width, height,
        boxstyle="round,pad=0.05",
        linewidth=1.6,
        edgecolor=color,
        facecolor="none",
        alpha=0.55,
        zorder=0.5,
        clip_on=False
    )
    ax.add_patch(rect_border)


def main():
    ap = argparse.ArgumentParser(
        description="Interleaved K27me3+active vs active-only combined plot (v2)."
    )
    ap.add_argument("--k27_inputs",    nargs="+", required=True)
    ap.add_argument("--active_inputs", nargs="+", required=True)
    ap.add_argument("--outprefix",     required=True)
    ap.add_argument("--remove_outliers", action="store_true")
    args = ap.parse_args()

    k27_df = pd.concat(
        [pd.read_csv(p, sep="\t") for p in args.k27_inputs],
        ignore_index=True
    )
    k27_df["Source"] = "K27me3"

    act_df = pd.concat(
        [pd.read_csv(p, sep="\t") for p in args.active_inputs],
        ignore_index=True
    )
    act_df["Source"] = "Active"

    rows = []

    for grp_k, label_k in [("K27me3 only", "K27me3\nonly")]:
        sub = k27_df.loc[k27_df["Group"] == grp_k, ["Score", "CellLine"]].copy()
        sub["PlotGroup"] = label_k
        rows.append(sub)

    for grp_a, label_a in [("Active only", "Active\nonly")]:
        sub = k27_df.loc[k27_df["Group"] == grp_a, ["Score", "CellLine"]].copy()
        sub["PlotGroup"] = label_a
        rows.append(sub)

    for n in range(1, 5):
        k27_grp = f"K27me3 + {n} active"
        sub = k27_df.loc[k27_df["Group"] == k27_grp, ["Score", "CellLine"]].copy()
        sub["PlotGroup"] = f"K27me3\n+{n} active"
        rows.append(sub)

        act_grp = f"{n} active mark" if n == 1 else f"{n} active marks"
        sub = act_df.loc[act_df["Group"] == act_grp, ["Score", "CellLine"]].copy()
        sub["PlotGroup"] = f"{n} active\n(no K27)"
        rows.append(sub)

    df = pd.concat(rows, ignore_index=True)

    pair_definitions = [
        ("K27me3\nonly",     "Active\nonly"),
    ] + [
        (f"K27me3\n+{n} active", f"{n} active\n(no K27)")
        for n in range(1, 5)
    ]

    pair_colors = [
        "#6b5b95",
        "#388e3c",
        "#1976d2",
        "#f57c00",
        "#c62828",
    ]

    gap = 1.2
    box_spacing = 1.0

    x_positions = {}
    pair_x_ranges = []
    cur = 0.0
    order = []
    for i, (g1, g2) in enumerate(pair_definitions):
        x1 = cur
        x2 = cur + box_spacing
        x_positions[g1] = x1
        x_positions[g2] = x2
        pair_x_ranges.append((x1, x2))
        order.extend([g1, g2])
        cur = x2 + gap

    order = [g for g in order if g in df["PlotGroup"].unique()]

    palette = {}
    for g in order:
        if "K27me3" in g:
            palette[g] = "#8da0cb"
        else:
            palette[g] = "#fc8d62"

    df_plot = (remove_outliers_iqr_per_group(df, col="Score", group_col="PlotGroup")
               if args.remove_outliers else df.copy())

    cell_lines = sorted(df["CellLine"].unique())
    cl_colors  = sns.color_palette("tab10", n_colors=max(len(cell_lines), 3))
    cl_palette = {cl: cl_colors[i] for i, cl in enumerate(cell_lines)}

    medians = (df.groupby(["PlotGroup", "CellLine"])["Score"]
                 .median().reset_index())
    medians.columns = ["PlotGroup", "CellLine", "Median"]

    stats_rows = []
    pair_stats = {}

    for i, (g1, g2) in enumerate(pair_definitions):
        s1 = df.loc[df["PlotGroup"] == g1, "Score"].values
        s2 = df.loc[df["PlotGroup"] == g2, "Score"].values
        label = f"{g1.replace(chr(10),' ')} vs {g2.replace(chr(10),' ')}"

        if len(s1) < 2 or len(s2) < 2:
            stats_rows.append({
                "Comparison": label,
                "Group1": g1.replace("\n", " "), "Group2": g2.replace("\n", " "),
                "N1": len(s1), "N2": len(s2),
                "Median1": np.nan, "Median2": np.nan,
                "U_stat": np.nan, "P_value": np.nan,
                "EffectSize_RBC": np.nan, "Significance": "NA"
            })
            pair_stats[i] = (np.nan, "NA")
            continue

        stat, pval = mannwhitneyu(s1, s2, alternative="two-sided")
        n1, n2 = len(s1), len(s2)
        rbc = 1 - (2 * stat) / (n1 * n2)
        stars = sig_label(pval)

        stats_rows.append({
            "Comparison": label,
            "Group1": g1.replace("\n", " "), "Group2": g2.replace("\n", " "),
            "N1": n1, "N2": n2,
            "Median1": round(np.median(s1), 4),
            "Median2": round(np.median(s2), 4),
            "U_stat": round(stat, 4),
            "P_value": pval,
            "EffectSize_RBC": round(rbc, 4),
            "Significance": stars
        })
        pair_stats[i] = (pval, stars)

    stats_df = pd.DataFrame(stats_rows)

    total_width = cur - gap
    fig_width = max(14, total_width * 1.5)
    fig, ax = plt.subplots(figsize=(fig_width, 7))

    df_plot2 = df_plot.copy()
    df_plot2["x_pos"] = df_plot2["PlotGroup"].map(x_positions)
    df_plot2 = df_plot2.dropna(subset=["x_pos"])

    def whisker_bounds(data):
        q1, q3 = np.percentile(data, 25), np.percentile(data, 75)
        iqr = q3 - q1
        lo = q1 - 1.5 * iqr; hi = q3 + 1.5 * iqr
        actual_lo = data[data >= lo].min() if (data >= lo).any() else q1
        actual_hi = data[data <= hi].max() if (data <= hi).any() else q3
        return actual_lo, actual_hi

    group_whiskers = {}
    for grp in order:
        vals = df_plot2.loc[df_plot2["PlotGroup"] == grp, "Score"].values
        if len(vals) > 0:
            group_whiskers[grp] = whisker_bounds(vals)

    PAD_TOP    = df_plot2["Score"].std() * 0.25
    PAD_BOTTOM = df_plot2["Score"].std() * 0.18

    pair_ymin_map = {}
    pair_ymax_map = {}
    for i, (g1, g2) in enumerate(pair_definitions):
        lo_vals = [group_whiskers[g][0] for g in [g1, g2] if g in group_whiskers]
        hi_vals = [group_whiskers[g][1] for g in [g1, g2] if g in group_whiskers]
        if lo_vals:
            pair_ymin_map[i] = min(lo_vals)
            pair_ymax_map[i] = max(hi_vals)

    for i, (g1, g2) in enumerate(pair_definitions):
        if g1 not in x_positions or g2 not in x_positions:
            continue
        if i not in pair_ymin_map:
            continue
        x1 = x_positions[g1]
        x2 = x_positions[g2]
        draw_pair_box(ax, x1, x2,
                      pair_ymin_map[i] - PAD_BOTTOM,
                      pair_ymax_map[i] + PAD_TOP,
                      color=pair_colors[i], alpha=0.07)

    box_width = 0.60
    for grp in order:
        sub = df_plot2[df_plot2["PlotGroup"] == grp]["Score"].values
        if len(sub) == 0:
            continue
        xp = x_positions[grp]
        color = palette[grp]

        bp = ax.boxplot(
            sub,
            positions=[xp],
            widths=box_width,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color="black", linewidth=2.0),
            boxprops=dict(facecolor=color, edgecolor="black", linewidth=0.9),
            whiskerprops=dict(color="black", linewidth=0.9),
            capprops=dict(color="black", linewidth=0.9),
            zorder=2
        )

    np.random.seed(42)
    for _, row in medians.iterrows():
        grp, cl, med = row["PlotGroup"], row["CellLine"], row["Median"]
        if grp not in x_positions:
            continue
        xp = x_positions[grp] + np.random.uniform(-0.13, 0.13)
        ax.scatter(xp, med, color=cl_palette[cl], s=55,
                   edgecolor="black", linewidth=0.6, zorder=4,
                   label=cl if grp == order[0] else None)

    score_std   = df_plot2["Score"].std()
    bracket_h   = score_std * 0.08

    for i, (g1, g2) in enumerate(pair_definitions):
        if g1 not in x_positions or g2 not in x_positions:
            continue
        if i not in pair_stats or i not in pair_ymax_map:
            continue
        pval, stars = pair_stats[i]
        if stars == "NA":
            continue
        x1 = x_positions[g1]
        x2 = x_positions[g2]
        color = pair_colors[i]
        bracket_y = pair_ymax_map[i] + PAD_TOP * 0.25
        draw_significance_bracket(ax, x1, x2,
                                  bracket_y,
                                  bracket_h,
                                  stars,
                                  color=color,
                                  fontsize=10)

    ax.axhline(0, linestyle=":", color="gray", alpha=0.45, linewidth=1)

    global_ymin = min(pair_ymin_map.values()) - PAD_BOTTOM - score_std * 0.3
    global_ymax = max(pair_ymax_map.values()) + PAD_TOP   + score_std * 0.35

    ax.set_xticks(list(x_positions.values()))
    ax.set_xticklabels(
        [g for g in x_positions.keys()],
        fontsize=8.5, rotation=0, ha="center"
    )
    ax.set_xlim(-0.75, total_width + 0.55)
    ax.set_ylim(global_ymin, global_ymax)

    ax.set_ylabel("CRUSH Compartment Score", fontsize=12)
    ax.set_xlabel("")
    ax.set_title(
        "K27me3 + Active vs Active Only: Matched Dose-Response Comparison",
        fontsize=13, pad=14
    )
    ax.tick_params(axis="x", bottom=False)
    sns.despine(ax=ax)

    for i, (g1, g2) in enumerate(pair_definitions):
        if g1 not in x_positions or g2 not in x_positions:
            continue
        x_mid = (x_positions[g1] + x_positions[g2]) / 2
        pair_label = "Pair 0\n(K27 vs Active)" if i == 0 else f"Pair {i}\n(+{i} active marks)"
        ax.text(x_mid, global_ymin + score_std * 0.05,
                pair_label,
                ha="center", va="bottom",
                fontsize=7.5, color=pair_colors[i],
                fontstyle="italic")

    handles, labels_leg = ax.get_legend_handles_labels()
    seen = {}
    uh, ul = [], []
    for h, l in zip(handles, labels_leg):
        if l not in seen:
            seen[l] = True
            uh.append(h)
            ul.append(l)

    uh.append(mpatches.Patch(facecolor="#8da0cb", edgecolor="black", label="With K27me3"))
    ul.append("With K27me3")
    uh.append(mpatches.Patch(facecolor="#fc8d62", edgecolor="black", label="No K27me3"))
    ul.append("No K27me3")

    for i, pc in enumerate(pair_colors):
        pname = "K27 only vs Active only" if i == 0 else f"K27+{i} vs {i} active"
        uh.append(mpatches.Patch(facecolor=pc, edgecolor=pc, alpha=0.5, label=pname))
        ul.append(pname)

    ax.legend(uh, ul, title="Cell Line / Group", bbox_to_anchor=(1.02, 1),
              loc="upper left", frameon=False, fontsize=7.5, title_fontsize=8.5)

    plt.tight_layout()
    svg_out = f"{args.outprefix}_interleaved_boxplot_v2.svg"
    png_out = f"{args.outprefix}_interleaved_boxplot_v2.png"
    fig.savefig(svg_out, format="svg", bbox_inches="tight")
    fig.savefig(png_out, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Interleaved plot  → {png_out}")

    stats_out = f"{args.outprefix}_matched_pairwise_stats.tsv"
    stats_df.to_csv(stats_out, sep="\t", index=False)
    print(f"Pairwise stats    → {stats_out}")

    summary = df.groupby("PlotGroup")["Score"].agg(
        N="count", Median="median", Mean="mean", Std="std",
        Q1=lambda x: x.quantile(0.25), Q3=lambda x: x.quantile(0.75),
        Pct_positive=lambda x: round((x > 0).mean() * 100, 2)
    ).reindex(order)
    summary_out = f"{args.outprefix}_pooled_summary.tsv"
    summary.to_csv(summary_out, sep="\t")
    print(f"Pooled summary    → {summary_out}")

    print("\nDone.")


if __name__ == "__main__":
    main()
