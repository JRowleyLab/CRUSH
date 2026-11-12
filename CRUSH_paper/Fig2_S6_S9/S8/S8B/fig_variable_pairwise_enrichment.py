import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

BASE    = "fig2_Encode"
OUT_DIR = f"{BASE}/hmm/Feb18_2026/Feb23_Outie"
MASTER  = f"{OUT_DIR}/master_table_7lines.tsv.gz"

STATES = [
    "1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk",
    "6_EnhG","7_Enh","8_ZNF/Rpts","9_Het",
    "10_TssBiv","11_BivFlnk","12_EnhBiv",
    "13_ReprPC","14_ReprPCWk","15_Quies"
]

STATE_COLORS = {
    "1_TssA":     "#FF0000",
    "2_TssAFlnk": "#FF6E00",
    "3_TxFlnk":   "#32CD32",
    "4_Tx":       "#008000",
    "5_TxWk":     "#006400",
    "6_EnhG":     "#C2E105",
    "7_Enh":      "#FFFF00",
    "8_ZNF/Rpts": "#66CDAA",
    "9_Het":      "#8A91D0",
    "10_TssBiv":  "#CD5C5C",
    "11_BivFlnk": "#E9967A",
    "12_EnhBiv":  "#BDB76B",
    "13_ReprPC":  "#808080",
    "14_ReprPCWk":"#C0C0C0",
    "15_Quies":   "#DDDDDD",
}

TC = "#111111"
SC = "#AAAAAA"
GC = "#E8E8E8"

print(f"Loading {MASTER}...")
df = pd.read_csv(MASTER, sep="\t", compression="gzip")
cats = df["category"].value_counts()
print(f"  {len(df):,} rows  |  {cats.to_dict()}")

var  = df[df["category"] == "Variable"]
stA  = df[df["category"] == "StableA"]
stB  = df[df["category"] == "StableB"]
stAB = df[df["category"].isin(["StableA", "StableB"])]

n_var  = len(var)
n_stA  = len(stA)
n_stB  = len(stB)
n_stAB = len(stAB)

print(f"  Variable: {n_var:,} | StableA: {n_stA:,} | StableB: {n_stB:,} | StableA+B: {n_stAB:,}")

def log2_enr(query_df, ref_df, states=STATES):
    """
    For each state: log2( freq_in_query / freq_in_ref )
    Positive = enriched in query (Variable)
    Negative = depleted in query (enriched in reference)
    """
    enr = {}
    for s in states:
        fq = (query_df["coarse_state"] == s).sum() / max(len(query_df), 1)
        fr = (ref_df["coarse_state"]   == s).sum() / max(len(ref_df),   1)
        enr[s] = np.log2((fq + 1e-6) / (fr + 1e-6))
    return pd.Series(enr).reindex(states)

enr_A   = log2_enr(var, stA)
enr_B   = log2_enr(var, stB)
enr_AB  = log2_enr(var, stAB)

print("\nEnrichment computed:")
print(pd.DataFrame({"vs_StableA": enr_A, "vs_StableB": enr_B,
                    "vs_StableAB": enr_AB}).round(3).to_string())

fig, axes = plt.subplots(1, 3, figsize=(20, 8), facecolor="white")
fig.patch.set_facecolor("white")
fig.subplots_adjust(wspace=0.08, left=0.10, right=0.97, top=0.88, bottom=0.06)

panels = [
    (axes[0], enr_A,  "Variable vs StableA",
     f"log2(Variable / StableA)\n"
     f"Variable n={n_var:,} | StableA n={n_stA:,}",
     "Enriched in Variable →",
     "← Enriched in StableA"),

    (axes[1], enr_B,  "Variable vs StableB",
     f"log2(Variable / StableB)\n"
     f"Variable n={n_var:,} | StableB n={n_stB:,}",
     "Enriched in Variable →",
     "← Enriched in StableB"),

    (axes[2], enr_AB, "Variable vs StableA + StableB",
     f"log2(Variable / [StableA+StableB])\n"
     f"Variable n={n_var:,} | Stable n={n_stAB:,}",
     "Enriched in Variable →",
     "← Enriched in Stable"),
]

for ax, enr_s, title, xlabel, lbl_pos, lbl_neg in panels:
    ax.set_facecolor("white")

    pass

    y_pos = np.arange(len(STATES))
    vals  = enr_s.values

    bar_colors = []
    for s, v in zip(STATES, vals):
        base = np.array(mcolors.to_rgb(STATE_COLORS[s]))
        if v > 0:
            mix = np.array(mcolors.to_rgb("#CC2222"))
            c   = tuple(0.5 * base + 0.5 * mix)
        else:
            mix = np.array(mcolors.to_rgb("#1A4FA0"))
            c   = tuple(0.45 * base + 0.55 * mix)
        bar_colors.append(c)

    bars = ax.barh(y_pos, vals,
                   color=bar_colors,
                   edgecolor="white", linewidth=0.5,
                   height=0.70, zorder=3)

    for i, (v, bar) in enumerate(zip(vals, bars)):
        if abs(v) < 0.05:
            continue
        ha   = "left"  if v > 0 else "right"
        xoff = 0.05 if v > 0 else -0.05
        ax.text(v + xoff, i, f"{v:+.2f}",
                va="center", ha=ha, fontsize=7,
                color=TC, fontweight="bold", zorder=5)

    ax.axvline(0, color="#444444", lw=1.1, ls="--", zorder=4, alpha=0.8)

    ax.set_yticks(y_pos)
    if ax is axes[0]:
        ax.set_yticklabels(STATES, fontsize=9, color=TC)
    else:
        ax.set_yticklabels([""] * len(STATES))

    ax.invert_yaxis()

    ax.set_xlabel(xlabel, fontsize=8.5, color="#444444", labelpad=5)
    ax.set_title(title, fontsize=11.5, fontweight="bold", color=TC, pad=10)

    ax.grid(axis="x", color=GC, lw=0.6, zorder=0)
    for sp in ax.spines.values():
        sp.set_edgecolor(SC); sp.set_linewidth(0.7)
    ax.tick_params(axis="x", colors=TC, labelsize=8)
    ax.tick_params(axis="y", colors=TC)

    xlim = ax.get_xlim()
    ax.text(xlim[1] * 0.98, -0.8, lbl_pos,
            ha="right", va="top", fontsize=7.5, color="#AA2222",
            style="italic", transform=ax.transData)
    ax.text(xlim[0] * 0.98, -0.8, lbl_neg,
            ha="left", va="top", fontsize=7.5, color="#1A4FA0",
            style="italic", transform=ax.transData)

    xmax = max(abs(xlim[0]), abs(xlim[1]))
    ax.axvspan(0,     xmax, alpha=0.025, color="#CC2222", zorder=0)
    ax.axvspan(-xmax, 0,    alpha=0.025, color="#1A4FA0", zorder=0)
    ax.set_xlim(-xmax, xmax)

all_xlim_max = max(max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1])) for ax in axes)
for ax in axes:
    ax.set_xlim(-all_xlim_max, all_xlim_max)

group_boundaries = [4.5, 8.5, 12.5]
for ax in axes:
    for y in group_boundaries:
        ax.axhline(y, color="#CCCCCC", lw=0.8, ls=":", zorder=2)

group_labels = [
    (2.0,  "Active\nPromoters\n& TxFlnk"),
    (6.5,  "Enhancers\n& ZNF/Het"),
    (10.5, "Bivalent"),
    (13.5, "Polycomb\n& Quies"),
]
for y, lbl in group_labels:
    axes[0].text(-all_xlim_max - 0.15, y, lbl,
                 ha="right", va="center", fontsize=7, color="#777777",
                 style="italic", transform=axes[0].transData)

fig.suptitle(
    "Chromatin State Enrichment in Variable Bins — Pairwise Comparisons with Stable Compartments",
    fontsize=13, fontweight="bold", color=TC, y=0.96
)
fig.text(
    0.5, 0.005,
    "log2(freq_Variable / freq_Reference) per chromatin state  |  "
    "Positive = enriched in Variable bins  |  Negative = enriched in reference category  |  "
    "7 cell lines, 1kb bins, Roadmap 15-state model",
    ha="center", va="bottom", fontsize=7.5, color="#888888", style="italic"
)

out_svg = os.path.join(OUT_DIR, "fig_variable_pairwise_enrichment.svg")
fig.savefig(out_svg, bbox_inches="tight", facecolor="white", dpi=180)
plt.close()
print(f"\nSaved → {out_svg}")

tbl = pd.DataFrame({
    "Variable_vs_StableA":       enr_A.round(4),
    "Variable_vs_StableB":       enr_B.round(4),
    "Variable_vs_StableA_and_B": enr_AB.round(4),
})
tbl_path = os.path.join(OUT_DIR, "fig_variable_pairwise_enrichment_values.tsv")
tbl.to_csv(tbl_path, sep="\t")
print(f"Values → {tbl_path}")
