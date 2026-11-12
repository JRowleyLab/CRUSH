import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick

BASE    = "fig2_Encode"
IN_DIR  = f"{BASE}/hmm/Feb18_2026"
OUT_DIR = f"{BASE}/hmm/Feb18_2026/Feb23_Outie"
os.makedirs(OUT_DIR, exist_ok=True)

CELL_LINES = ["A549", "HeLa", "HepG2", "HUVEC", "IMR90", "K562", "HCT116"]

CL_FILES = {
    "A549":   "a549_CR_HMM_1kb.bed",
    "HeLa":   "helaS3_CR_HMM_1kb.bed",
    "HepG2":  "hepg2_CR_HMM_1kb.bed",
    "HUVEC":  "huvec_CR_HMM_1kb.bed",
    "IMR90":  "imr90_CR_HMM_1kb.bed",
    "K562":   "k562_CR_HMM_1kb.bed",
    "HCT116": "hct116_CR_HMM_1kb_re.bed",
}

COARSE_MAP = {s: s for s in [
    "1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk",
    "6_EnhG","7_Enh","8_ZNF/Rpts","9_Het",
    "10_TssBiv","11_BivFlnk","12_EnhBiv",
    "13_ReprPC","14_ReprPCWk","15_Quies"
]}
COARSE_STATES = list(COARSE_MAP.keys())
STATES_NO_QUIES = [s for s in COARSE_STATES if s != "15_Quies"]

COARSE_COLORS = {
    "1_TssA":    "#FF0000", "2_TssAFlnk":"#FF6E00", "3_TxFlnk":  "#32CD32",
    "4_Tx":      "#008000", "5_TxWk":    "#006400", "6_EnhG":    "#C2E105",
    "7_Enh":     "#FFFF00", "8_ZNF/Rpts":"#66CDAA", "9_Het":     "#8A91D0",
    "10_TssBiv": "#CD5C5C", "11_BivFlnk":"#E9967A", "12_EnhBiv": "#BDB76B",
    "13_ReprPC": "#808080", "14_ReprPCWk":"#C0C0C0","15_Quies":  "#EBEBEB"
}

CL_COLORS = {
    "A549":"#E63946","HeLa":"#457B9D","HepG2":"#2A9D8F",
    "HUVEC":"#E9C46A","IMR90":"#F4A261","K562":"#6A4C93","HCT116":"#FF6B9D"
}

CAT_COLORS  = {"StableA":"#1565C0","StableB":"#C62828","Variable":"#2E7D32"}
CATEGORIES  = ["StableA","StableB","Variable"]

A_EXPECTED = {"1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk",
              "6_EnhG","7_Enh","10_TssBiv","11_BivFlnk","12_EnhBiv"}
B_EXPECTED = {"8_ZNF/Rpts","9_Het","13_ReprPC","14_ReprPCWk","15_Quies"}

BG   = "white"
TC   = "#111111"
GC   = "#DDDDDD"
SC   = "#AAAAAA"

def style_ax(ax):
    ax.set_facecolor(BG)
    ax.tick_params(colors=TC, labelsize=8)
    for sp in ax.spines.values():
        sp.set_edgecolor(SC); sp.set_linewidth(0.7)
    ax.xaxis.label.set_color(TC)
    ax.yaxis.label.set_color(TC)
    ax.title.set_color(TC)
    ax.grid(axis="y", color=GC, lw=0.5, zorder=0)

def save_fig(fig, name):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, bbox_inches="tight", facecolor="white", dpi=180)
    plt.close(fig)
    print(f"  Saved → {path}")

def normalize_state(s):
    s = str(s).strip()
    if s in COARSE_MAP: return s
    lookup = {k.split("_",1)[1]: k for k in COARSE_MAP}
    return lookup.get(s, None)

def build_master(in_dir):
    master_path = os.path.join(OUT_DIR, "master_table_7lines.tsv.gz")
    if os.path.exists(master_path):
        df = pd.read_csv(master_path, sep="\t", compression="gzip")
        if df["coarse_state"].nunique() <= 7:
            print("Old 7-state master — rebuilding...")
            os.remove(master_path)
        else:
            print(f"Loaded master: {len(df):,} rows | {df['bin'].nunique():,} bins")
            print(df["category"].value_counts())
            return df

    print("Building master table (7 cell lines)...")
    records = []
    for cl in CELL_LINES:
        path = os.path.join(in_dir, CL_FILES[cl])
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found — skipping {cl}")
            continue
        with open(path) as fh:
            first = fh.readline()
        has_header = first.strip().split()[0] == "chrom"
        cdf = pd.read_csv(path, sep="\t",
                          header=0 if has_header else None,
                          names=None if has_header else
                                ["chrom","start","end","eigen","hmm_state"])
        cdf.columns = ["chrom","start","end","eigen","hmm_state"]
        cdf = cdf[~cdf["chrom"].astype(str).str.startswith("track")]
        cdf = cdf.dropna(subset=["eigen","hmm_state"])
        cdf = cdf[cdf["hmm_state"].astype(str) != "."]
        cdf["hmm_state"]    = cdf["hmm_state"].apply(normalize_state)
        cdf = cdf[cdf["hmm_state"].notna()]
        cdf["coarse_state"] = cdf["hmm_state"].map(COARSE_MAP)
        cdf = cdf[cdf["coarse_state"].notna()]
        cdf["eigen"]        = pd.to_numeric(cdf["eigen"], errors="coerce")
        cdf = cdf.dropna(subset=["eigen"])
        cdf["comp_call"]    = cdf["eigen"].apply(lambda e: "A" if e > 0 else "B")
        cdf["cell_line"]    = cl
        cdf["bin"]          = cdf["chrom"].astype(str) + ":" + cdf["start"].astype(str)
        records.append(cdf[["bin","cell_line","eigen","comp_call","hmm_state","coarse_state"]])
        print(f"  {cl}: {len(cdf):,} bins | "
              f"A={(cdf['comp_call']=='A').sum():,} | B={(cdf['comp_call']=='B').sum():,}")

    df = pd.concat(records, ignore_index=True)
    a_counts = (df.groupby("bin")["comp_call"]
                  .apply(lambda x: (x=="A").sum()).rename("a_count"))
    df = df.merge(a_counts, on="bin")
    df["category"] = "Variable"
    df.loc[df["a_count"] == len(CELL_LINES), "category"] = "StableA"
    df.loc[df["a_count"] == 0,               "category"] = "StableB"
    print(f"\nTotal: {len(df):,} rows | {df['bin'].nunique():,} unique bins")
    print(df["category"].value_counts())
    df.to_csv(master_path, sep="\t", index=False, compression="gzip")
    print(f"Saved → {master_path}")
    return df

def fig1_composition(df, states, suffix="", title_suffix=""):
    tag = "noQ" if "noQ" in suffix else "full"
    print(f"\n[Fig 1{suffix}] Composition ({tag})...")
    fig, axes = plt.subplots(1, 3, figsize=(22, 8), facecolor=BG)
    fig.patch.set_facecolor(BG)
    fig.subplots_adjust(wspace=0.06)

    dfp = df if not states else df[df["coarse_state"].isin(states)]
    n_cl = len(CELL_LINES)

    for ax, cat in zip(axes, CATEGORIES):
        style_ax(ax)
        ax.grid(visible=False)
        sub  = dfp[dfp["category"]==cat]
        bw   = 0.65

        frac = (sub["coarse_state"].value_counts(normalize=True)
                  .reindex(states, fill_value=0))
        bot  = 0.0
        for s in states:
            h = frac[s]
            if h < 0.001: bot += h; continue
            ax.bar(0, h, bw, bottom=bot, color=COARSE_COLORS[s],
                   edgecolor="#888888", lw=0.3, zorder=3)
            if h > 0.045:
                ax.text(0, bot+h/2, s.split("_",1)[1][:5],
                        ha="center", va="center", fontsize=5.5,
                        fontweight="bold", color="black", zorder=5)
            bot += h

        for i, cl in enumerate(CELL_LINES):
            xi   = 2 + i * 0.72
            frc  = (sub[sub["cell_line"]==cl]["coarse_state"]
                      .value_counts(normalize=True)
                      .reindex(states, fill_value=0))
            bot  = 0.0
            for s in states:
                h = frc[s]
                if h < 0.001: bot += h; continue
                base  = np.array(matplotlib.colors.to_rgb(COARSE_COLORS[s]))
                tint  = np.array(matplotlib.colors.to_rgb(CL_COLORS[cl]))
                blend = tuple(0.72*base + 0.28*tint)
                ax.bar(xi, h, bw*0.85, bottom=bot, color=blend,
                       edgecolor=CL_COLORS[cl], lw=0.9, zorder=3)
                bot += h

        x_max = 2 + n_cl*0.72
        ax.set_xlim(-0.65, x_max)
        ax.set_ylim(0, 1.05)
        ax.set_xticks([0] + [2+i*0.72 for i in range(n_cl)])
        ax.set_xticklabels(["All"]+CELL_LINES, rotation=45, ha="right",
                            fontsize=8, color=TC)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
        ax.tick_params(axis="both", colors=TC)
        for sp in ax.spines.values(): sp.set_edgecolor(SC); sp.set_linewidth(0.7)
        ax.set_title(cat, color=CAT_COLORS[cat], fontsize=14, fontweight="bold", pad=8)
        ax.axvline(0.48, color=SC, lw=1.0, ls="--", zorder=1)
        if cat == "StableA":
            ylabel = "Fraction of bins" if not title_suffix else f"Fraction of bins {title_suffix}"
            ax.set_ylabel(ylabel, fontsize=10, color=TC)
        else:
            ax.set_yticklabels([])

    p1 = [mpatches.Patch(color=COARSE_COLORS[s], label=s, 
                         edgecolor="#888", linewidth=0.5) for s in states]
    p2 = [mpatches.Patch(color=CL_COLORS[c], label=c,
                         edgecolor="#888", linewidth=0.5) for c in CELL_LINES]
    l1 = fig.legend(handles=p1, title="ChromHMM State", loc="center right",
                    bbox_to_anchor=(0.99, 0.68), fontsize=7, title_fontsize=8,
                    facecolor=BG, labelcolor=TC, edgecolor=SC, framealpha=1,
                    ncol=1)
    l1.get_title().set_color(TC)
    l2 = fig.legend(handles=p2, title="Cell Line", loc="center right",
                    bbox_to_anchor=(0.99, 0.22), fontsize=8, title_fontsize=9,
                    facecolor=BG, labelcolor=TC, edgecolor=SC, framealpha=1)
    l2.get_title().set_color(TC)
    ttl = f"ChromHMM State Composition by Compartment Category{title_suffix}"
    fig.suptitle(ttl, color=TC, fontsize=13, fontweight="bold", y=0.97)
    save_fig(fig, f"fig1{suffix}_composition.svg")

def fig2_enrichment(df, states, suffix="", title_suffix=""):
    print(f"\n[Fig 2{suffix}] Enrichment heatmap...")
    dfp = df if not suffix.endswith("noQ") else \
          df[df["coarse_state"] != "15_Quies"]
    rows = []
    for cl in CELL_LINES:
        sub = dfp[dfp["cell_line"]==cl]
        a   = sub[sub["comp_call"]=="A"]
        b   = sub[sub["comp_call"]=="B"]
        for s in states:
            fA = (a["coarse_state"]==s).sum() / max(len(a),1)
            fB = (b["coarse_state"]==s).sum() / max(len(b),1)
            rows.append({"cell_line":cl,"state":s,
                         "log2FC": np.log2((fA+1e-6)/(fB+1e-6))})
    hdf = (pd.DataFrame(rows)
             .pivot(index="state", columns="cell_line", values="log2FC")
             .reindex(index=states, columns=CELL_LINES))

    fig_h = max(6, len(states)*0.52)
    fig, ax = plt.subplots(figsize=(10, fig_h), facecolor=BG)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    vmax = np.abs(hdf.values).max()
    im   = ax.imshow(hdf.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(CELL_LINES)))
    ax.set_xticklabels(CELL_LINES, rotation=45, ha="right", color=TC, fontsize=10)
    ax.set_yticks(range(len(states)))
    ax.set_yticklabels(states, color=TC, fontsize=9)
    for i in range(len(states)):
        for j in range(len(CELL_LINES)):
            v = hdf.values[i,j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    fontweight="bold",
                    color="black" if abs(v)<vmax*0.6 else "white")
    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("log2(A/B enrichment)", color=TC, fontsize=9)
    plt.setp(cb.ax.yaxis.get_ticklabels(), color=TC, fontsize=8)
    cb.ax.yaxis.set_tick_params(color=TC)
    cb.outline.set_edgecolor(SC)
    ttl = f"Chromatin State Enrichment: A vs B per Cell Line{title_suffix}"
    ax.set_title(ttl, color=TC, fontsize=12, fontweight="bold", pad=12)
    for sp in ax.spines.values(): sp.set_edgecolor(SC); sp.set_linewidth(0.7)
    fig.tight_layout()
    save_fig(fig, f"fig2{suffix}_enrichment.svg")
    hdf.to_csv(os.path.join(OUT_DIR, f"fig2{suffix}_enrichment_table.tsv"), sep="\t")

def fig3_switching(df, states, suffix="", title_suffix=""):
    print(f"\n[Fig 3{suffix}] Switching analysis...")
    var = df[df["category"]=="Variable"]
    if suffix.endswith("noQ"):
        var = var[var["coarse_state"] != "15_Quies"]
    sta = df[df["category"]=="StableA"]
    if suffix.endswith("noQ"):
        sta = sta[sta["coarse_state"] != "15_Quies"]

    fig, axes = plt.subplots(1, 3, figsize=(22, 7), facecolor=BG)
    fig.patch.set_facecolor(BG)

    ax = axes[0]; style_ax(ax)
    spec = (var.drop_duplicates(["bin","cell_line"])
               .groupby("bin")["a_count"].first()
               .value_counts().sort_index()
               .reindex(range(1, len(CELL_LINES)), fill_value=0))
    cols = plt.cm.RdYlGn(np.linspace(0.15, 0.85, len(spec)))
    ax.bar(spec.index, spec.values, color=cols, edgecolor="white", lw=0.6, zorder=3)
    for x,y in zip(spec.index, spec.values):
        ax.text(x, y + max(spec.values)*0.01, f"{y:,}", ha="center",
                va="bottom", color=TC, fontsize=8)
    ax.set_xlabel(f"# cell lines calling A (of {len(CELL_LINES)})", color=TC, fontsize=10)
    ax.set_ylabel("Number of bins", color=TC, fontsize=10)
    ax.set_title(f"Switching Spectrum{title_suffix}",
                 color=CAT_COLORS["Variable"], fontsize=11, fontweight="bold")

    ax = axes[1]; style_ax(ax); ax.grid(visible=False)
    afr = (var[var["comp_call"]=="A"]["coarse_state"]
             .value_counts(normalize=True).reindex(states, fill_value=0))
    bfr = (var[var["comp_call"]=="B"]["coarse_state"]
             .value_counts(normalize=True).reindex(states, fill_value=0))
    x   = np.arange(len(states))
    w   = 0.35
    ax.bar(x-w/2, afr.values, w, label="A calls",
           color=[COARSE_COLORS[s] for s in states],
           edgecolor="#2166AC", lw=1.0, alpha=0.95)
    ax.bar(x+w/2, bfr.values, w, label="B calls",
           color=[COARSE_COLORS[s] for s in states],
           edgecolor="#D6604D", lw=1.0, alpha=0.55)
    ax.set_xticks(x)
    ax.set_xticklabels([s.split("_",1)[1][:6] for s in states],
                       rotation=45, ha="right", color=TC, fontsize=7.5)
    ax.set_ylabel("Fraction of Variable bins", color=TC, fontsize=10)
    ax.set_title(f"HMM States: A vs B in Variable Bins{title_suffix}",
                 color=CAT_COLORS["Variable"], fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, facecolor=BG, edgecolor=SC, framealpha=0.9)
    ax.tick_params(axis="x", colors=TC, pad=6, rotation=45)

    ax = axes[2]; style_ax(ax); ax.grid(axis="x", color=GC, lw=0.5)
    ax.grid(axis="y", visible=False)
    enr = {}
    for s in states:
        fv = (var["coarse_state"]==s).sum() / max(len(var),1)
        fa = (sta["coarse_state"]==s).sum() / max(len(sta),1)
        enr[s] = np.log2((fv+1e-6)/(fa+1e-6))
    enr_s = pd.Series(enr).reindex(states)
    cols  = ["#C62828" if v>0 else "#1565C0" for v in enr_s.values]
    ax.barh(range(len(states)), enr_s.values, color=cols,
            edgecolor="white", lw=0.3, zorder=3)
    ax.set_yticks(range(len(states)))
    ax.set_yticklabels(states, color=TC, fontsize=8)
    ax.axvline(0, color=TC, lw=0.8, ls="--", zorder=4)
    ax.set_xlabel("log2(Variable / StableA)", color=TC, fontsize=10)
    ax.set_title(f"State Enrichment: Variable vs StableA{title_suffix}",
                 color=CAT_COLORS["Variable"], fontsize=11, fontweight="bold")

    fig.suptitle(f"Variable Bin Switching Analysis{title_suffix}",
                 color=TC, fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    save_fig(fig, f"fig3{suffix}_switching.svg")

def fig4_drivers(df):
    print("\n[Fig 4] Switching drivers...")
    var    = df[df["category"]=="Variable"]
    n_cl   = len(CELL_LINES)
    uniq_a = var[(var["a_count"]==1)        & (var["comp_call"]=="A")]
    uniq_b = var[(var["a_count"]==n_cl-1)   & (var["comp_call"]=="B")]

    fig, axes = plt.subplots(1, 2, figsize=(16, 7), facecolor=BG)
    fig.patch.set_facecolor(BG)
    for ax, data, label in [
        (axes[0], uniq_a, f"Unique A caller (1 of {n_cl} lines)"),
        (axes[1], uniq_b, f"Unique B caller ({n_cl-1} of {n_cl} lines A)")
    ]:
        style_ax(ax); ax.grid(visible=False)
        cnt  = data.groupby(["cell_line","coarse_state"]).size().unstack(fill_value=0)
        cnt  = cnt.reindex(columns=COARSE_STATES, fill_value=0)
        for cl in CELL_LINES:
            if cl not in cnt.index:
                cnt.loc[cl] = 0
        cnt  = cnt.reindex(CELL_LINES, fill_value=0)
        nrm  = cnt.div(cnt.sum(axis=1).replace(0,1), axis=0)
        x    = np.arange(len(CELL_LINES))
        bot  = np.zeros(len(CELL_LINES))
        for s in COARSE_STATES:
            vals = nrm[s].values if s in nrm.columns else np.zeros(len(x))
            ax.bar(x, vals, 0.6, bottom=bot, color=COARSE_COLORS[s],
                   edgecolor="white", lw=0.3, zorder=3)
            bot += vals
        ax.set_xticks(x)
        ax.set_xticklabels(CELL_LINES, rotation=45, ha="right",
                           fontsize=9, color=TC, rotation_mode="anchor")
        ax.tick_params(axis="x", pad=6, colors=TC)
        ax.tick_params(axis="y", colors=TC)
        ax.set_ylabel("Fraction of switching bins", color=TC, fontsize=10)
        ax.set_title(label, color=CAT_COLORS["Variable"], fontsize=11, fontweight="bold")
        ax.set_ylim(0, 1.05)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
        for sp in ax.spines.values(): sp.set_edgecolor(SC); sp.set_linewidth(0.7)

    patches = [mpatches.Patch(color=COARSE_COLORS[s], label=s,
                               edgecolor="#888", lw=0.5) for s in COARSE_STATES]
    leg = fig.legend(handles=patches, title="ChromHMM State", loc="lower center",
                     bbox_to_anchor=(0.5, -0.04), ncol=8, fontsize=7.5,
                     title_fontsize=9, facecolor=BG, labelcolor=TC,
                     edgecolor=SC, framealpha=1)
    leg.get_title().set_color(TC)
    fig.suptitle("HMM States Driving Cell-Line-Specific Compartment Switching",
                 color=TC, fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    save_fig(fig, "fig4_switching_drivers.svg")

def fig5_jaccard(df):
    print("\n[Fig 5] Pairwise Jaccard...")
    a_bins = {cl: set(df[(df["cell_line"]==cl)&(df["comp_call"]=="A")]["bin"])
              for cl in CELL_LINES}
    jacc = pd.DataFrame(index=CELL_LINES, columns=CELL_LINES, dtype=float)
    for c1 in CELL_LINES:
        for c2 in CELL_LINES:
            i = len(a_bins[c1] & a_bins[c2])
            u = len(a_bins[c1] | a_bins[c2])
            jacc.loc[c1,c2] = i/u if u else 0.0

    fig, ax = plt.subplots(figsize=(8, 7), facecolor=BG)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    vals = jacc.values.astype(float)
    im   = ax.imshow(vals, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(CELL_LINES)))
    ax.set_xticklabels(CELL_LINES, rotation=45, ha="right", color=TC, fontsize=10)
    ax.set_yticks(range(len(CELL_LINES)))
    ax.set_yticklabels(CELL_LINES, color=TC, fontsize=10)
    for i in range(len(CELL_LINES)):
        for j in range(len(CELL_LINES)):
            v = vals[i,j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=9,
                    fontweight="bold", color="black" if v>0.45 else "white")
    cb = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
    cb.set_label("Jaccard (A-bin sets)", color=TC, fontsize=9)
    plt.setp(cb.ax.yaxis.get_ticklabels(), color=TC, fontsize=8)
    cb.ax.yaxis.set_tick_params(color=TC)
    cb.outline.set_edgecolor(SC)
    for sp in ax.spines.values(): sp.set_edgecolor(SC); sp.set_linewidth(0.7)
    ax.set_title("Pairwise A-Compartment Similarity (Jaccard)",
                 color=TC, fontsize=12, fontweight="bold", pad=12)
    fig.tight_layout()
    save_fig(fig, "fig5_jaccard.svg")
    jacc.to_csv(os.path.join(OUT_DIR,"fig5_jaccard_table.tsv"),sep="\t")

def fig6_concordance(df):
    print("\n[Fig 6] Concordance...")
    df = df.copy()
    df["concordant"] = (
        ((df["comp_call"]=="A") & (df["coarse_state"].isin(A_EXPECTED))) |
        ((df["comp_call"]=="B") & (df["coarse_state"].isin(B_EXPECTED)))
    )
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor=BG)
    fig.patch.set_facecolor(BG)

    ax = axes[0]; style_ax(ax); ax.grid(visible=False)
    disc = (df.groupby(["cell_line","category"])["concordant"]
              .apply(lambda x: (~x).mean()).reset_index())
    disc.columns = ["cell_line","category","disc_rate"]
    x  = np.arange(len(CELL_LINES))
    w  = 0.24
    for i, cat in enumerate(CATEGORIES):
        vals = [disc[(disc["cell_line"]==cl)&(disc["category"]==cat)]["disc_rate"].values
                for cl in CELL_LINES]
        vals = [v[0] if len(v)>0 else 0 for v in vals]
        ax.bar(x+(i-1)*w, vals, w, label=cat, color=CAT_COLORS[cat],
               edgecolor="white", lw=0.4, alpha=0.88, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(CELL_LINES, rotation=45, ha="right", color=TC,
                       fontsize=9, rotation_mode="anchor")
    ax.tick_params(axis="x", pad=8, colors=TC)
    ax.tick_params(axis="y", colors=TC)
    ax.set_ylabel("Discordance rate", color=TC, fontsize=10)
    ax.set_title("HMM–Compartment Discordance per Cell Line",
                 color=TC, fontsize=11, fontweight="bold")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
    ax.legend(fontsize=9, facecolor=BG, edgecolor=SC, framealpha=0.9)
    for sp in ax.spines.values(): sp.set_edgecolor(SC); sp.set_linewidth(0.7)

    ax = axes[1]; style_ax(ax); ax.grid(visible=False)
    dv  = df[(df["category"]=="Variable")&(~df["concordant"])]
    cnt = dv.groupby(["coarse_state","comp_call"]).size().unstack(fill_value=0)
    nrm = cnt.div(cnt.sum(axis=1).replace(0,1), axis=0)
    x2  = np.arange(len(cnt.index))
    w2  = 0.35
    for j,(comp,col) in enumerate([("A","#2166AC"),("B","#D6604D")]):
        vals = [nrm.loc[s,comp] if comp in nrm.columns and s in nrm.index else 0
                for s in cnt.index]
        ax.bar(x2+j*w2-w2/2, vals, w2, color=col, edgecolor="white",
               lw=0.4, alpha=0.85, label=f"{comp} call", zorder=3)
    ax.set_xticks(x2)
    ax.set_xticklabels(cnt.index, rotation=45, ha="right", color=TC,
                       fontsize=8, rotation_mode="anchor")
    ax.tick_params(axis="x", pad=8, colors=TC)
    ax.tick_params(axis="y", colors=TC)
    ax.set_ylabel("Fraction of discordant bins", color=TC, fontsize=10)
    ax.set_title("Discordant States in Variable Bins",
                 color=TC, fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, facecolor=BG, edgecolor=SC, framealpha=0.9)
    for sp in ax.spines.values(): sp.set_edgecolor(SC); sp.set_linewidth(0.7)

    fig.suptitle("Chromatin State – Compartment Concordance",
                 color=TC, fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    save_fig(fig, "fig6_concordance.svg")

if __name__ == "__main__":
    print("="*60)
    print("CRUSH Analysis — 7 Cell Lines | White Background")
    print("="*60)

    df = build_master(IN_DIR)

    print("\n── Full 15-state figures ──")
    fig1_composition(df, COARSE_STATES)
    fig2_enrichment(df, COARSE_STATES)
    fig3_switching(df, COARSE_STATES)
    fig4_drivers(df)
    fig5_jaccard(df)
    fig6_concordance(df)

    print("\n── No-Quiescent figures ──")
    fig1_composition(df, STATES_NO_QUIES, suffix="_noQ",
                     title_suffix=" (Quiescent excluded)")
    fig2_enrichment(df, STATES_NO_QUIES, suffix="_noQ",
                    title_suffix=" (Quiescent excluded)")
    fig3_switching(df, STATES_NO_QUIES, suffix="_noQ",
                   title_suffix="\n(Quiescent excluded)")

    print(f"\n✓ All 9 figures saved → {OUT_DIR}")
