#!/usr/bin/env python3
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

import fig4E_new as f
import regulon_rank as rr

HERE = Path(__file__).resolve().parent
THR = [0, 25, 50, 75, 100, 150]


def adjusted():
    g = pd.read_csv(rr.CT / "genes_median", sep="\t", header=None, names=["gene", "e", "m"])
    g = g.drop_duplicates("gene").set_index("gene").dropna()
    ref = (np.sort(g.e.to_numpy()) + np.sort(g.m.to_numpy())) / 2
    for c in ("e", "m"):
        g[c + "q"] = ref[g[c].rank(method="first").astype(int) - 1]
    pos = pd.read_csv(rr.CT / "Loxafr3.0_gtf_lifted_to_LoxAfr3.0_HiC_capitalized.bed", sep="\t", header=None,
                      names=["chr", "start", "end", "gene"]).drop_duplicates("gene").set_index("gene")
    reg = rr.regulons(rr.LIBS)
    u = g[g.index.isin(set().union(*reg.values())) & g.index.isin(pos.index)]
    d = u["mq"] - u["eq"]
    deg = Counter(x for t, s in reg.items() for x in s - {t})
    b = pd.qcut(pd.Series([deg[x] for x in d.index], index=d.index), rr.NBIN, duplicates="drop")
    return d - d.groupby(b, observed=True).transform("median"), reg


def main():
    d, reg = adjusted()
    sign = {tf: (-1 if r.combined < 0 else 1)
            for tf, r in pd.read_csv(HERE / "fig4E_rank.tsv", sep="\t").set_index("TF").loc[f.TFS].iterrows()}
    order = [t for t in pd.read_csv(HERE / "fig4E_targets.tsv", sep="\t").TF]
    rows = []
    for tf in order:
        tg = d.index.isin(reg[tf] - {tf})
        x = sign[tf] * d
        for T in THR:
            conc, disc = (x > 0) if T == 0 else (x >= T), (x < 0) if T == 0 else (x <= -T)
            a, b_, c, e = (conc & tg).sum(), (disc & tg).sum(), (conc & ~tg).sum(), (disc & ~tg).sum()
            rows.append({"TF": tf, "direction": "elephant" if sign[tf] < 0 else "mammoth", "threshold": T,
                         "n_targets": int(tg.sum()), "n_past_T": a + b_, "concordant": a, "discordant": b_,
                         "pct_concordant": 100 * a / (a + b_), "bg_pct_concordant": 100 * c / (c + e),
                         "odds_ratio": (a * e) / (b_ * c) if b_ * c else np.nan,
                         "p_fisher": fisher_exact([[a, b_], [c, e]], alternative="greater")[1]})
    t = pd.DataFrame(rows)
    t.to_csv(HERE / "fig4E_targets_thr.tsv", sep="\t", index=False, float_format="%.4g")
    print(t.to_string(index=False))

    for notext in (False, True):
        fig, axes = plt.subplots(2, 3, figsize=(5.4, 3.4), sharex=True, sharey=True,
                                 gridspec_kw={"hspace": 0.55, "wspace": 0.15})
        for ax, tf in zip(axes.flat, order):
            x = t[t.TF == tf]
            col = f.E_COL if sign[tf] < 0 else f.M_COL
            ax.plot(x.threshold, x.pct_concordant, "-o", color=col, ms=2.5, lw=1)
            ax.plot(x.threshold, x.bg_pct_concordant, "--", color="grey", lw=0.8)
            ax.set_ylim(0, 100)
            ax.set_xticks(THR)
            ax.tick_params(labelsize=6)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if notext:
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                continue
            ax.set_title(tf, fontsize=7, style="italic", color=col)
            for _, s in x.iterrows():
                ax.text(s.threshold, 4, f"{int(s.n_past_T)}", ha="center", fontsize=4.5, color="dimgrey")
                if s.p_fisher < 0.05:
                    ax.text(s.threshold, s.pct_concordant + 4, "*", ha="center", fontsize=7)
        if not notext:
            fig.supxlabel("|ΔCRUSH| threshold", fontsize=8, y=0.0)
            fig.supylabel("Targets in the TF's direction (%)", fontsize=8, x=0.04)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"fig4E_targets_thr{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
