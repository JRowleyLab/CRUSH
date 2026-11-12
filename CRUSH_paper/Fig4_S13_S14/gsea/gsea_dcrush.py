#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "pylib"))

import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from scipy.special import expit

CT = Path("comptracks")
COLL = {"go": ["GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
               "MSigDB_Hallmark_2020"],
        "pathways": ["KEGG_2021_Human", "Reactome_Pathways_2024", "WikiPathways_2024_Human"],
        "gofull": ["GO_BP_full"]}
WHICH = sys.argv[1] if len(sys.argv) > 1 else "go"
LIBS = COLL[WHICH]
TAG = "" if WHICH == "go" else f"{WHICH}_"
VALUES = sys.argv[2].split(",") if len(sys.argv) > 2 else None
BIN, CLUST, NTOP, TAU, TAU_NARROW = 5000, 0.25, 10, 50, 15
E_COL, M_COL = "#b5307a", "#2f6fe0"


def qvals(e, m):
    ref = (np.sort(e.to_numpy()) + np.sort(m.to_numpy())) / 2
    return ref[e.rank(method="first").astype(int) - 1], ref[m.rank(method="first").astype(int) - 1]


def qnorm(e, m):
    eq, mq = qvals(e, m)
    return mq - eq


def gmt(lib):
    sets = {}
    for line in open(HERE / "ref" / f"{lib}.gmt"):
        f = line.rstrip("\n").split("\t")
        sets[f[0]] = [g for g in f[2:] if g]
    return sets


def values(pos):
    g = pd.read_csv(CT / "genes_median", sep="\t", header=None, names=["gene", "e", "m"])
    g = g.drop_duplicates("gene").set_index("gene").dropna()
    out = {"gene": pd.Series(qnorm(g.e, g.m), index=g.index)}
    eq, mq = qvals(g.e, g.m)
    out["state"] = pd.Series(expit(mq / TAU) - expit(eq / TAU), index=g.index)
    out["stateA"] = out["state"][(eq > 0) | (mq > 0)]
    out["stateAnarrow"] = pd.Series(expit(mq / TAU_NARROW) - expit(eq / TAU_NARROW), index=g.index)[(eq > 0) | (mq > 0)]
    bw = {k: pyBigWig.open(str(CT / f"{k}_5kb_CR.bw")) for k in ("elephant", "mammoth")}
    rows = {}
    for gene, (c, s, e) in pos[["chr", "start", "end"]].iterrows():
        v = []
        for k in ("elephant", "mammoth"):
            try:
                ends = [bw[k].stats(c, max(0, p // BIN * BIN), p // BIN * BIN + BIN)[0] for p in (s, e - 1)]
            except RuntimeError:
                ends = [None]
            v.append(np.nan if None in ends else np.mean(ends))
        rows[gene] = v
    t = pd.DataFrame.from_dict(rows, orient="index", columns=["e", "m"]).dropna()
    out["tss"] = pd.Series(qnorm(t.e, t.m), index=t.index)
    return out


def max_window(genes, pos, w=1_000_000):
    p = pos.reindex(genes).dropna()
    best = 0
    for c, x in p.groupby("chr"):
        s = np.sort(x.start.to_numpy())
        best = max(best, (np.searchsorted(s, s + w) - np.arange(len(s))).max())
    return best / max(len(genes), 1)


def main():
    pos = pd.read_csv(CT / "Loxafr3.0_gtf_lifted_to_LoxAfr3.0_HiC_capitalized.bed", sep="\t", header=None,
                      names=["chr", "start", "end", "gene"]).drop_duplicates("gene").set_index("gene")
    sets = {lib: gmt(lib) for lib in LIBS}
    univ = set().union(*(set().union(*s.values()) for s in sets.values()))
    vals = values(pos)
    for name, d in vals.items():
        if VALUES and name not in VALUES:
            continue
        d = d[d.index.isin(univ)].sort_values(ascending=False)
        print(f"{name}: {len(d)} genes ranked; median dCRUSH {d.median():.1f}")
        allsig = []
        for lib in LIBS:
            res = gp.prerank(rnk=d, gene_sets=sets[lib], min_size=15, max_size=500, permutation_num=1000,
                             weight=1, seed=0, threads=8, outdir=None, verbose=False).res2d
            res = res.rename(columns={"Term": "term", "ES": "ES", "NES": "NES", "NOM p-val": "p",
                                      "FDR q-val": "FDR", "FWER p-val": "FWER", "Tag %": "tag_pct",
                                      "Gene %": "gene_pct", "Lead_genes": "leading_edge"})
            res[["NES", "p", "FDR", "FWER"]] = res[["NES", "p", "FDR", "FWER"]].astype(float)
            res["direction"] = np.where(res.NES > 0, "A in mammoth", "A in elephant")
            le = res.leading_edge.str.split(";")
            res["n_leading_edge"] = le.str.len()
            res["max_LE_1Mb"] = [max_window(x, pos) for x in le]
            res["clustered"] = np.where(res.max_LE_1Mb >= CLUST, "yes", "")
            res.insert(0, "library", lib)
            res = res.drop(columns=["Name"]).sort_values(["FDR", "p"])
            res.to_csv(HERE / f"gsea_{name}_{lib}.tsv", sep="\t", index=False, float_format="%.4g")
            s = res[res.FDR < 0.25]
            allsig.append(s)
            print(f"  {lib}: {len(res)} sets; FDR<0.25 {len(s)} (mammoth {(s.NES > 0).sum()}, "
                  f"elephant {(s.NES < 0).sum()}); FDR<0.05 {(res.FDR < 0.05).sum()}; "
                  f"clustered among FDR<0.25 {(s.clustered == 'yes').sum()}")
        a = pd.concat(allsig)
        a.to_csv(HERE / f"gsea_{TAG}{name}_all.tsv", sep="\t", index=False, float_format="%.4g")
        print(a[["library", "term", "NES", "p", "FDR", "n_leading_edge", "max_LE_1Mb"]]
              .head(40).to_string(index=False))

    r = pd.read_csv(HERE / f"gsea_{TAG}gene_all.tsv", sep="\t")
    r = r[r.clustered != "yes"].drop_duplicates("leading_edge")
    top = pd.concat([r[r.NES < 0].nsmallest(NTOP, "NES"), r[r.NES > 0].nlargest(NTOP, "NES")]).sort_values("NES")
    if top.empty:
        return
    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(3.0, 0.18 * len(top) + 0.6))
        y = np.arange(len(top))
        ax.barh(y, top.NES, color=np.where(top.NES > 0, M_COL, E_COL), height=0.7, lw=0)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_yticks(y)
        ax.tick_params(labelsize=6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if notext:
            ax.set_yticklabels([])
            ax.set_xticklabels([])
        else:
            ax.set_yticklabels([f"{t.split(' (GO')[0]} ({q:.2g})" for t, q in zip(top.term, top.FDR)], fontsize=5.5)
            for lab, lib in zip(ax.get_yticklabels(), top.library):
                lab.set_fontstyle("italic" if lib.startswith("MSigDB") else "normal")
            ax.set_xlabel("NES (< 0 A in elephant, > 0 A in mammoth)", fontsize=7)
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"gsea_{TAG}top{'_notext' if notext else ''}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
