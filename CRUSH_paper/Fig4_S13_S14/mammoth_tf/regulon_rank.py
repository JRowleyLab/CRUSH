#!/usr/bin/env python3
from pathlib import Path

import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, spearmanr

HERE = Path(__file__).resolve().parent
REF = HERE / "ref"
CT = Path("comptracks")
S9 = HERE / "tableS9_published.tsv"
LIBS = {"ARCHS4_Coexpression": "coexp", "GTEx_Coexpression": "coexp", "Enrichr_Queries": "coocc",
        "ENCODE_ChIP-seq": "ChIP", "ReMap_ChIP-seq": "ChIP", "Literature_ChIP-seq": "ChIP"}
FIG4 = ["BHLHE40", "THAP5", "LHX2", "POU6F2", "FOXR2"]
MIN, NROT, SEED, NBIN = 15, 2000, 0, 20
ADJUST = "--noadjust" not in sys.argv
TAG = "" if ADJUST else "_noadjust"


def bh(p):
    p = np.asarray(p)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty_like(q)
    out[o] = np.minimum.accumulate(q[::-1])[::-1].clip(max=1)
    return out


def regulons(keep):
    reg = {}
    for lib in keep:
        for line in open(REF / f"{lib}.gmt"):
            f = line.rstrip("\n").split("\t")
            reg.setdefault(f[0].split("_")[0], set()).update(g for g in f[1:] if g)
    return reg


def main():
    g = pd.read_csv(CT / "genes_median", sep="\t", header=None, names=["gene", "e", "m"])
    g = g.drop_duplicates("gene").set_index("gene").dropna()
    ref = (np.sort(g.e.to_numpy()) + np.sort(g.m.to_numpy())) / 2
    for c in ("e", "m"):
        g[c + "q"] = ref[g[c].rank(method="first").astype(int) - 1]
    g["d"] = g["mq"] - g["eq"]

    pos = pd.read_csv(CT / "Loxafr3.0_gtf_lifted_to_LoxAfr3.0_HiC_capitalized.bed", sep="\t", header=None,
                      names=["chr", "start", "end", "gene"]).drop_duplicates("gene").set_index("gene")
    allreg = regulons(LIBS)
    univ = set().union(*allreg.values())
    u = g[g.index.isin(univ) & g.index.isin(pos.index)].join(pos).sort_values(["chr", "start"])
    n = len(u)
    idx = pd.Series(np.arange(n), index=u.index)
    print(f"universe {n} genes")

    rng = np.random.default_rng(SEED)
    shifts = rng.integers(n // 20, n - n // 20, NROT)
    out = []
    for name, keep in (("all", LIBS), ("ChIP", {k for k, t in LIBS.items() if t == "ChIP"})):
        reg = regulons(keep)
        deg = Counter(x for t, s in reg.items() for x in s - {t})
        d = u.d.copy()
        if ADJUST:
            b = pd.qcut(pd.Series([deg[x] for x in u.index], index=u.index), NBIN, duplicates="drop")
            d = d - d.groupby(b, observed=True).transform("median")
        rk = d.rank().to_numpy()
        med = d.median()
        tfs = [t for t in sorted(reg) if len((reg[t] - {t}) & set(idx.index)) >= MIN]
        rows, cols = [], []
        for i, t in enumerate(tfs):
            j = idx[sorted((reg[t] - {t}) & set(idx.index))].to_numpy()
            rows += [i] * len(j)
            cols += list(j)
        M = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(tfs), n))
        k = np.asarray(M.sum(1)).ravel()

        def rbis(rsum):
            U = rsum - k * (k + 1) / 2
            return 2 * U / (k * (n - k)) - 1

        obs = rbis(M @ rk)
        z = (M @ rk - k * (n + 1) / 2) / np.sqrt(k * (n - k) * (n + 1) / 12)
        ge = np.zeros(len(tfs))
        for c in np.array_split(shifts, 20):
            R = np.stack([np.roll(rk, s) for s in c], axis=1)
            null = rbis((M @ R).T).T
            ge += (np.abs(null) >= np.abs(obs)[:, None]).sum(1)
        up = np.asarray(M @ (d.to_numpy() > med)).ravel()
        r = pd.DataFrame({"TF": tfs, "regulon": name, "n_targets": k.astype(int),
                          "n_more_A_mammoth": up.astype(int), "n_more_A_elephant": (k - up).astype(int),
                          "rank_biserial": obs, "z": z, "p_wilcoxon": 2 * norm.sf(np.abs(z)),
                          "p_rotation": (ge + 1) / (NROT + 1)})
        r["FDR_rotation"] = bh(r.p_rotation)
        out.append(r)
    r = pd.concat(out)
    r["TF_dCRUSH_norm"] = r.TF.map(g.d)
    r["TF_dCRUSH_raw"] = r.TF.map(g.m - g.e)
    r["concordant"] = np.where(r.TF_dCRUSH_norm.isna(), "", np.where(
        np.sign(r.TF_dCRUSH_norm) == np.sign(r.rank_biserial), "yes", "no"))
    s9 = pd.read_csv(S9, sep="\t").set_index("TF")
    r["TableS9_set"] = r.TF.map(s9.set).fillna("")
    r["fig4"] = np.where(r.TF.isin(FIG4), "yes", "")
    r.to_csv(HERE / f"regulon_rank{TAG}.tsv", sep="\t", index=False, float_format="%.4g")

    cols = ["TF", "TableS9_set", "TF_dCRUSH_raw", "TF_dCRUSH_norm", "n_targets", "n_more_A_elephant",
            "n_more_A_mammoth", "rank_biserial", "p_wilcoxon", "p_rotation", "FDR_rotation", "concordant"]
    for name in ("all", "ChIP"):
        x = r[r.regulon == name]
        print(f"\n=== regulon: {name} ({len(x)} TFs) ===")
        print(x[x.fig4 == "yes"][cols].to_string(index=False))
        sig = x[x.FDR_rotation < 0.05]
        print(f"FDR<0.05: {len(sig)} TFs ({(sig.rank_biserial > 0).sum()} toward mammoth A, "
              f"{(sig.rank_biserial < 0).sum()} toward elephant A); p_wilcoxon<0.05: {(x.p_wilcoxon < 0.05).sum()}")
        xx = x.dropna(subset=["TF_dCRUSH_norm"])
        rho, p = spearmanr(xx.TF_dCRUSH_norm, xx.rank_biserial)
        print(f"TF own dCRUSH vs regulon r: Spearman {rho:.3f} (p {p:.2g}, n {len(xx)}); "
              f"concordant {(xx.concordant == 'yes').mean():.3f} of all, "
              f"{(xx[xx.FDR_rotation < 0.05].concordant == 'yes').mean():.3f} of FDR<0.05")
        for st in ("elephant", "mammoth"):
            y = x[x.TableS9_set == st]
            print(f"Table S9 {st} TFs: {len(y)}; regulon r median {y.rank_biserial.median():.3f}; "
                  f"r>0 {(y.rank_biserial > 0).sum()}, r<0 {(y.rank_biserial < 0).sum()}; "
                  f"FDR<0.05 {(y.FDR_rotation < 0.05).sum()}")


if __name__ == "__main__":
    main()
