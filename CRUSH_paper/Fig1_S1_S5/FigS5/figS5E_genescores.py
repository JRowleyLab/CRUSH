#!/usr/bin/env python3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
REV = FIG.parent

RES = 1_000
WINDOW_BP = 10_000
RMIN = 0.05
REF = Path("hg19_genes_refseq_nochr.bed")

TNBC = ["hcc70", "bt549"]
RECEPTOR_POS = ["t47D", "zr7530", "hcc1954"]
NORMAL = ["hmec"]
LINES = NORMAL + RECEPTOR_POS + TNBC
CHROMS = [str(i) for i in range(1, 23)] + ["X"]


def refseq_tss():
    """One outermost TSS per symbol, single-locus genes only."""
    d = pd.read_csv(REF, sep="\t", header=None, usecols=[0, 1, 2, 4, 5],
                    names=["chrom", "start", "end", "strand", "gene"],
                    dtype={"chrom": str, "strand": str})
    d = d[d["chrom"].isin(CHROMS)]
    nun = d.groupby("gene")[["chrom", "strand"]].nunique()
    ok = nun[(nun["chrom"] == 1) & (nun["strand"] == 1)].index
    d = d[d["gene"].isin(ok)]
    plus = d[d["strand"] == "+"].groupby(["gene", "chrom"], as_index=False)["start"].min()
    plus = plus.rename(columns={"start": "tss"})
    minus = d[d["strand"] == "-"].groupby(["gene", "chrom"], as_index=False)["end"].max()
    minus = minus.rename(columns={"end": "tss"})
    t = pd.concat([plus, minus], ignore_index=True)
    t = t[~t["gene"].duplicated(keep=False)]
    return t.reset_index(drop=True)


def gene_density(t):
    b = t.assign(b=(t["tss"] // 100_000).astype(np.int64))
    return b.groupby(["chrom", "b"]).size().rename("n").reset_index()


def flip_map(bw_path, dens):
    """Chromosomes whose CRUSH sign disagrees with gene density."""
    bw = pyBigWig.open(str(bw_path))
    sizes = bw.chroms()
    flip = set()
    for c in CHROMS:
        k = c if c in sizes else "chr" + c
        if k not in sizes:
            continue
        n = sizes[k] // 100_000 + 1
        v = bw.stats(k, 0, min(n * 100_000, sizes[k]), type="mean", nBins=n)
        v = np.array([np.nan if x is None else x for x in v], float)
        sub = dens[dens["chrom"] == c].set_index("b")["n"]
        d = np.zeros(len(v))
        d[sub.index[sub.index < len(v)]] = sub[sub.index < len(v)].to_numpy()
        m = np.isfinite(v) & (v != 0)
        if m.sum() < 50 or d[m].std() == 0:
            continue
        r = np.corrcoef(v[m], d[m])[0, 1]
        if np.isfinite(r) and abs(r) >= RMIN and r < 0:
            flip.add(c)
    bw.close()
    return flip


def scores(line, t, dens):
    p = REV / "tracks" / "tnbc" / f"{line}_r1000" / f"{line}_mergedCrush_{RES}.bw"
    fl = flip_map(p, dens)
    bw = pyBigWig.open(str(p))
    sizes = bw.chroms()
    out = np.full(len(t), np.nan)
    for i, (c, tss) in enumerate(zip(t["chrom"], t["tss"])):
        k = c if c in sizes else "chr" + c
        if k not in sizes:
            continue
        lo, hi = max(0, int(tss) - WINDOW_BP), min(sizes[k], int(tss) + WINDOW_BP)
        if hi <= lo:
            continue
        try:
            q = bw.stats(k, lo, hi, type="mean")
        except RuntimeError:
            continue
        if q and q[0] is not None and np.isfinite(q[0]):
            out[i] = q[0] * (-1 if c in fl else 1)
    bw.close()
    print(f"    {line:<9}{int(np.isfinite(out).sum()):>7,} genes scored "
          f"({len(fl)} chroms flipped)", flush=True)
    return out


def main() -> None:
    t = refseq_tss()
    print(f"  {len(t):,} single-locus RefSeq symbols with an unambiguous TSS")
    dens = gene_density(t)
    sc = t.copy()
    for line in LINES:
        sc[line] = scores(line, t, dens)
    sc = sc.dropna(subset=LINES).reset_index(drop=True)
    print(f"\n  {len(sc):,} genes scored in all {len(LINES)} lines")

    OUT.mkdir(exist_ok=True)
    dest = OUT / "figS5E_genescores_tss.tsv"
    sc.to_csv(dest, sep="\t", index=False)
    print(f"  wrote {dest.name}")

    for tag, plus, minus in (("vs_rpos", TNBC, RECEPTOR_POS),
                             ("vs_hmec", TNBC, NORMAL),
                             ("bt549_vs_hcc70", ["bt549"], ["hcc70"])):
        delta = sc[plus].mean(axis=1) - sc[minus].mean(axis=1)
        r = pd.DataFrame({"gene": sc["gene"].str.upper(), "delta": delta})
        r = r.groupby("gene", as_index=False)["delta"].mean()
        rnk = OUT / f"figS5E_tnbc_{tag}.rnk"
        r.sort_values("delta", ascending=False).to_csv(rnk, sep="\t", header=False,
                                                       index=False)
        print(f"  {tag:<9} delta mean {delta.mean():+.2f}  sd {delta.std():.2f}"
              f"  -> {rnk.name}")

    go = sc[["gene", "chrom", "tss"]].copy()
    go["gene"] = go["gene"].str.upper()
    go.to_csv(OUT / "figS5E_geneorder.tsv", sep="\t", index=False)
    print(f"  wrote figS5E_geneorder.tsv")


if __name__ == "__main__":
    main()
