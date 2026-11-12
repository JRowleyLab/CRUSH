#!/usr/bin/env python3
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
S4C = Path(__file__).resolve().parent / "figS4C_rankrank.py"
src = S4C.read_text().split("def main(")[0]
m = type(sys)("m"); m.__file__ = str(S4C)
exec(compile(src, "m", "exec"), m.__dict__)

OUT = HERE.parent / "results"
GMT = OUT / "KEGG_2021_Human.gmt"
RES = 1000
MIN_TERM, MAX_TERM = 5, 500
TOPN = 12


WINDOW = sys.argv[1] if len(sys.argv) > 1 else "tss"

REFGTF = Path("galGal6.refGene.gtf")


def refgene_tss():
    """One outermost TSS per symbol, from refGene transcript lines."""
    import re
    rows = []
    pat = re.compile(r'gene_name "([^"]+)"')
    with REFGTF.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            g = pat.search(f[8])
            if not g:
                continue
            rows.append((g.group(1), f[0], int(f[3]) - 1, int(f[4]), f[6]))
    d = pd.DataFrame(rows, columns=["gene", "chrom", "start", "end", "strand"])
    nun = d.groupby("gene")[["chrom", "strand"]].nunique()
    ok = nun[(nun["chrom"] == 1) & (nun["strand"] == 1)].index
    d = d[d["gene"].isin(ok)]
    plus = d[d["strand"] == "+"].groupby(["gene", "chrom"], as_index=False)["start"].min()
    plus = plus.rename(columns={"start": "tss"})
    minus = d[d["strand"] == "-"].groupby(["gene", "chrom"], as_index=False)["end"].max()
    minus = minus.rename(columns={"end": "tss"})
    t = pd.concat([plus, minus], ignore_index=True)
    t = t[~t["gene"].duplicated(keep=False)]
    span = d.groupby(["gene", "chrom"], as_index=False).agg(lo=("start", "min"),
                                                            hi=("end", "max"))
    return t.merge(span, on=["gene", "chrom"]).reset_index(drop=True)


def gene_scores():
    import pyBigWig
    dens = m.gene_density()
    t = refgene_tss()
    print(f"    {len(t):,} refGene symbols with a single locus")
    out = {}
    for cell in ("Fibro", "Erythro", "eRBC", "aRBC"):
        bw_p = m.arm_bw("CRUSH 1 kb", cell)
        fl = m.flip_map(bw_p, dens)
        bw = pyBigWig.open(str(bw_p)); sizes = bw.chroms()
        vals = np.full(len(t), np.nan)
        for i, (c, tss, glo, ghi) in enumerate(zip(t["chrom"], t["tss"],
                                                   t["lo"], t["hi"])):
            if c not in sizes:
                continue
            if WINDOW == "body":
                lo, hi = max(0, int(glo)), min(sizes[c], int(ghi))
            else:
                lo, hi = max(0, tss - m.WINDOW_BP), min(sizes[c], tss + m.WINDOW_BP)
            if hi <= lo:
                continue
            try:
                q = bw.stats(c, int(lo), int(hi), type="mean")
            except RuntimeError:
                continue
            if q and q[0] is not None and np.isfinite(q[0]):
                vals[i] = q[0] * (-1 if c.removeprefix("chr") in fl else 1)
        bw.close()
        out[cell] = pd.Series(vals, index=t["gene"])
        print(f"    {cell:<9} {int(np.isfinite(vals).sum()):>7,} genes scored "
              f"({len(fl)} chroms oriented)", flush=True)
    return pd.DataFrame(out).dropna()


if __name__ == "__main__":
    import sys as _sys
    for w in ("tss", "body"):
        WINDOW = w
        globals()["WINDOW"] = w
        sc = gene_scores()
        dest = OUT / f"figS4_genescores_{w}.tsv"
        sc.to_csv(dest, sep="\t")
        print(f"  {w}: {len(sc):,} genes scored in all four lines -> {dest.name}")
