#!/usr/bin/env python3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results"
sys.path.insert(0, str(HERE))
NROT = 1000
TOPN = 8

RANKS, GSEA, LIB, DEST = sys.argv[1:5]
ORDER = sys.argv[5] if len(sys.argv) > 5 else None


def load_gmt(p):
    out = {}
    for line in Path(p).read_text().rstrip("\n").split("\n"):
        f = line.split("\t")
        if len(f) > 2:
            out[f[0]] = {g.split(",")[0].strip().upper() for g in f[2:] if g.strip()}
    return out


def es(order_idx, member, weights):
    """Weighted KS running-sum enrichment score for one ranking."""
    hit = member[order_idx]
    w = np.abs(weights[order_idx]) * hit
    tot = w.sum()
    if tot <= 0:
        return 0.0
    nmiss = (~hit).sum()
    run = np.cumsum(np.where(hit, w / tot, -1.0 / nmiss))
    return run[np.argmax(np.abs(run))]


def main():
    if ORDER:
        t = pd.read_csv(ORDER, sep="\t", dtype={"chrom": str})
    else:
        sys.argv = ["figS4_genescores.py", "body"]
        import figS4_genescores as mod
        t = mod.refgene_tss()
    t = t.assign(gene=t["gene"].str.upper())[["gene", "chrom", "tss"]]
    terms = load_gmt(LIB)

    for tag in ("only",):
        r = pd.read_csv(RANKS, sep="\t", header=None, names=["gene", "delta"])
        d = r.merge(t, on="gene").dropna(subset=["delta"])
        d = d.sort_values(["chrom", "tss"]).reset_index(drop=True)
        genes = d["gene"].to_numpy()
        vals = d["delta"].to_numpy(float)
        n = len(vals)
        obs_order = np.argsort(-vals)
        gsea = pd.read_csv(GSEA, sep="\t")
        gsea["abs"] = gsea["NES"].abs()
       
        sel = pd.concat([gsea[gsea.NES > 0].nlargest(TOPN, "NES"),
                         gsea[gsea.NES < 0].nsmallest(TOPN, "NES")])
        rng = np.random.default_rng(0)
        shifts = rng.integers(n // 20, n - n // 20, size=NROT)
        print(f"  {n:,} genes in genomic order, {NROT} rotations")
        print(f"  {'term':<48}{'NES':>8}{'n':>5}{'ES':>8}{'p(rot)':>9}")
        rows = []
        for _, g in sel.iterrows():
            member = np.isin(genes, list(terms.get(g["pathway"], ())))
            if member.sum() < 5:
                continue
            obs = es(obs_order, member, vals)
            null = np.empty(NROT)
            for k, s in enumerate(shifts):
                v = np.roll(vals, int(s))
                null[k] = es(np.argsort(-v), member, v)
            p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (NROT + 1)
            rows.append({"term": g["pathway"], "NES": g["NES"], "ES": obs, "p_rotation": p,
                         "n_genes": int(member.sum())})
            print(f"  {g['pathway'].split(' (GO:')[0][:48]:<48}"
                  f"{g['NES']:>8.2f}{int(member.sum()):>5}{obs:>8.3f}{p:>9.3f}")
        pd.DataFrame(rows).to_csv(DEST, sep="\t", index=False)
        print(f"  wrote {Path(DEST).name}")


if __name__ == "__main__":
    main()
