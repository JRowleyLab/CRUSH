#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REV = HERE.parents[1]
sys.path.insert(0, str(REV))
import MSD_tracks as msd

CRUSH = HERE / "crush"
GENES = HERE / "refs" / "galGal5_genes.bed"
SADDLE = HERE / "results" / "saddle.json"
OUT = HERE / "results"
RES = 25_000

ORDER = ["WT-90G2-R1", "WT-90G2-R2", "WT-NS-R1", "WT-92Pphase-10m-R1",
         "WT-87PMphase-30m-R1", "WT-90PMphase-30m-R2", "WT-92PMphase-60m-R1"]
LABEL = {"WT-90G2-R1": "G2 (rep1)", "WT-90G2-R2": "G2 (rep2)",
         "WT-NS-R1": "non-sync", "WT-92Pphase-10m-R1": "prophase 10m",
         "WT-87PMphase-30m-R1": "prometa 30m (r1)",
         "WT-90PMphase-30m-R2": "prometa 30m (r2)",
         "WT-92PMphase-60m-R1": "prometa 60m"}
REFERENCE = "WT-90G2-R1"


def gene_density(res: int) -> pd.Series:
    """Genes per bin, keyed chrom_bin. TSS is strand-aware."""
    g = pd.read_csv(GENES, sep="\t", header=None, usecols=[0, 1, 2, 5],
                    names=["c", "s", "e", "st"], dtype={"c": str})
    tss = np.where(g["st"].to_numpy() == "-", g["e"].to_numpy(), g["s"].to_numpy())
    key = g["c"].to_numpy() + "_" + ((tss // res) * res).astype(str)
    return pd.Series(1, index=key).groupby(level=0).sum()


def track(sample: str) -> pd.Series:
    p = CRUSH / f"{sample}_r{RES}" / f"{sample}_mergedCrush_{RES}.bedgraph"
    if not p.exists():
        return None
    t = msd.load(str(p)).copy()
    t["Chromosome"] = msd.norm_chrom(t["Chromosome"]).to_numpy()
    k = t["Chromosome"] + "_" + t["Start"].astype(np.int64).astype(str)
    return pd.Series(t["Score"].to_numpy(), index=k.to_numpy())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dens = gene_density(RES)
    saddle = {}
    if SADDLE.exists():
        saddle = {k: v.get("mean_strength") for k, v in
                  json.loads(SADDLE.read_text()).items() if isinstance(v, dict)}

    tracks = {s: track(s) for s in ORDER}
    ref = tracks.get(REFERENCE)
    rows = []
    print(f"{'sample':<20}{'saddle':>9}{'r_gene':>10}{'r_G2':>9}{'n bins':>10}",
          flush=True)
    for s in ORDER:
        t = tracks.get(s)
        if t is None:
            print(f"{LABEL[s]:<20}   (not finished)"); continue
        d = dens.reindex(t.index).fillna(0.0)
        ok = np.isfinite(t.to_numpy()) & np.isfinite(d.to_numpy())
        r_gene = float(np.corrcoef(t.to_numpy()[ok], d.to_numpy()[ok])[0, 1])
        r_ref = float("nan")
        if ref is not None and s != REFERENCE:
            j = t.index.intersection(ref.index)
            if len(j) > 100:
                r_ref = float(np.corrcoef(t.reindex(j), ref.reindex(j))[0, 1])
        rows.append({"sample": s, "label": LABEL[s], "r_gene": r_gene,
                     "r_G2": r_ref, "saddle": saddle.get(s), "n": int(ok.sum())})
        sv = f"{saddle.get(s):.3f}" if saddle.get(s) else "-"
        rv = f"{r_ref:+.3f}" if np.isfinite(r_ref) else "  ref"
        print(f"{LABEL[s]:<20}{sv:>9}{r_gene:>+10.3f}{rv:>9}{ok.sum():>10,}",
              flush=True)

    if rows:
        inter = [r for r in rows if r["sample"].startswith(("WT-90G2", "WT-NS"))]
        mito = [r for r in rows if "PMphase" in r["sample"]]
        if inter and mito:
            gi = float(np.mean([r["r_gene"] for r in inter]))
            gm = float(np.mean([r["r_gene"] for r in mito]))
            print(f"\n  interphase mean r_gene   {gi:+.3f}")
            print(f"  prometaphase mean r_gene {gm:+.3f}")
            print(f"  retained: {100*gm/gi:.0f}% of the interphase correlation"
                  if gi else "")
    (OUT / "crush_vs_prior.json").write_text(json.dumps(rows, indent=1))
    print(f"\nwrote {OUT / 'crush_vs_prior.json'}")


if __name__ == "__main__":
    main()
