#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import iqr

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
sys.path.insert(0, str(FIG.parent))
import MSD_tracks as msd
REF = Path("hg19_genes_refseq_nochr.bed")
RECALL = Path("tnbc")
PUB = Path("Remapped_1kb")

SUBTYPE = [("Normal", ["hmec"]), ("Luminal", ["t47D", "zr7530"]),
           ("HER2+", ["hcc1954"]), ("TNBC", ["hcc70", "bt549"])]
KEYS = [k for _, ks in SUBTYPE for k in ks]

SOURCES = {
    "crush": lambda k: RECALL / f"{k}_r1000" / f"{k}_mergedCrush_1000.bedgraph",
    "eigen": lambda k: PUB / "Eigen_KR" / k / f"Eigen_{k}_100kb.bedgraph",
    "crush_published": lambda k: PUB / "crush" / f"{k}_combined"
                                 / "mergedCrush_1000.bedgraph",
}

_cache: dict = {}


def track(src: str, key: str):
    if (src, key) not in _cache:
        df = msd.load(str(SOURCES[src](key))).dropna(subset=["Score"])
        _cache[(src, key)] = (df, float(iqr(df["Score"].to_numpy())),
                              msd.norm_chrom(df["Chromosome"]).to_numpy())
    return _cache[(src, key)]


def locus(symbol: str, promoter):
    rows = [l.rstrip("\n").split("\t") for l in REF.open()]
    hit = [r for r in rows if r[5] == symbol]
    if not hit:
        raise SystemExit(f"{symbol}: not in {REF.name}")
    ch = hit[0][0]
    if len({r[0] for r in hit}) > 1:
        raise SystemExit(f"{symbol}: transcripts on multiple chromosomes")
    lo, hi = min(int(r[1]) for r in hit), max(int(r[2]) for r in hit)
    if promoter is None:
        return ch, lo, hi, "gene body"
    strand = hit[0][4]
    tss = hi if strand == "-" else lo
    return (ch, tss - promoter, tss + promoter,
            f"TSS+/-{promoter//1000}kb (strand {strand}, TSS {tss:,})")


def score(src: str, key: str, ch: str, lo: int, hi: int):
    df, gi, c = track(src, key)
    sub = df[(c == ch) & (df["End"].to_numpy() > lo) & (df["Start"].to_numpy() < hi)]
    return (None, 0) if sub.empty else (float(sub["Score"].mean()) / gi, len(sub))


def main() -> None:
    argv = sys.argv[1:]
    promoter, suffix = None, ""
    raw = "--raw" in argv
    if raw:
        argv.remove("--raw"); suffix = "_raw"
    if "--promoter" in argv:
        i = argv.index("--promoter")
        promoter = int(argv[i + 1]); suffix = "_promoter" + suffix
        del argv[i:i + 2]
    genes = argv or ["RSU1"]
    OUT.mkdir(parents=True, exist_ok=True)
    for g in genes:
        ch, lo, hi, how = locus(g, promoter)
        print(f"\n{'='*72}\n{g}  chr{ch}:{lo:,}-{hi:,}  ({(hi-lo)/1000:.1f} kb)  [{how}]")
        rec: dict = {"gene": g, "locus": f"chr{ch}:{lo:,}-{hi:,} hg19", "window": how}
        for src in ("crush", "eigen", "crush_published"):
            per = {k: score(src, k, ch, lo, hi) for k in KEYS}
            base = 0.0 if raw else per["hmec"][0]
            if base is None:
                print(f"  {src}: HMEC has no bin over the locus -- cannot centre")
                rec[src] = None
                continue
            d = {k: (None if v is None else v - base) for k, (v, _) in per.items()}
            sub = {lab: (None if any(d[k] is None for k in ks)
                         else float(np.mean([d[k] for k in ks])))
                   for lab, ks in SUBTYPE}
            rec[src] = sub
            rec[f"{src}_per_line"] = d
            rec[f"{src}_n_bins"] = {k: per[k][1] for k in KEYS}
            print(f"  {src:<16}" + "".join(f"{lab}={sub[lab]:+.3f}  "
                                          if sub[lab] is not None else f"{lab}=NA  "
                                          for lab, _ in SUBTYPE))
            print(f"  {'':<16}per line: " + "  ".join(
                f"{k}={d[k]:+.3f}" if d[k] is not None else f"{k}=NA" for k in KEYS))
            print(f"  {'':<16}bins:     " + "  ".join(
                f"{k}={per[k][1]}" for k in KEYS))
        p = OUT / f"fig1G_{g}{suffix}_subtype.json"
        p.write_text(json.dumps(rec, indent=1))
        print(f"  -> {p}")


if __name__ == "__main__":
    main()
