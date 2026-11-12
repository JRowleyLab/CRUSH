#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REV = HERE.parents[1]
sys.path.insert(0, str(REV))
import MSD_tracks as msd

V2 = HERE / "haplotype"
V1 = Path("GSE99991")
EIG = Path("GSE99991")
ARMS = ["WT_Xa", "WT_Xi", "KO_Xa", "KO_Xi"]
OUT = HERE / "results"


def load(p: Path) -> pd.DataFrame:
    t = msd.load(str(p)).copy()
    t["Chromosome"] = msd.norm_chrom(t["Chromosome"]).to_numpy()
    return t[t["Chromosome"] == "X"].sort_values("Start")


def best_split(pos: np.ndarray, val: np.ndarray, step: int = 500_000):
    """Scan every candidate boundary; return (position, R^2) of the best.

    One-way ANOVA with two groups -- the model a megadomain implies. Scanning
    beats assuming a coordinate: it finds the split if one exists, and reports
    a weak best-R^2 if none does.
    """
    lo, hi = pos.min() + 5_000_000, pos.max() - 5_000_000
    grand = val.mean()
    tot = float(np.sum((val - grand) ** 2))
    if tot == 0:
        return None, float("nan")
    bestp, bestr = None, -1.0
    for b in range(int(lo), int(hi), step):
        left = pos < b
        nl, nr = int(left.sum()), int((~left).sum())
        if nl < 100 or nr < 100:
            continue
        fit = (nl * (val[left].mean() - grand) ** 2
               + nr * (val[~left].mean() - grand) ** 2)
        r2 = fit / tot
        if r2 > bestr:
            bestr, bestp = r2, b
    return bestp, float(bestr)


def agree(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    ka = msd.keyed(a, literal=True)
    kb = msd.keyed(b, literal=True)
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        m = msd.align(ka, kb, True, "t")
    x, y = m["Score_x"].to_numpy(), m["Score_y"].to_numpy()
    r = float(np.corrcoef(x, y)[0, 1])
    flipped = r < 0
    ya = -y if flipped else y
    return {"n": int(len(m)), "r_raw": r, "r_aligned": abs(r),
            "sign_agreement_raw": float(np.mean(np.sign(x) == np.sign(y))),
            "sign_agreement_aligned": float(np.mean(np.sign(x) == np.sign(ya))),
            "flipped": bool(flipped)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    res: dict = {}

    v2 = {a: load(V2 / a / f"{a}_mergedCrush_5000.bedgraph") for a in ARMS}
    v1 = {}
    for a in ARMS:
        p = V1 / a / "mergedCrush_5000.bedgraph"
        if p.exists():
            v1[a] = load(p)

    print("=== 1. does CRUSH separate the two alleles? (5 kb, shared bins) ===",
          flush=True)
    print(f"{'comparison':<22}{'v0.2 r':>9}{'v0.2 sign':>11}"
          f"{'v1.0 r':>9}{'v1.0 sign':>11}", flush=True)
    pairs = [("WT_Xa", "WT_Xi"), ("KO_Xa", "KO_Xi"),
             ("WT_Xa", "KO_Xa"), ("WT_Xi", "KO_Xi")]
    res["allele_agreement"] = {}
    for a, b in pairs:
        g2 = agree(v2[a], v2[b])
        g1 = agree(v1[a], v1[b]) if a in v1 and b in v1 else None
        res["allele_agreement"][f"{a} vs {b}"] = {"v2": g2, "v1": g1}
        s1r = f"{g1['r_aligned']:>9.3f}" if g1 else f"{'-':>9}"
        s1s = f"{100*g1['sign_agreement_aligned']:>10.1f}%" if g1 else f"{'-':>11}"
        flag = " (sign-flipped)" if g2["flipped"] else ""
        print(f"{a+' vs '+b:<22}{g2['r_aligned']:>9.3f}"
              f"{100*g2['sign_agreement_aligned']:>10.1f}%{s1r}{s1s}{flag}",
              flush=True)

    print("\n=== 2. strongest two-block split on each allele ===", flush=True)
    print(f"{'track':<26}{'best boundary':>16}{'R^2':>9}", flush=True)
    res["best_split"] = {}
    for a in ARMS:
        t = v2[a]
        p, r2 = best_split(t["Start"].to_numpy(), t["Score"].to_numpy())
        res["best_split"][f"CRUSH {a}"] = {"pos": p, "r2": r2}
        print(f"{'CRUSH ' + a:<26}{p/1e6 if p else float('nan'):>13.1f} Mb{r2:>9.4f}",
              flush=True)
    for a, f in [("WT_Xa", "WT1.merged.cas.eigen_500kb.bedgraph"),
                 ("WT_Xi", "WT1.merged.mus.eigen_500kb.bedgraph"),
                 ("KO_Xa", "KO1.merged.cas.eigen_500kb.bedgraph"),
                 ("KO_Xi", "KO1.merged.mus.eigen_500kb.bedgraph")]:
        p_ = EIG / f
        if not p_.exists():
            continue
        t = load(p_)
        if t.empty:
            continue
        pos, r2 = best_split(t["Start"].to_numpy(), t["Score"].to_numpy())
        res["best_split"][f"eigenvector {a}"] = {"pos": pos, "r2": r2}
        print(f"{'eigenvector ' + a:<26}{pos/1e6 if pos else float('nan'):>13.1f} Mb"
              f"{r2:>9.4f}", flush=True)

    print("\n=== 3. Xist-KO control: is Xi drawn toward Xa when Xist is lost? ===",
          flush=True)
    wt = res["allele_agreement"]["WT_Xa vs WT_Xi"]["v2"]
    ko = res["allele_agreement"]["KO_Xa vs KO_Xi"]["v2"]
    print(f"  WT  Xa vs Xi : r={wt['r_aligned']:.3f}  "
          f"sign={100*wt['sign_agreement_aligned']:.1f}%", flush=True)
    print(f"  KO  Xa vs Xi : r={ko['r_aligned']:.3f}  "
          f"sign={100*ko['sign_agreement_aligned']:.1f}%", flush=True)
    d = ko["r_aligned"] - wt["r_aligned"]
    print(f"  -> Xist loss moves the two alleles {'together' if d > 0 else 'apart'}"
          f" by dr = {d:+.3f}", flush=True)

    (OUT / "haplotype_x.json").write_text(json.dumps(res, indent=1))
    print(f"\nwrote {OUT / 'haplotype_x.json'}")


if __name__ == "__main__":
    main()
