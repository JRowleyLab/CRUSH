#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, "revision")
import MSD_tracks as msd

HIC = "full.hic"
TRACK = "full_Crush_1000.bw"
EIGEN = "e_full.bedgraph"
RES = 1000
SPAN = 10_000_000


def track(path: str, chrom: str):
    d = msd.load(path).dropna(subset=["Score"])
    c = msd.norm_chrom(d["Chromosome"]).to_numpy()
    d = d[c == chrom]
    n = int(d["Start"].max()) // RES + 1
    v = np.zeros(n)
    v[(d["Start"].to_numpy() // RES).astype(int)] = d["Score"].to_numpy()
    return v


def pick_window(v, span_bins):
    """Window with the most balanced A/B split among well-covered windows."""
    cov = (v != 0).astype(float)
    a = (v > 0).astype(float)
    k = np.ones(span_bins)
    csum = lambda x: np.convolve(x, k, "valid")
    c, ap = csum(cov), csum(a)
    ok = c > 0.9 * span_bins
    if not ok.any():
        ok = c > 0.5 * span_bins
    frac = np.where(c > 0, ap / np.maximum(c, 1), 0)
    score = np.where(ok, -np.abs(frac - 0.5), -np.inf)
    return int(np.argmax(score))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="20")
    ap.add_argument("--start", type=int, default=None, help="window start, bp")
    args = ap.parse_args()

    v = track(TRACK, args.chrom)
    ev = track(EIGEN, args.chrom)
    nb = SPAN // RES
    s = args.start // RES if args.start is not None else pick_window(v, nb)
    lo, hi = s, s + nb
    sub = v[lo:hi]
    esub = np.zeros(nb)
    m = min(hi, len(ev))
    if m > lo:
        esub[: m - lo] = ev[lo:m]
    print(f"chr{args.chrom}:{lo*RES:,}-{hi*RES:,}  ({nb:,} bins at {RES} bp)")
    print(f"  covered {100*np.mean(sub!=0):.1f}%   A fraction "
          f"{100*np.mean(sub[sub!=0]>0):.1f}%", flush=True)

    import hicstraw
    reg = f"{args.chrom}:{lo*RES}:{hi*RES}"
    recs = []
    for norm in ("KR", "SCALE", "VC", "NONE"):
        try:
            recs = hicstraw.straw("observed", norm, HIC, reg, reg, "BP", RES)
            if recs:
                print(f"  norm={norm}  {len(recs):,} records", flush=True)
                break
        except Exception as e:
            print(f"  {norm}: {str(e)[:60]}", flush=True)
    if not recs:
        raise SystemExit("no contacts")

    M = np.zeros((nb, nb), np.float32)
    x = np.fromiter((r.binX // RES - lo for r in recs), np.int32, len(recs))
    y = np.fromiter((r.binY // RES - lo for r in recs), np.int32, len(recs))
    c = np.fromiter((r.counts for r in recs), np.float32, len(recs))
    ok = (x >= 0) & (x < nb) & (y >= 0) & (y < nb) & np.isfinite(c)
    x, y, c = x[ok], y[ok], c[ok]
    M[x, y] = c
    M[y, x] = c

    d = np.abs(np.subtract.outer(np.arange(nb), np.arange(nb)))
    tot = np.bincount(d.ravel(), weights=M.ravel(), minlength=nb)
    cnt = np.bincount(d.ravel(), minlength=nb).astype(np.float64)
    exp = np.divide(tot, cnt, out=np.zeros_like(tot), where=cnt > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        OE = np.where(exp[d] > 0, M / exp[d], np.nan).astype(np.float32)

    keep = (sub != 0) & (esub != 0) & np.isfinite(OE).sum(0).astype(bool)
    OE = OE[np.ix_(keep, keep)]
    sc, esc = sub[keep], esub[keep]
    base = HERE / f"win_chr{args.chrom}_{lo*RES}"
    np.savez_compressed(str(base) + ".npz", OE=OE, crush=sc, eigen=esc,
                        chrom=args.chrom, start=lo * RES, res=RES)
    print(f"  kept {keep.sum():,} bins called by both methods", flush=True)

    def strength(scores, label):
        o = np.argsort(scores)
        M = OE[np.ix_(o, o)]
        n = len(scores); q = n // 4
        bb = np.nanmean(M[:q, :q]); aa = np.nanmean(M[-q:, -q:])
        ab = np.nanmean(M[:q, -q:])
        print(f"  {label:<12} A-A {aa:.3f}  B-B {bb:.3f}  A-B {ab:.3f}   "
              f"(AA+BB)/2AB = {(aa+bb)/(2*ab):.3f}")
        return (aa + bb) / (2 * ab)

    print()
    strength(sc, "CRUSH")
    strength(esc, "eigenvector")


if __name__ == "__main__":
    main()
