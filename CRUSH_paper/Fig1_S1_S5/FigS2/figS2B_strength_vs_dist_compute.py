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
CRUSH_BW = "full_Crush_1000.bw"
EIGEN_BG = "e_full.bedgraph"
RES = 1000
BLOCK = 10_000_000


def track(path, chrom, n):
    d = msd.load(path).dropna(subset=["Score"])
    c = msd.norm_chrom(d["Chromosome"]).to_numpy()
    d = d[c == chrom]
    v = np.zeros(n)
    b = (d["Start"].to_numpy() // RES).astype(int)
    ok = b < n
    v[b[ok]] = d["Score"].to_numpy()[ok]
    return v


def npairs(a, b, maxd):
    """Exact count of pairs (i, i+d) with a[i] and b[i+d], for every d, via FFT."""
    n = len(a)
    L = 1 << int(np.ceil(np.log2(2 * n)))
    fa = np.fft.rfft(a.astype(np.float64), L)
    fb = np.fft.rfft(b.astype(np.float64), L)
    cc = np.fft.irfft(np.conj(fa) * fb, L)[:maxd + 1]
    return np.maximum(cc, 0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="20")
    ap.add_argument("--maxsep", type=int, default=50_000_000)
    args = ap.parse_args()
    import hicstraw
    hf = hicstraw.HiCFile(HIC)
    length = {c.name: c.length for c in hf.getChromosomes()}.get(args.chrom)
    n = length // RES + 1
    maxd = min(args.maxsep, length) // RES

    cr = track(CRUSH_BW, args.chrom, n)
    ei = track(EIGEN_BG, args.chrom, n)
    cov = (cr != 0) & (ei != 0)
    print(f"chr{args.chrom}: {n:,} bins, {cov.sum():,} called by both, "
          f"max separation {maxd*RES/1e6:.0f} Mb", flush=True)

    sums = {k: np.zeros(maxd + 1) for k in ("all", "cr_AA", "cr_BB", "cr_AB",
                                            "ei_AA", "ei_BB", "ei_AB")}
    idx = {}
    for tag, sc in (("cr", cr), ("ei", ei)):
        v = sc[cov]
        lo, hi = np.quantile(v, 0.25), np.quantile(v, 0.75)
        top = np.zeros(n, bool); bot = np.zeros(n, bool)
        w = np.where(cov)[0]
        top[w[v >= hi]] = True; bot[w[v <= lo]] = True
        idx[tag] = (top, bot)
        print(f"  {tag}: {top.sum():,} top-quartile, {bot.sum():,} bottom", flush=True)

    nb = (n * RES + BLOCK - 1) // BLOCK
    for i in range(nb):
        for j in range(i, nb):
            if (j - i) * BLOCK > args.maxsep + BLOCK:
                continue
            r1 = f"{args.chrom}:{i*BLOCK}:{min((i+1)*BLOCK, length)}"
            r2 = f"{args.chrom}:{j*BLOCK}:{min((j+1)*BLOCK, length)}"
            try:
                recs = hicstraw.straw("observed", "KR", HIC, r1, r2, "BP", RES)
            except Exception:
                continue
            if not recs:
                continue
            x = np.fromiter((r.binX // RES for r in recs), np.int64, len(recs))
            y = np.fromiter((r.binY // RES for r in recs), np.int64, len(recs))
            c = np.fromiter((r.counts for r in recs), np.float64, len(recs))
            ok = np.isfinite(c) & cov[x] & cov[y]
            x, y, c = x[ok], y[ok], c[ok]
            d = np.abs(y - x)
            keep = d <= maxd
            x, y, c, d = x[keep], y[keep], c[keep], d[keep]
            np.add.at(sums["all"], d, c)
            for tag in ("cr", "ei"):
                top, bot = idx[tag]
                aa = top[x] & top[y]; bb = bot[x] & bot[y]
                ab = (top[x] & bot[y]) | (bot[x] & top[y])
                np.add.at(sums[f"{tag}_AA"], d[aa], c[aa])
                np.add.at(sums[f"{tag}_BB"], d[bb], c[bb])
                np.add.at(sums[f"{tag}_AB"], d[ab], c[ab])
        print(f"    block row {i+1}/{nb}", flush=True)

    npr = {"all": npairs(cov, cov, maxd)}
    for tag in ("cr", "ei"):
        top, bot = idx[tag]
        npr[f"{tag}_AA"] = npairs(top, top, maxd)
        npr[f"{tag}_BB"] = npairs(bot, bot, maxd)
        npr[f"{tag}_AB"] = npairs(top, bot, maxd) + npairs(bot, top, maxd)

    with np.errstate(invalid="ignore", divide="ignore"):
        exp = np.where(npr["all"] > 0, sums["all"] / npr["all"], np.nan)
    out = {}
    for tag in ("cr", "ei"):
        for k in ("AA", "BB", "AB"):
            m = np.where(npr[f"{tag}_{k}"] > 0,
                         sums[f"{tag}_{k}"] / np.maximum(npr[f"{tag}_{k}"], 1), np.nan)
            out[f"{tag}_{k}"] = m / exp

    steps = np.arange(0, args.maxsep + 1, 1_000_000)
    rows = []
    print(f"\n  {'min sep':>9}{'CRUSH':>9}{'eigen':>9}{'ratio':>8}")
    for s in steps:
        b = s // RES
        if b > maxd - 10:
            break
        r = {}
        for tag in ("cr", "ei"):
            v = [np.nanmean(out[f"{tag}_{k}"][b:maxd + 1]) for k in ("AA", "BB", "AB")]
            r[tag] = (v[0] + v[1]) / (2 * v[2])
        rows.append((s, r["cr"], r["ei"]))
        if s % 5_000_000 == 0:
            print(f"  {s/1e6:>7.0f}M{r['cr']:>9.3f}{r['ei']:>9.3f}"
                  f"{r['cr']/r['ei']:>8.3f}")
    np.save(HERE / f"strength_vs_dist_chr{args.chrom}.npy", np.array(rows))
    print(f"\nwrote strength_vs_dist_chr{args.chrom}.npy ({len(rows)} points)")


if __name__ == "__main__":
    main()
