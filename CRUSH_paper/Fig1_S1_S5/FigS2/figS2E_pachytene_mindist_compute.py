#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DS = [100_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000]


def load(p: Path):
    if not p.exists() or p.stat().st_size == 0:
        return None
    d = pd.read_csv(p, sep="\t", skiprows=1, header=None,
                    names=["c", "s", "e", "v"], dtype={"c": str})
    d = d[d["v"] != 0]
    return d.set_index(d["c"] + "_" + d["s"].astype(str))[["c", "v"]]


def agree(a, b):
    j = a.join(b, how="inner", lsuffix="_a", rsuffix="_b")
    sa, sb = np.sign(j["v_a"].to_numpy()), np.sign(j["v_b"].to_numpy())
    chrom = j["c_a"].to_numpy()
    aligned = np.empty(len(j), bool)
    flipped = 0
    for c in np.unique(chrom):
        m = chrom == c
        same = (sa[m] == sb[m]).mean()
        aligned[m] = (sa[m] != sb[m]) if same < 0.5 else (sa[m] == sb[m])
        flipped += same < 0.5
    return len(j), float((sa == sb).mean()), float(aligned.mean()), int(flipped)


def main() -> None:
    out = {}
    for S, label in (("monkey_pac_250bp_intraonly", "pachytene"),
                     ("monkey_fibro", "fibroblast (somatic)")):
        ref = load(HERE / S / f"{S}_mergedCrush_1000.bedgraph")
        if ref is None:
            print(f"{label}: reference missing"); continue
        print(f"\n=== {label} ===")
        print(f"  {'-d':>8}{'bins':>12}{'raw':>9}{'aligned':>10}{'disagree':>12}")
        out[label] = {}
        for d in DS:
            v = load(HERE / f"{S}_d{d}" / "d_mergedCrush_1000.bedgraph")
            if v is None:
                print(f"  {d/1e6:>6.1f}M   (pending)"); continue
            n, raw, al, fl = agree(ref, v)
            out[label][str(d)] = {"n": n, "raw": raw, "aligned": al,
                                  "disagree": int(round((1 - al) * n))}
            print(f"  {d/1e6:>6.1f}M{n:>12,}{100*raw:>8.1f}%{100*al:>9.2f}%"
                  f"{int(round((1-al)*n)):>12,}")
    (HERE / "mindist_agreement.json").write_text(json.dumps(out, indent=1))
    print("\nreference for comparison -- somatic LCL (chr20, same test):")
    print("   1.0M 99.07%   2.0M 98.11%   5.0M 92.64%   10.0M 94.07%")


if __name__ == "__main__":
    main()
