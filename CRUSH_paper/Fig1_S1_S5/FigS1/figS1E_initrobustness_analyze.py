#!/usr/bin/env python3
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
DATA = FIG / "init_robustness"
OUT = FIG / "results"
ARMS = ["std", "gconly", "rand", "rand2", "rand3", "rand4", "rand5"]


def load(arm: str):
    p = DATA / arm / f"{arm}_mergedCrush_1000.bedgraph"
    if not p.exists() or p.stat().st_size == 0:
        return None
    d = pd.read_csv(p, sep="\t", skiprows=1, header=None,
                    names=["c", "s", "e", "v"], dtype={"c": str})
    d = d[d["v"] != 0]
    return d.set_index(d["c"] + "_" + d["s"].astype(str))[["c", "v"]]


def compare(a, b):
    j = a.join(b, how="inner", lsuffix="_a", rsuffix="_b")
    sa, sb = np.sign(j["v_a"].to_numpy()), np.sign(j["v_b"].to_numpy())
    raw = float((sa == sb).mean())
    flipped, aligned = [], np.empty(len(j), bool)
    chrom = j["c_a"].to_numpy()
    for c in np.unique(chrom):
        m = chrom == c
        same = (sa[m] == sb[m]).mean()
        if same < 0.5:
            flipped.append(c)
            aligned[m] = sa[m] != sb[m]
        else:
            aligned[m] = sa[m] == sb[m]
    return {"n_bins": int(len(j)), "raw": raw, "aligned": float(aligned.mean()),
            "disagreeing_bins": int((~aligned).sum()),
            "chromosomes_flipped": sorted(flipped)}


def main() -> None:
    tracks = {a: load(a) for a in ARMS}
    have = [a for a in ARMS if tracks[a] is not None]
    print(f"arms present: {', '.join(have)}\n")
    res = {}
    print(f"{'comparison':<24}{'1 kb bins':>12}{'raw':>9}{'aligned':>10}"
          f"{'disagree':>11}{'chr flipped':>13}")
    for a, b in [("std", x) for x in have if x != "std"]:
        r = compare(tracks[a], tracks[b])
        res[f"{a} vs {b}"] = r
        print(f"{a+' vs '+b:<20}{r['n_bins']:>12,}{100*r['raw']:>8.1f}%"
              f"{100*r['aligned']:>9.3f}%{r['disagreeing_bins']:>11,}"
              f"{len(r['chromosomes_flipped']):>13}")
    rand = [a for a in have if a.startswith("rand")]
    vs_std = [res[f"std vs {a}"]["aligned"] for a in rand]
    dis = [res[f"std vs {a}"]["disagreeing_bins"] for a in rand]
    if vs_std:
        print(f"\n  {len(vs_std)} independent random initialisations vs standard:")
        print(f"    aligned agreement {100*min(vs_std):.3f}-{100*max(vs_std):.3f}%")
        print(f"    disagreeing 1 kb bins {min(dis):,}-{max(dis):,} "
              f"of {res['std vs '+rand[0]]['n_bins']:,}")
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "figS1E_initrobustness_agreement.json"
    dest.write_text(json.dumps(res, indent=1))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
