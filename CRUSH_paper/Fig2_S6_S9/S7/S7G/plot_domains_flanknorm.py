#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

import domain_classes as dc
import loop_vs_compartment_domains as lv

HERE = Path(__file__).resolve().parent
NB, FL = lv.NB, lv.FL


def meta_flank(ins, doms):
    rows = []
    for c, arr in doms.items():
        v = ins.get(c)
        if v is None:
            continue
        for s, e in arr:
            L = e - s
            a, b = s - L // 2, e + L // 2
            if a < 0 or b >= len(v):
                continue
            parts = np.array_split(v[a:s], FL) + np.array_split(v[s:e], NB) + np.array_split(v[e:b], FL)
            row = np.array([np.nanmean(x) if x.size else np.nan for x in parts])
            flank = np.nanmean(np.r_[row[:FL], row[FL + NB:]])
            if np.isfinite(flank):
                rows.append(row - flank)
    rows = np.array(rows)
    return np.nanmean(rows, axis=0), len(rows)


def main():
    sizes = {l.split()[0]: int(l.split()[1]) for l in open(dc.SIZES)}
    bw = pyBigWig.open(lv.INS)
    ins = {}
    for c in dc.CHROMS:
        n = bw.chroms()[c] // dc.RES
        v = np.array([np.nan if x is None else x for x in bw.stats(c, 0, n * dc.RES, nBins=n)])
        v[v <= 0] = np.nan
        ins[c] = np.log2(v / np.nanmean(v))
    bw.close()
    lp = pd.read_csv(dc.LOOPS, sep="\t", usecols=["chr1", "x1", "x2", "y1", "y2"], dtype={"chr1": str})
    lp["chr1"] = lp["chr1"].str.replace("chr", "", regex=False)
    lp = lp[lp.chr1.isin(dc.CHROMS)]
    loops = {c: np.c_[((g.x1 + g.x2) // 2 // dc.RES).to_numpy(), ((g.y1 + g.y2) // 2 // dc.RES).to_numpy()]
             for c, g in lp.groupby("chr1")}

    curves, ns = {}, {}
    curves["loop domains"], ns["loop domains"] = meta_flank(ins, loops)
    for name in ["CRUSH v2 1 kb", "eigenvector 1 kb", "eigenvector 25 kb"]:
        tr = dc.load_track(*dc.TRACKS[name], sizes)
        A = {}
        for c, a in tr.items():
            d = [(s, e) for s, e, sg in lv.comp_domains(a) if sg > 0]
            if d:
                A[c] = np.array(d)
        key = f"A compartment domains, {name}"
        curves[key], ns[key] = meta_flank(ins, A)
    res = {k: {"n": ns[k], "boundary_dip_start": float(curves[k][FL]), "boundary_dip_end": float(curves[k][FL + NB - 1]),
               "body_mean": float(np.nanmean(curves[k][FL:FL + NB])), "curve": [float(y) for y in curves[k]]}
           for k in curves}
    (HERE / "domain_insulation_flanknorm.json").write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        print(f"{k:45} n={v['n']:>6,}  start {v['boundary_dip_start']:+.3f}  end {v['boundary_dip_end']:+.3f}  body {v['body_mean']:+.3f}")

    lines = [("loop domains", "#7570b3", "-"),
             ("A compartment domains, CRUSH v2 1 kb", "#b2182b", "-"),
             ("A compartment domains, eigenvector 1 kb", "#2166ac", "--"),
             ("A compartment domains, eigenvector 25 kb", "#67a9cf", "--")]
    x = np.arange(NB + 2 * FL)
    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(4.6, 3.3))
        for key, col, ls in lines:
            ax.plot(x, curves[key], color=col, ls=ls, lw=1.3, label=key.replace("CRUSH v2 1 kb", "CRUSH 1 kb"))
        for b in (FL, FL + NB):
            ax.axvline(b, color="grey", lw=0.5, ls=":")
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xticks([0, FL, FL + NB, NB + 2 * FL])
        ax.tick_params(labelsize=7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if notext:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(["-0.5L", "start", "end", "+0.5L"], fontsize=7)
            ax.set_ylabel("insulation relative to flanks\n(log2, domain flank mean = 0)", fontsize=8)
            ax.legend(frameon=False, fontsize=6.5, loc="upper left")
        fig.tight_layout()
        stem = "domain_insulation_flanknorm" + ("_notext" if notext else "")
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"{stem}.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
