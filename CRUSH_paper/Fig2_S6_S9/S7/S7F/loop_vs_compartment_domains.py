#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from matplotlib_venn import venn2

import domain_classes as dc

HERE = Path(__file__).resolve().parent
INS = "highResInsulation_Total.bw"
USE = ["CRUSH v2 1 kb", "eigenvector 1 kb", "eigenvector 25 kb"]
MIN_COMP = 10
RECIP = 0.8
NB, FL = 50, 25


def comp_domains(a: np.ndarray) -> list[tuple[int, int, int]]:
    """(start_bin, end_bin, sign) runs of the smoothed sign, NaN gaps bridged."""
    sg = dc.sign_track(a)
    idx = np.flatnonzero(np.isfinite(sg) & (sg != 0))
    if idx.size == 0:
        return []
    v = sg[idx]
    cut = np.flatnonzero(v[1:] != v[:-1]) + 1
    out = []
    for lo, hi in zip(np.r_[0, cut], np.r_[cut, len(idx)]):
        s, e = idx[lo], idx[hi - 1] + 1
        if e - s >= MIN_COMP:
            out.append((s, e, int(v[lo])))
    return out


def matches(loops: np.ndarray, comps: np.ndarray) -> np.ndarray:
    """Boolean per loop domain: any compartment domain with reciprocal overlap >= RECIP."""
    hit = np.zeros(len(loops), bool)
    if len(comps) == 0:
        return hit
    order = np.argsort(comps[:, 0])
    cs, ce = comps[order, 0], comps[order, 1]
    for i, (s, e) in enumerate(loops):
        lo = np.searchsorted(cs, s - (e - s) * (1 / RECIP))
        hi = np.searchsorted(cs, e)
        for j in range(lo, hi):
            ov = min(e, ce[j]) - max(s, cs[j])
            if ov > 0 and ov >= RECIP * (e - s) and ov >= RECIP * (ce[j] - cs[j]):
                hit[i] = True
                break
    return hit


def meta(ins: dict, doms: dict) -> np.ndarray:
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
            body = np.array_split(v[s:e], NB)
            lf = np.array_split(v[a:s], FL)
            rf = np.array_split(v[e:b], FL)
            rows.append([np.nanmean(x) if x.size else np.nan for x in lf + body + rf])
    return np.nanmean(np.array(rows), axis=0) if rows else np.full(NB + 2 * FL, np.nan)


def main():
    sizes = {l.split()[0]: int(l.split()[1]) for l in open(dc.SIZES)}
    lp = pd.read_csv(dc.LOOPS, sep="\t", usecols=["chr1", "x1", "x2", "y1", "y2"], dtype={"chr1": str})
    lp["chr1"] = lp["chr1"].str.replace("chr", "", regex=False)
    lp = lp[lp.chr1.isin(dc.CHROMS)]
    loops = {c: np.c_[((g.x1 + g.x2) // 2 // dc.RES).to_numpy(), ((g.y1 + g.y2) // 2 // dc.RES).to_numpy()]
             for c, g in lp.groupby("chr1")}
    n_loop = int(sum(len(v) for v in loops.values()))

    bw = pyBigWig.open(INS)
    ins = {}
    for c in dc.CHROMS:
        name = next((k for k in (c, "chr" + c) if k in bw.chroms()), None)
        if name:
            n = bw.chroms()[name] // dc.RES
            v = np.array([np.nan if x is None else x for x in bw.stats(name, 0, n * dc.RES, nBins=n)])
            v[v <= 0] = np.nan
            ins[c] = np.log2(v / np.nanmean(v))
    bw.close()

    res = {"n_loop_domains": n_loop, "loop_domain_median_kb": float(np.median(np.concatenate(
        [v[:, 1] - v[:, 0] for v in loops.values()]))), "tracks": {}}
    meta_curves = {"loop domains": meta(ins, loops)}
    fig_e, ax_e = plt.subplots(1, len(USE), figsize=(11, 3.6))
    for k, name in enumerate(USE):
        kind, path = dc.TRACKS[name]
        tr = dc.load_track(kind, path, sizes)
        comps = {c: np.array([(s, e, sg) for s, e, sg in comp_domains(a)]) for c, a in tr.items()}
        comps = {c: v for c, v in comps.items() if len(v)}
        n_comp = int(sum(len(v) for v in comps.values()))
        matched_loops = 0
        matched_comps = set()
        for c, L in loops.items():
            C = comps.get(c, np.empty((0, 3), int))
            hit = matches(L, C[:, :2]) if len(C) else np.zeros(len(L), bool)
            matched_loops += int(hit.sum())
            if len(C):
                for j, (s, e, _) in enumerate(C):
                    for (ls, le) in L[(L[:, 0] < e) & (L[:, 1] > s)]:
                        ov = min(e, le) - max(s, ls)
                        if ov >= RECIP * (le - ls) and ov >= RECIP * (e - s):
                            matched_comps.add((c, j))
                            break
        sizes_kb = np.concatenate([v[:, 1] - v[:, 0] for v in comps.values()])
        res["tracks"][name] = {"n_compartment_domains": n_comp,
                               "compartment_domain_median_kb": float(np.median(sizes_kb)),
                               "loop_domains_matching_a_compartment_domain": matched_loops,
                               "compartment_domains_matching_a_loop_domain": len(matched_comps),
                               "pct_loop_domains_matched": 100 * matched_loops / n_loop}
        print(name, res["tracks"][name], flush=True)
        both = len(matched_comps)
        v = venn2(subsets=(n_loop - matched_loops, n_comp - both, both),
                  set_labels=("loop domains", "compartment\ndomains"), ax=ax_e[k],
                  set_colors=("#7570b3", "#d95f02"), alpha=0.55)
        for t in v.set_labels:
            if t:
                t.set_fontsize(8)
        for t in v.subset_labels:
            if t:
                t.set_fontsize(7)
        ax_e[k].set_title(f"{name}\n(matched: reciprocal overlap ≥ {int(RECIP * 100)}%)", fontsize=8)
        A = {c: v[v[:, 2] > 0, :2] for c, v in comps.items()}
        B = {c: v[v[:, 2] < 0, :2] for c, v in comps.items()}
        meta_curves[f"A compartment domains, {name}"] = meta(ins, A)
        meta_curves[f"B compartment domains, {name}"] = meta(ins, B)
    fig_e.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig_e.savefig(HERE / f"domain_euler.{ext}", dpi=300, bbox_inches="tight")

    x = np.arange(NB + 2 * FL)
    fig, ax = plt.subplots(1, len(USE), figsize=(11, 3.2), sharey=True)
    for k, name in enumerate(USE):
        ax[k].plot(x, meta_curves["loop domains"], color="#7570b3", lw=1.2, label="loop domains")
        ax[k].plot(x, meta_curves[f"A compartment domains, {name}"], color="#b2182b", lw=1.2, label="A compartment domains")
        ax[k].plot(x, meta_curves[f"B compartment domains, {name}"], color="#2166ac", lw=1.2, label="B compartment domains")
        for b in (FL, FL + NB):
            ax[k].axvline(b, color="grey", lw=0.5, ls="--")
        ax[k].set_xticks([0, FL, FL + NB, NB + 2 * FL])
        ax[k].set_xticklabels(["-0.5L", "start", "end", "+0.5L"], fontsize=7)
        ax[k].set_title(name, fontsize=8)
        ax[k].tick_params(labelsize=7)
        for s in ("top", "right"):
            ax[k].spines[s].set_visible(False)
    ax[0].set_ylabel("insulation score\nlog2(score / chromosome mean)", fontsize=8)
    ax[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"domain_insulation_log2.{ext}", dpi=300, bbox_inches="tight")
    res["meta_insulation"] = {k: [float(y) for y in v] for k, v in meta_curves.items()}
    (HERE / "loop_vs_compartment_domains_log2.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
