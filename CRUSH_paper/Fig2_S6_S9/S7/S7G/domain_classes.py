#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

HERE = Path(__file__).resolve().parent
RES = 1000
CHROMS = [str(i) for i in range(1, 23)]
SIZES = "hg19.sizes"
DOMAINS = "GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
LOOPS = "GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt"
HIC = "full.hic"
STRAW = "straw"
TRACKS = {
    "CRUSH v1 1 kb": ("bedgraph", "c_full.bedgraph"),
    "CRUSH v2 1 kb": ("bw", "full_mergedCrush_1000.bw"),
    "eigenvector 1 kb": ("bw", "EV_full_1kb.bw"),
    "eigenvector 25 kb": ("bw", "EV_full_25kb.bw"),
}
MIN_LEN, SMOOTH, BTOL, LTOL = 80_000, 5, 10, 25
INS_RES, INS_HALF, NNULL = 5000, 12, 20


def load_track(kind: str, path: str, sizes: dict) -> dict[str, np.ndarray]:
    out = {}
    if kind == "bedgraph":
        t = pd.read_csv(path, sep="\t", header=None, skiprows=1, names=["c", "s", "e", "v"], dtype={"c": str})
        t["c"] = t["c"].str.replace("chr", "", regex=False)
        for c, g in t.groupby("c"):
            if c in CHROMS:
                a = np.full(sizes[c] // RES + 1, np.nan)
                a[(g.s // RES).to_numpy()] = g.v.to_numpy()
                out[c] = a
    else:
        b = pyBigWig.open(path)
        for c in CHROMS:
            name = next((k for k in (c, "chr" + c) if k in b.chroms()), None)
            a = np.full(sizes[c] // RES + 1, np.nan)
            if name:
                for s, e, v in b.intervals(name) or []:
                    a[s // RES:e // RES] = v
            out[c] = a
        b.close()
    for c, a in out.items():
        a[a == 0] = np.nan
    return out


def sign_track(a: np.ndarray) -> np.ndarray:
    """Sign after a centred rolling median; NaN where no data."""
    s = pd.Series(a).rolling(SMOOTH, center=True, min_periods=1).median().to_numpy()
    return np.sign(s)


def switches(sg: np.ndarray) -> np.ndarray:
    """Bin indices i where the sign differs between the last called bin before
    i and bin i (gaps of NaN are bridged)."""
    idx = np.flatnonzero(np.isfinite(sg) & (sg != 0))
    v = sg[idx]
    return idx[1:][v[1:] != v[:-1]]


def classify(dom: pd.DataFrame, sw: dict, loops: dict) -> pd.DataFrame:
    rows = []
    for c, g in dom.groupby("c"):
        s_sw = np.sort(sw.get(c, np.array([], int)))
        lp = loops.get(c, np.empty((0, 2)))
        for s, e in zip(g.s.to_numpy(), g.e.to_numpy()):
            bs, be = s // RES, e // RES
            near = lambda b: np.any(np.abs(s_sw[np.searchsorted(s_sw, b - BTOL):np.searchsorted(s_sw, b + BTOL + 1)] - b) <= BTOL)
            comp = near(bs) and near(be)
            corner = bool(lp.size and np.any((np.abs(lp[:, 0] - bs) <= LTOL) & (np.abs(lp[:, 1] - be) <= LTOL)))
            inner = s_sw[np.searchsorted(s_sw, bs + BTOL + 1):np.searchsorted(s_sw, be - BTOL)]
            rows.append((c, s, e, comp, corner, inner.size))
    out = pd.DataFrame(rows, columns=["c", "s", "e", "comp", "corner", "n_internal"])
    out["cls"] = np.select([out.comp & out.corner, out.corner, out.comp], ["both", "loop", "compartment"], "neither")
    return out


def insulation(dom: pd.DataFrame, sizes: dict) -> pd.Series:
    """Per-domain mean insulation of its two boundaries (same-side / cross)."""
    def one_chrom(c):
        g = dom[dom.c == c]
        if g.empty:
            return {}
        cmd = (f"{STRAW} oe VC_SQRT {HIC} {c} {c} BP {INS_RES} | "
               f"awk '($2-$1) <= {2 * INS_HALF * INS_RES}'")
        txt = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True).stdout
        a = np.loadtxt(txt.splitlines(), ndmin=2) if txt else np.empty((0, 3))
        n = sizes[c] // INS_RES + 1
        band = {}
        for x, y, v in a:
            if np.isfinite(v):
                band[(int(x) // INS_RES, int(y) // INS_RES)] = v
        get = lambda i, j: band.get((min(i, j), max(i, j)), np.nan)
        res = {}
        for k, (s, e) in enumerate(zip(g.s, g.e)):
            vals = []
            for b in (s // INS_RES, e // INS_RES):
                if b - INS_HALF < 0 or b + INS_HALF >= n:
                    continue
                L = range(b - INS_HALF, b)
                R = range(b, b + INS_HALF)
                same = [get(i, j) for i in L for j in L if i < j] + [get(i, j) for i in R for j in R if i < j]
                cross = [get(i, j) for i in L for j in R]
                ms, mc = np.nanmean(same), np.nanmean(cross)
                if np.isfinite(ms) and np.isfinite(mc) and mc > 0:
                    vals.append(ms / mc)
            res[g.index[k]] = np.mean(vals) if vals else np.nan
        return res
    with ThreadPoolExecutor(6) as ex:
        parts = list(ex.map(one_chrom, CHROMS))
    out = {}
    for p in parts:
        out.update(p)
    return pd.Series(out)


def main():
    sizes = {l.split()[0]: int(l.split()[1]) for l in open(SIZES)}
    dom = pd.read_csv(DOMAINS, sep="\t", usecols=["chr1", "x1", "x2"], dtype={"chr1": str})
    dom.columns = ["c", "s", "e"]
    dom["c"] = dom["c"].str.replace("chr", "", regex=False)
    dom = dom[dom.c.isin(CHROMS) & (dom.e - dom.s > MIN_LEN)].reset_index(drop=True)
    lp = pd.read_csv(LOOPS, sep="\t", usecols=["chr1", "x1", "x2", "y1", "y2"], dtype={"chr1": str})
    lp["chr1"] = lp["chr1"].str.replace("chr", "", regex=False)
    loops = {c: np.c_[((g.x1 + g.x2) // 2 // RES).to_numpy(), ((g.y1 + g.y2) // 2 // RES).to_numpy()]
             for c, g in lp.groupby("chr1")}

    rng = np.random.default_rng(0)
    res = {"n_domains": int(len(dom)), "tracks": {}}
    classes = {}
    for name, (kind, path) in TRACKS.items():
        tr = load_track(kind, path, sizes)
        sw = {c: switches(sign_track(a)) for c, a in tr.items()}
        cl = classify(dom, sw, loops)
        classes[name] = cl
        null_comp = []
        for _ in range(NNULL):
            sh = dom.copy()
            for c, g in sh.groupby("c"):
                L = sizes[c]
                off = rng.integers(0, L, size=len(g))
                sh.loc[g.index, "s"] = (g.s + off) % L
                sh.loc[g.index, "e"] = sh.loc[g.index, "s"] + (g.e - g.s)
            sh = sh[sh.e < sh.c.map(sizes)]
            null_comp.append(classify(sh, sw, {}).comp.mean())
        loopdom = cl.corner
        res["tracks"][name] = {
            "switches_per_Mb": float(sum(v.size for v in sw.values()) / (sum(sizes[c] for c in CHROMS) / 1e6)),
            "class_counts": cl.cls.value_counts().to_dict(),
            "frac_both_boundaries_switch": float(cl.comp.mean()),
            "null_frac_both_boundaries_switch": float(np.mean(null_comp)),
            "loop_domains_with_internal_switch": float((cl.loc[loopdom, "n_internal"] > 0).mean()),
            "loop_domains_internal_switches_per_Mb": float(
                cl.loc[loopdom, "n_internal"].sum() / ((cl.loc[loopdom, "e"] - cl.loc[loopdom, "s"]).sum() / 1e6)),
        }
        print(name, json.dumps(res["tracks"][name]), flush=True)

    ins = insulation(dom, sizes)
    res["insulation_by_class"] = {}
    for name, cl in classes.items():
        cl = cl.assign(ins=ins.reindex(cl.index).to_numpy())
        res["insulation_by_class"][name] = cl.groupby("cls")["ins"].agg(["count", "median", "mean"]).round(4).to_dict(orient="index")
    (HERE / "domain_classes.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res["insulation_by_class"], indent=1))

    names = list(TRACKS)
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.3))
    x = np.arange(len(names))
    obs = [res["tracks"][n]["frac_both_boundaries_switch"] * 100 for n in names]
    nul = [res["tracks"][n]["null_frac_both_boundaries_switch"] * 100 for n in names]
    ax[0].bar(x - 0.2, obs, 0.4, color="#b2182b", label="observed")
    ax[0].bar(x + 0.2, nul, 0.4, color="#bbbbbb", label="shifted domains")
    ax[0].set_ylabel("domains with a switch at both\nboundaries (%)", fontsize=8)
    ax[0].legend(frameon=False, fontsize=7)
    ax[1].bar(x, [res["tracks"][n]["loop_domains_with_internal_switch"] * 100 for n in names],
              color=["#b2182b", "#ef8a62", "#2166ac", "#67a9cf"])
    ax[1].set_ylabel("corner-loop domains with an\ninternal compartment switch (%)", fontsize=8)
    for a in ax[:2]:
        a.set_xticks(x)
        a.set_xticklabels([n.replace(" ", "\n", 1) for n in names], fontsize=7)
    cl = classes["CRUSH v1 1 kb"].assign(ins=ins.reindex(classes["CRUSH v1 1 kb"].index).to_numpy())
    order = ["both", "loop", "compartment", "neither"]
    data = [cl.loc[cl.cls == k, "ins"].dropna().to_numpy() for k in order]
    ax[2].boxplot(data, showfliers=False, widths=0.6)
    ax[2].set_xticks(range(1, 5))
    ax[2].set_xticklabels([f"{k}\n(n={len(d):,})" for k, d in zip(order, data)], fontsize=7)
    ax[2].set_ylabel("boundary insulation\n(same side / cross, 5 kb O/E)", fontsize=8)
    ax[2].set_title("classes from CRUSH v1", fontsize=8)
    for a in ax:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.tick_params(labelsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"domain_classes.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
