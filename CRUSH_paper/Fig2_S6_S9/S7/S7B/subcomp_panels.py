#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
FIGS = HERE.parent / "figures"
RES = 1000
SUBC = Path("GSE63525_GM12878_subcompartments.bed")
LOOPS = Path("GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt")
DOMAINS = Path("GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt")
V1 = Path("c_full.bedgraph")
V2 = Path("full_mergedCrush_1000.bw")
EIG = Path("EV_full_1kb.bw")
GRO = Path("GMgro.bigwig")
POL2 = Path("LCL_RNAPII_GSM935386.bw")
CHROMS = [str(i) for i in range(1, 23)]
STATES = ["A1", "A2", "B1", "B2", "B3"]
COL = {"A1": "#b2182b", "A2": "#ef8a62", "B1": "#67a9cf", "B2": "#2166ac", "B3": "#053061"}
FLANK = 250


def bw_bins(path: Path, chrom: str, n: int, stat: str = "mean") -> np.ndarray:
    b = pyBigWig.open(str(path))
    name = next((k for k in (chrom, "chr" + chrom) if k in b.chroms()), None)
    out = np.full(n, np.nan)
    if name is not None:
        m = min(n, b.chroms()[name] // RES)
        v = b.stats(name, 0, m * RES, type=stat, nBins=m)
        out[:m] = [np.nan if x is None else x for x in v]
    b.close()
    return out


def load_v1() -> dict[str, np.ndarray]:
    t = pd.read_csv(V1, sep="\t", header=None, skiprows=1, names=["c", "s", "e", "v"], dtype={"c": str})
    t["c"] = t["c"].str.replace("chr", "", regex=False)
    out = {}
    for c, g in t.groupby("c"):
        a = np.full(g["s"].max() // RES + 1, np.nan)
        a[(g["s"] // RES).to_numpy()] = g["v"].to_numpy()
        out[c] = a
    return out


def build() -> pd.DataFrame:
    sc = pd.read_csv(SUBC, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["c", "s", "e", "state"])
    sc["c"] = sc["c"].str.replace("chr", "", regex=False)
    v1 = load_v1()
    frames = []
    for c in CHROMS:
        n = len(v1.get(c, []))
        if n == 0:
            continue
        state = np.full(n, "", dtype=object)
        for s, e, st in sc.loc[sc["c"] == c, ["s", "e", "state"]].itertuples(index=False):
            state[s // RES:min(n, e // RES)] = st
        frames.append(pd.DataFrame({
            "c": c, "bin": np.arange(n), "state": state, "v1": v1[c],
            "v2": bw_bins(V2, c, n), "eig": bw_bins(EIG, c, n),
            "gro": bw_bins(GRO, c, n), "pol2": bw_bins(POL2, c, n),
        }))
    d = pd.concat(frames, ignore_index=True)
    d.loc[d["v1"] == 0, "v1"] = np.nan
    d.loc[d["v2"] == 0, "v2"] = np.nan
    return d


def profile(d: pd.DataFrame, col: str, sites: dict[str, np.ndarray]) -> np.ndarray:
    rows = []
    for c, g in d.groupby("c"):
        v = g[col].to_numpy()
        for p in sites.get(c, []):
            if FLANK <= p < len(v) - FLANK:
                rows.append(v[p - FLANK:p + FLANK + 1])
    return np.nanmean(np.array(rows), axis=0), len(rows)


def main() -> None:
    d = build()
    res: dict = {"track_versions": {"v1": str(V1), "v2": str(V2)}}

    res["1a"] = {}
    for ver in ("v1", "v2"):
        res["1a"][ver] = {}
        for st in STATES:
            x = d.loc[d["state"] == st, ver].dropna()
            res["1a"][ver][st] = {"n": int(len(x)), "median": float(x.median()), "pct_A": float((x > 0).mean() * 100)}

    res["1b"] = {}
    for ver in ("v1", "v2"):
        isB = d["state"].isin(["B1", "B2", "B3"])
        isA = d["state"].isin(["A1", "A2"])
        ok = d[ver].notna()
        groups = {
            "A in A subcomp": isA & ok & (d[ver] > 0),
            "A in B subcomp": isB & ok & (d[ver] > 0),
            "B in B subcomp": isB & ok & (d[ver] < 0),
        }
        out = {"frac_of_B_subcomp_bins_called_A": float((d.loc[isB & ok, ver] > 0).mean())}
        for k, m in groups.items():
            out[k] = {
                "n": int(m.sum()),
                "gro_mean": float(d.loc[m, "gro"].mean()),
                "gro_frac_nonzero": float((d.loc[m, "gro"] > 0).mean()),
                "pol2_mean": float(d.loc[m, "pol2"].mean()),
                "eig_calls_A": float((d.loc[m, "eig"] > 0).mean()),
            }
        runs = []
        for _, g in d[isB].groupby("c"):
            a = (g[ver].to_numpy() > 0).astype(int)
            edges = np.flatnonzero(np.diff(np.r_[0, a, 0]))
            runs += list(edges[1::2] - edges[::2])
        out["A_in_B_runs"] = {"n_runs": len(runs), "median_len_kb": float(np.median(runs)) if runs else 0.0}
        res["1b"][ver] = out

    lp = pd.read_csv(LOOPS, sep="\t", usecols=["chr1", "x1", "x2", "y1", "y2"], dtype={"chr1": str})
    lp["chr1"] = lp["chr1"].str.replace("chr", "", regex=False)
    anchors = {c: np.r_[((g.x1 + g.x2) // 2 // RES).to_numpy(), ((g.y1 + g.y2) // 2 // RES).to_numpy()]
               for c, g in lp.groupby("chr1")}
    dm = pd.read_csv(DOMAINS, sep="\t", usecols=["chr1", "x1", "x2"], dtype={"chr1": str})
    dm["chr1"] = dm["chr1"].str.replace("chr", "", regex=False)
    bounds = {c: np.unique(np.r_[(g.x1 // RES).to_numpy(), (g.x2 // RES).to_numpy()]) for c, g in dm.groupby("chr1")}
    res["2"] = {}
    prof = {}
    for ver in ("v1", "v2"):
        for name, sites in (("loop anchors", anchors), ("domain boundaries", bounds)):
            p, n = profile(d, ver, sites)
            prof[(ver, name)] = p
            res["2"][f"{ver} {name}"] = {"n": n, "centre": float(p[FLANK]),
                                        "flank_mean_200_250kb": float(np.nanmean(np.r_[p[:50], p[-50:]]))}
    (RESULTS / "subcomp_panels.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "track_versions"}, indent=1))

    fig, ax = plt.subplots(1, 4, figsize=(13, 3.2))
    data = [d.loc[d["state"] == st, "v1"].dropna().to_numpy() for st in STATES]
    bp = ax[0].boxplot(data, showfliers=False, patch_artist=True, widths=0.6, medianprops=dict(color="black"))
    for b, st in zip(bp["boxes"], STATES):
        b.set_facecolor(COL[st])
    ax[0].axhline(0, ls="--", lw=0.6, color="grey")
    ax[0].set_xticks(range(1, 6))
    ax[0].set_xticklabels(STATES)
    ax[0].set_ylabel("CRUSH 1 kb (v1)")
    ax[0].set_title("Rao 2014 subcompartments", fontsize=9)

    r = res["1b"]["v1"]
    keys = ["A in A subcomp", "A in B subcomp", "B in B subcomp"]
    cols = ["#b2182b", "#ef8a62", "#2166ac"]
    ax[1].bar(range(3), [r[k]["gro_frac_nonzero"] * 100 for k in keys], color=cols)
    ax[1].set_xticks(range(3))
    ax[1].set_xticklabels(["A in A", "A in B", "B in B"])
    ax[1].set_ylabel("bins with GRO-seq signal (%)")
    ax[1].set_title(f"A bins inside B subcomp.: {r['frac_of_B_subcomp_bins_called_A'] * 100:.1f}% of B", fontsize=9)
    ax[2].bar(range(3), [r[k]["eig_calls_A"] * 100 for k in keys], color=cols)
    ax[2].set_xticks(range(3))
    ax[2].set_xticklabels(["A in A", "A in B", "B in B"])
    ax[2].set_ylabel("called A by 1st eigenvector (%)")
    ax[2].set_title("Does the eigenvector see them?", fontsize=9)
    x = np.arange(-FLANK, FLANK + 1)
    for name, ls in (("loop anchors", "-"), ("domain boundaries", "--")):
        ax[3].plot(x, prof[("v1", name)], ls, color="#b2182b", label=f"{name} (n={res['2'][f'v1 {name}']['n']:,})")
    ax[3].axvline(0, color="grey", lw=0.5)
    ax[3].set_xlabel("distance (kb)")
    ax[3].set_ylabel("mean CRUSH 1 kb (v1)")
    ax[3].legend(frameon=False, fontsize=7)
    for a in ax:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.tick_params(labelsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIGS / f"subcomp_panels.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
