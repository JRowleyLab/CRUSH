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
WORK = HERE / "work"
LCL = Path("LCL")
EIG = LCL / "EIGENtracks" / "e_full.bedgraph"
CR = LCL / "CRUSHtracks" / "full_mergedCrush_1000.bw"
HIC = LCL / "hiclinks" / "full.hic"
GRO = Path("GMgro.bigwig")
STRAW = "straw"
RES = 1000
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
CLASSES = ["B->B", "B->A", "A->B", "A->A"]
COLORS = {"B->B": "#3b8fd9", "B->A": "#d62728", "A->B": "#2ca02c", "A->A": "#ff7f0e"}
NJOBS = 8

AWK = r'''
FNR==NR { s[$1]=$2; next }
{
  x=$1; y=$2; n=$3
  if (y in s) { if (s[y]=="A") a[x]+=n; else b[x]+=n }
  if (x!=y && (x in s)) { if (s[x]=="A") a[y]+=n; else b[y]+=n }
}
END { for (k in s) print k "\t" (a[k]+0) "\t" (b[k]+0) }
'''


def bw_chrom(bw, c: str) -> str | None:
    return next((k for k in (c, "chr" + c) if k in bw.chroms()), None)


def load_calls() -> pd.DataFrame:
    e = pd.read_csv(EIG, sep="\t", header=None, usecols=[0, 1, 3],
                    names=["c", "s", "eig"], dtype={"c": str})
    e["c"] = e["c"].str.replace("chr", "", regex=False)
    bw = pyBigWig.open(str(CR))
    parts = []
    for c in CHROMS:
        name = bw_chrom(bw, c)
        if name is None:
            continue
        a = np.array(bw.intervals(name))
        parts.append(pd.DataFrame({"c": c, "s": (a[:, 0] // RES * RES).astype(np.int64),
                                   "cr": a[:, 2]}))
    bw.close()
    cr = pd.concat(parts)
    d = e[e["c"].isin(CHROMS)].merge(cr, on=["c", "s"])
    d = d[(d["eig"] != 0) & (d["cr"] != 0)].reset_index(drop=True)
    d["cls"] = (np.where(d["eig"] > 0, "A", "B").astype(object) + "->" +
                np.where(d["cr"] > 0, "A", "B").astype(object))
    return d, cr[cr["cr"] != 0]


def gro(d: pd.DataFrame) -> np.ndarray:
    bw = pyBigWig.open(str(GRO))
    out = np.full(len(d), np.nan)
    for c, g in d.groupby("c"):
        name = bw_chrom(bw, c)
        if name is None:
            continue
        nb = bw.chroms()[name] // RES
        v = np.array([x if x is not None else 0.0
                      for x in bw.stats(name, 0, nb * RES, type="mean", nBins=nb)])
        idx = (g["s"].to_numpy() // RES)
        ok = idx < len(v)
        vals = np.full(len(g), np.nan)
        vals[ok] = v[idx[ok]]
        out[g.index.to_numpy()] = vals
    bw.close()
    return (out - np.nanmean(out)) / np.nanstd(out)


def contacts_one(c: str, signs: pd.DataFrame) -> Path:
    out = WORK / f"ab_{c}.tsv"
    if out.exists() and out.stat().st_size:
        return out
    sf = WORK / f"sign_{c}.tsv"
    signs.to_csv(sf, sep="\t", header=False, index=False)
    cmd = (f"{STRAW} NONE {HIC} {c} {c} BP {RES} | "
           f"mawk -F'\\t' '{AWK}' {sf} - > {out}.tmp && mv {out}.tmp {out}")
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    return out


def contact_sums(cr: pd.DataFrame) -> pd.DataFrame:
    WORK.mkdir(exist_ok=True)
    jobs = {}
    with ThreadPoolExecutor(NJOBS) as ex:
        for c in CHROMS:
            g = cr[cr["c"] == c]
            if len(g):
                signs = pd.DataFrame({"s": g["s"], "sign": np.where(g["cr"] > 0, "A", "B")})
                jobs[c] = ex.submit(contacts_one, c, signs)
        frames = []
        for c, f in jobs.items():
            t = pd.read_csv(f.result(), sep="\t", header=None, names=["s", "toA", "toB"])
            t["c"] = c
            frames.append(t)
    return pd.concat(frames)


def main() -> None:
    d, cr = load_calls()
    d["gro"] = gro(d)
    ab = contact_sums(cr)
    d = d.merge(ab, on=["c", "s"], how="left")

    res = {}
    for k in CLASSES:
        s = d[d["cls"] == k]
        res[k] = {"n_bins": int(len(s)),
                  "gro_z_median": float(s["gro"].median()),
                  "contacts_to_A": float(s["toA"].sum()),
                  "contacts_to_B": float(s["toB"].sum()),
                  "contact_ratio_log2": float(np.log2(s["toA"].sum() / s["toB"].sum()))}
    (HERE / "fig2B.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))

    fig, ax = plt.subplots(3, 1, figsize=(3.0, 6.4),
                           gridspec_kw={"height_ratios": [0.8, 1.2, 1.2], "hspace": 0.55})
    n = {k: res[k]["n_bins"] for k in CLASSES}
    tot = sum(n.values())
    left = 0
    for k in CLASSES:
        ax[0].barh(0, n[k] / tot, left=left, color=COLORS[k], height=0.6)
        if n[k] / tot > 0.05:
            ax[0].text(left + n[k] / tot / 2, 0, f"{n[k] / tot:.0%}", ha="center",
                       va="center", fontsize=7, color="white")
        left += n[k] / tot
    ax[0].set_xlim(0, 1)
    ax[0].set_yticks([])
    ax[0].set_xticks([])
    ax[0].set_title(f"bins (n = {tot:,})", fontsize=8)
    for s in ("top", "right", "left", "bottom"):
        ax[0].spines[s].set_visible(False)

    data = [d.loc[d["cls"] == k, "gro"].dropna().to_numpy() for k in CLASSES]
    bp = ax[1].boxplot(data, showfliers=False, patch_artist=True, widths=0.6,
                       medianprops=dict(color="black"))
    for b, k in zip(bp["boxes"], CLASSES):
        b.set_facecolor(COLORS[k])
    ax[1].set_ylabel("GRO-seq\n(z-score)", fontsize=8)
    ax[1].set_xticks([])

    r = [res[k]["contact_ratio_log2"] for k in CLASSES]
    ax[2].bar(range(1, 5), r, color=[COLORS[k] for k in CLASSES], width=0.6)
    ax[2].axhline(0, color="black", lw=0.6)
    ax[2].set_ylabel("Contact ratio\n(log2 A/B)", fontsize=8)
    ax[2].set_xticks(range(1, 5))
    ax[2].set_xticklabels(CLASSES, fontsize=8)
    ax[2].set_xlabel("Eigen. -> CRUSH", fontsize=8)
    for a in ax[1:]:
        a.tick_params(labelsize=7)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(HERE / f"fig2B.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
