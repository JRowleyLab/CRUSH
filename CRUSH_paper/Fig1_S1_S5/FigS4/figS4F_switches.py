#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
sys.path.insert(0, str(FIG.parent))
import MSD_tracks as msd

T = FIG.parent / "tracks" / "chicken"
GENES = Path("chicken_nochr_genes.bed")
PAIRS = [("Fibroblast", "Erythroblast", "Fibro→Erythro"),
         ("Erythroblast", "eRBC", "Erythro→eRBC\n(early)"),
         ("eRBC", "aRBC", "eRBC→aRBC\n(late)")]


def gene_density(res):
    g = pd.read_csv(GENES, sep="\t", header=None, usecols=[0, 1, 3],
                    names=["c", "s", "kind"], dtype={"c": str})
    g = g[g["kind"] == "transcript"]
    g["b"] = (g["s"] // res).astype(np.int64)
    return g.groupby(["c", "b"]).size().rename("n").reset_index()


def load(sample, res, dens):
    p = T / f"{sample}_r{res}" / f"{sample}_mergedCrush_{res}.bedgraph"
    d = msd.load(str(p)).dropna(subset=["Score"])
    t = pd.DataFrame({"c": msd.norm_chrom(d["Chromosome"]).to_numpy(),
                      "b": (d["Start"].to_numpy() // res).astype(np.int64),
                      "v": d["Score"].to_numpy()})
    t = t[t["v"] != 0]
    m = t.merge(dens, on=["c", "b"], how="left").fillna({"n": 0.0})
    flipped = []
    for c, sub in m.groupby("c"):
        if len(sub) < 50 or sub["n"].std() == 0:
            continue
        r = np.corrcoef(sub["v"], sub["n"])[0, 1]
        if np.isfinite(r) and r < 0:
            t.loc[t["c"] == c, "v"] *= -1
            flipped.append(c)
    return t, flipped


def counts(res):
    dens = gene_density(res)
    tr = {}
    for s in ("Fibroblast", "Erythroblast", "eRBC", "aRBC"):
        tr[s], fl = load(s, res, dens)
        print(f"    {s:<13} oriented: flipped {len(fl)} chrom"
              f"{'s' if len(fl) != 1 else ''}"
              f"{' (' + ','.join(fl[:8]) + ')' if fl else ''}", flush=True)
    out = []
    for a, b, lab in PAIRS:
        m = tr[a].merge(tr[b], on=["c", "b"], suffixes=("_a", "_b"))
        sa, sb = np.sign(m["v_a"]), np.sign(m["v_b"])
        ab = int(((sa > 0) & (sb < 0)).sum())
        ba = int(((sa < 0) & (sb > 0)).sum())
        out.append({"label": lab, "shared": len(m), "AtoB": ab, "BtoA": ba,
                    "pct": 100 * (ab + ba) / len(m)})
    return out


def plot(res=1000, labels=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = json.loads((OUT / "figS4F_switches.json").read_text())[str(res)]
    rows = [r for r in rows if "early" in r["label"] or "late" in r["label"]]
    AB, BA, INK, INK2 = "#45BFDE", "#E9483D", "#16181c", "#5a5f68"
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    x = np.arange(len(rows)); w = 0.52
    ab = np.array([r["AtoB"] / 1000 for r in rows])
    ba = np.array([r["BtoA"] / 1000 for r in rows])
    y = x[::-1]
    ax.barh(y, ab, w, color=AB, lw=0, label="A \u2192 B", zorder=3)
    ax.barh(y, ba, w, left=ab, color=BA, lw=0, label="B \u2192 A", zorder=3)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_yticks(y)
    if labels:
        for yi, a, b in zip(y, ab, ba):
            if a > 0:
                ax.text(a / 2, yi, f"{a:,.0f}", ha="center", va="center",
                        fontsize=9.5, color="white", zorder=4)
            if b > 0:
                ax.text(a + b / 2, yi, f"{b:,.0f}", ha="center", va="center",
                        fontsize=9.5, color="white", zorder=4)
            ax.text(a + b + (ab + ba).max() * 0.02, yi, f"{a + b:,.0f}",
                    ha="left", va="center", fontsize=9.5, color=INK2, zorder=4)
        ax.set_yticklabels([r["label"].replace(chr(10), " ") for r in rows],
                           fontsize=10, color=INK)
        ax.set_xlabel("compartment switches (thousands of 1 kb bins)",
                      fontsize=10, color=INK, labelpad=8)
        ax.legend(frameon=False, fontsize=9.5, handlelength=1.2, loc="lower right")
        ax.tick_params(colors=INK2, labelsize=9.5, length=0)
        ax.set_axisbelow(True); ax.xaxis.grid(True, color="#f0f0ee", lw=0.8)
        ax.set_xlim(0, (ab + ba).max() * 1.16)
    else:
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.tick_params(length=0)
        ax.set_xlim(0, (ab + ba).max() * 1.02)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color("#d7d7d3" if labels else "none")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    stem = "figS4F_switches" + ("_labeled" if labels else "")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT / (stem + '.png')} (+ .pdf)")
    import csv
    q = OUT / "figS4F_switches_sourcedata.csv"
    with q.open("w", newline="") as fh:
        wcsv = csv.writer(fh)
        wcsv.writerow(["transition", "resolution_bp", "shared_bins",
                       "A_to_B", "B_to_A", "total_switches", "percent_of_shared"])
        for r in rows:
            wcsv.writerow([r["label"].replace(chr(10), " "), res, r["shared"],
                           r["AtoB"], r["BtoA"], r["AtoB"] + r["BtoA"],
                           f"{r['pct']:.2f}"])
    print(f"wrote {q}")


if __name__ == "__main__":
    res_all = {}
    for res in (1000, 50000):
        print(f"\n=== {res // 1000} kb ===")
        res_all[res] = counts(res)
        print(f"  {'transition':<22}{'shared':>10}{'A->B':>10}{'B->A':>10}"
              f"{'total':>10}{'%':>7}")
        for r in res_all[res]:
            print(f"  {r['label'].replace(chr(10), ' '):<22}{r['shared']:>10,}"
                  f"{r['AtoB']:>10,}{r['BtoA']:>10,}"
                  f"{r['AtoB'] + r['BtoA']:>10,}{r['pct']:>7.1f}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figS4F_switches.json").write_text(
        json.dumps({str(k): v for k, v in res_all.items()}, indent=1))
    print(f"\nwrote {OUT / 'figS4F_switches.json'}")
    plot(1000, labels=False)
    plot(1000, labels=True)
