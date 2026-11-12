#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import csv
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
REV = FIG.parent
sys.path.insert(0, str(REV))
import MSD_tracks as msd

OUT = FIG / "results"
XLSX = REV / "SourceData_revision_UPDATED.xlsx"
TAB = "DATA_S6_Fig.1D"

RED, BLUE = "#C6303A", "#2E77B5"
INK, INK2, GRID = "#16181c", "#5a5f68", "#e6e6e2"

SAMPLES = ["16_1", "16_2", "30_1", "30_2", "40_1", "40_2"]
BLOCKS = [("10 kb Eigen", "eigen", 10000, 23),
          ("10 kb CRUSH", "crush", 10000, 14),
          ("1 kb CRUSH",  "crush",  1000,  5)]


def crush_track(s: str, r: int) -> Path:
    return REV / "tracks" / "wi38" / f"Wi38_{s}_r{r}" / f"Wi38_{s}_mergedCrush_{r}.bedgraph"


def eigen_track(s: str, r: int) -> Path:
    return REV / "tracks" / "wi38_eigen" / f"Wi38_{s}_eigen_{r}.bedgraph"


_cache: dict[Path, object] = {}


def track(p: Path):
    if p not in _cache:
        if not p.exists():
            raise SystemExit(f"missing track: {p}")
        _cache[p] = msd.keyed(msd.load(str(p)), literal=False)
    return _cache[p]


def inv_sqrt_msd(a: Path, b: Path) -> float:
    """1/sqrt(MSD) on the bare IQR scale."""
    tag = f"{a.name} vs {b.name}"
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        m = msd.align(track(a), track(b), False, tag)
        m = msd.orient_to_control(m, False, 100, tag)
        xi, yi = msd.track_iqr(m)
    if xi == 0 or yi == 0:
        raise SystemExit(f"zero IQR for {tag}")
    d = ((m["Score_x"] / xi) - (m["Score_y"] / yi)) ** 2
    s = float(np.sqrt(d.mean()))
    return 1.0 / s if s else float("nan")


def measure():
    """One record per block: labelled pairs split into replicates and passages."""
    out = []
    for title, method, res, _row in BLOCKS:
        get = crush_track if method == "crush" else eigen_track
        pairs = []
        print(f"  {title}:")
        for i in range(len(SAMPLES)):
            for j in range(i + 1, len(SAMPLES)):
                a, b = SAMPLES[i], SAMPLES[j]
                v = inv_sqrt_msd(get(a, res), get(b, res))
                same = a.split("_")[0] == b.split("_")[0]
                pairs.append({"pair": f"{a} vs {b}",
                              "comparison": "replicates" if same else "passages",
                              "inv_sqrt_MSD": v})
        rep = [p["inv_sqrt_MSD"] for p in pairs if p["comparison"] == "replicates"]
        pas = [p["inv_sqrt_MSD"] for p in pairs if p["comparison"] == "passages"]
        print(f"    replicates n={len(rep):<3} mean {np.mean(rep):7.3f}  SD {np.std(rep, ddof=1):6.3f}")
        print(f"    passages   n={len(pas):<3} mean {np.mean(pas):7.3f}  SD {np.std(pas, ddof=1):6.3f}")
        print(f"    ratio rep/pass {np.mean(rep)/np.mean(pas):.2f}", flush=True)
        out.append({"title": title, "colour": RED if method == "crush" else BLUE,
                    "pairs": pairs, "rep": rep, "pas": pas})
    return out


def check_workbook(blocks):
    """DATA S6 holds MSD on the IQR x100 scale; bare 1/sqrt is 100/sqrt(v)."""
    if not XLSX.exists():
        print("  workbook not found -- skipping the cross-check")
        return
    ws = openpyxl.load_workbook(XLSX, data_only=True)[TAB]
    bad = tot = 0
    for blk, (title, _method, _res, r0) in zip(blocks, BLOCKS):
        got = {p["pair"]: p["inv_sqrt_MSD"] for p in blk["pairs"]}
        for i in range(5):
            for c in range(2 + i, 7):
                v = ws.cell(row=r0 + i, column=c).value
                if not isinstance(v, (int, float)) or v <= 0:
                    continue
                key = f"{SAMPLES[i]} vs {SAMPLES[c - 1]}"
                theirs = 100.0 / float(np.sqrt(float(v)))
                mine = got.get(key)
                tot += 1
                if mine is None or abs(mine - theirs) > 0.005 * max(abs(theirs), 1e-9):
                    print(f"  DIFFERS  {title} {key}: recomputed {mine!r} "
                          f"vs sheet {theirs:.4f}")
                    bad += 1
    print(f"  workbook cross-check: {tot - bad}/{tot} pairs agree")


def sourcedata(blocks):
    rows = []
    for blk in blocks:
        for p in blk["pairs"]:
            rows.append({"block": blk["title"], "pair": p["pair"],
                         "comparison": p["comparison"],
                         "inv_sqrt_MSD": p["inv_sqrt_MSD"]})
        for comparison, vals in (("replicates", blk["rep"]), ("passages", blk["pas"])):
            for stat, v in (("MEAN", float(np.mean(vals))),
                            ("SD", float(np.std(vals, ddof=1))),
                            ("n", len(vals))):
                rows.append({"block": blk["title"], "pair": stat,
                             "comparison": comparison, "inv_sqrt_MSD": v})
    p = OUT / "fig1D_bar_sourcedata.csv"
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["block", "pair", "comparison",
                                           "inv_sqrt_MSD"])
        w.writeheader(); w.writerows(rows)
    print(f"wrote {p}")


def figure(blocks):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    rng = np.random.default_rng(0)
    w, gap = 0.36, 0.10
    for k, blk in enumerate(blocks):
        col = blk["colour"]
        for m, (vals, solid) in enumerate(((blk["rep"], True), (blk["pas"], False))):
            x = k + (m - 0.5) * (w + gap)
            vals = np.asarray(vals, float)
            mu, sd = float(vals.mean()), float(vals.std(ddof=1))
            face = col if solid else (col + "33")
            ax.bar(x, mu, w, facecolor=face, edgecolor=col, linewidth=1.6, zorder=2)
            ax.errorbar(x, mu, yerr=sd, fmt="none", ecolor=INK2, elinewidth=2.0,
                        capsize=7, capthick=2.0, zorder=4)
            ax.scatter(x + rng.uniform(-0.10, 0.10, len(vals)), vals, s=46,
                       facecolor="white", edgecolor=INK, linewidth=1.5,
                       zorder=5, clip_on=False)
        ax.text(k, -0.9, blk["title"], ha="center", va="top", fontsize=10.5,
                color=INK, clip_on=False)
    ax.set_xlim(-0.6, len(blocks) - 0.4)
    ax.set_xticks([])
    ax.set_ylabel("1/sqrt(MSD)", fontsize=10.5, color=INK, labelpad=8)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#d7d7d3")
    ax.tick_params(colors=INK2, labelsize=9.5, length=0)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(bottom=0.16, top=0.96, left=0.11, right=0.98)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig1D_bar.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"\nwrote {OUT/'fig1D_bar.png'} (+ .pdf)")


def main() -> None:
    print("measuring with MSD_tracks.py:")
    blocks = measure()
    OUT.mkdir(parents=True, exist_ok=True)
    print()
    sourcedata(blocks)
    check_workbook(blocks)
    figure(blocks)


if __name__ == "__main__":
    main()
