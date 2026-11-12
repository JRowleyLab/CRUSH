#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import io
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
from matplotlib.ticker import FixedLocator, NullFormatter

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
REV = FIG.parent
sys.path.insert(0, str(REV))
import MSD_tracks as msd

OUT = FIG / "results"
XLSX = REV / "SourceData_revision_UPDATED.xlsx"
TAB = "DATA_S5_Fig.S1E"
CRUSHDEV = Path("CRUSHdev")

RED, BLUE = "#C6303A", "#2E77B5"
INK, INK2, GRID = "#16181c", "#5a5f68", "#e6e6e2"

LADDER = [100000, 25000, 10000, 5000, 1000]
RLAB = {100000: "100 kb", 25000: "25 kb", 10000: "10 kb",
        5000: "5 kb", 1000: "1 kb"}
CELLS = [("H1", "H1_insitu", "H1_microC"),
         ("HFFc6", "HFFc6_insitu", "HFFc6_microC")]
COL = {("H1", "CRUSH"): 1, ("H1", "Eigen"): 2,
       ("HFFc6", "CRUSH"): 3, ("HFFc6", "Eigen"): 4}


def crush_track(s: str, r: int) -> Path:
    return REV / "tracks" / "h1hff" / f"{s}_r{r}" / f"{s}_mergedCrush_{r}.bedgraph"


def eigen_track(s: str, r: int) -> Path:
    return CRUSHDEV / "h1hff" / "eigen" / f"{s}_eigen_{r}.bedgraph"


_cache: dict[Path, object] = {}


def track(p: Path):
    if p not in _cache:
        if not p.exists():
            raise SystemExit(f"missing track: {p}")
        _cache[p] = msd.keyed(msd.load(str(p)), literal=False)
    return _cache[p]


def inv_sqrt_msd(hic: Path, microc: Path) -> float:
    """1/sqrt(MSD), bare IQR scale. Hi-C is the control, Micro-C the sample."""
    tag = f"{hic.name} vs {microc.name}"
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        m = msd.align(track(hic), track(microc), False, tag)
        m = msd.orient_to_control(m, False, 100, tag)
        xi, yi = msd.track_iqr(m)
    if xi == 0 or yi == 0:
        raise SystemExit(f"zero IQR for {tag}")
    d = ((m["Score_x"] / xi) - (m["Score_y"] / yi)) ** 2
    s = float(np.sqrt(d.mean()))
    return 1.0 / s if s else float("nan")


def measure():
    rows = []
    for r in LADDER:
        rec = {"res": r, "label": RLAB[r]}
        for cell, hic, mic in CELLS:
            for method, get in (("CRUSH", crush_track), ("Eigen", eigen_track)):
                rec[f"{cell}_{method}"] = inv_sqrt_msd(get(hic, r), get(mic, r))
        rows.append(rec)
        print(f"  {RLAB[r]:>7}  " + "   ".join(
            f"{c}/{m} {rec[f'{c}_{m}']:8.4f}" for c, _, _ in CELLS
            for m in ("CRUSH", "Eigen")), flush=True)
    return rows


def check_workbook(rows):
    if not XLSX.exists():
        print("  workbook not found -- skipping the cross-check")
        return
    ws = openpyxl.load_workbook(XLSX, data_only=True)[TAB]
    pub = {}
    want = {RLAB[r] for r in LADDER}
    for row in ws.iter_rows(values_only=True):
        lab = str(row[0]).strip() if row and row[0] is not None else ""
        if lab in want:
            pub[lab] = list(row[1:5])
    missing = want - set(pub)
    if missing:
        print(f"  workbook rows not found: {sorted(missing)}")
    bad = tot = 0
    for rec in rows:
        got = pub.get(rec["label"])
        if got is None:
            continue
        for (cell, method), idx in COL.items():
            v = got[idx - 1]
            if not isinstance(v, (int, float)) or v <= 0:
                continue
            tot += 1
            theirs = 100.0 / float(v)
            mine = rec[f"{cell}_{method}"]
            if abs(mine - theirs) > 0.005 * max(abs(theirs), 1e-9):
                print(f"  DIFFERS  {rec['label']:>7} {cell} {method}: "
                      f"recomputed {mine:.5f} vs sheet {theirs:.5f}")
                bad += 1
    print(f"  workbook cross-check: {tot - bad}/{tot} points agree")


def figure(rows):
    x = np.array([r["res"] for r in rows], float)
    series = [("H1_CRUSH", RED, "-"), ("HFFc6_CRUSH", RED, (0, (5, 2))),
              ("H1_Eigen", BLUE, "-"), ("HFFc6_Eigen", BLUE, (0, (5, 2)))]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    for key, col, dash in series:
        y = np.array([r[key] for r in rows], float)
        ax.plot(x, y, color=col, linestyle=dash, linewidth=6.0, marker="o",
                markersize=11.0, markerfacecolor=col, markeredgecolor="white",
                markeredgewidth=2.0, zorder=3, clip_on=False,
                solid_capstyle="round", dash_capstyle="round")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
    ax.xaxis.set_major_locator(FixedLocator(x))
    ax.set_xticklabels([r["label"] for r in rows])
    ax.xaxis.set_minor_formatter(NullFormatter())
    yt = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]
    ax.yaxis.set_major_locator(FixedLocator(yt))
    ax.set_yticklabels([f"{v:g}" for v in yt])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(0.015, 11)
    ax.set_xlabel("Resolution", fontsize=10.5, color=INK, labelpad=8)
    ax.set_ylabel("1/sqrt(MSD)", fontsize=10.5, color=INK, labelpad=8)
    ax.grid(True, which="major", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d7d7d3")
    ax.tick_params(colors=INK2, labelsize=9.5, length=0)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"figS1F_line.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"\nwrote {OUT/'figS1F_line.png'} (+ .pdf)")


def main() -> None:
    print("measuring with MSD_tracks.py:")
    rows = measure()
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "figS1F_line_sourcedata.json"
    p.write_text(json.dumps(rows, indent=1))
    print(f"\nwrote {p}")
    check_workbook(rows)
    figure(rows)
    for key, _, _ in [("H1_CRUSH", 0, 0), ("HFFc6_CRUSH", 0, 0),
                      ("H1_Eigen", 0, 0), ("HFFc6_Eigen", 0, 0)]:
        a, b = rows[0][key], rows[-1][key]
        print(f"  {key:<14} 100 kb -> 1 kb   {a/b:8.1f}x fall")


if __name__ == "__main__":
    main()
