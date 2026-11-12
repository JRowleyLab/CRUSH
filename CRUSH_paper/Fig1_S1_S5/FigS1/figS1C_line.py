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

LCL = FIG / "LCL"
CRUSH_DIR, EIGEN_DIR = LCL / "CRUSHtracks", LCL / "EIGENtracks"
OUT = FIG / "results"
XLSX = REV / "SourceData_revision_UPDATED.xlsx"
TAB = "DATA_S2_Fig. S1C"

RED, BLUE, GREY = "#C6303A", "#2E77B5", "#9aa0a8"
INK, INK2, GRID = "#16181c", "#5a5f68", "#e6e6e2"
SIMILARITY = False
X100 = 100.0

DEPTHS = [("20.3 B", 20.3e9, "full_mergedCrush_1000", "e_full"),
          ("5 B",     5.0e9, "5B_mergedCrush_1000",   "e_5bil"),
          ("1 B",     1.0e9, "1B_mergedCrush_1000",   "e_1bil"),
          ("500 M",   5.0e8, "500M_mergedCrush_1000", "e_500mil"),
          ("250 M",   2.5e8, "250M_mergedCrush_1000", "e_250mil"),
          ("100 M",   1.0e8, "100M_mergedCrush_1000", "e_100mil")]

_cache: dict[Path, object] = {}


def track(p: Path):
    if p not in _cache:
        if not p.exists():
            raise SystemExit(f"missing track: {p}")
        _cache[p] = msd.keyed(msd.load(str(p)), literal=False)
    return _cache[p]


def sqrt_msd(sub: Path, ctl: Path) -> float:
    """sqrt(MSD) on the IQR x100 scale, control = the full-depth track."""
    tag = f"{sub.name} vs {ctl.name}"
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        m = msd.align(track(ctl), track(sub), False, tag)
        m = msd.orient_to_control(m, False, 100, tag)
        xi, yi = msd.track_iqr(m)
    if xi == 0 or yi == 0:
        raise SystemExit(f"zero IQR for {tag}")
    d = ((m["Score_x"] / xi) * X100 - (m["Score_y"] / yi) * X100) ** 2
    return float(np.sqrt(d.mean()))


def measure():
    cfull = CRUSH_DIR / f"{DEPTHS[0][2]}.bw"
    efull = EIGEN_DIR / f"{DEPTHS[0][3]}.bedgraph"
    rows = []
    for lab, contacts, cstem, estem in DEPTHS:
        c, e = CRUSH_DIR / f"{cstem}.bw", EIGEN_DIR / f"{estem}.bedgraph"
        self_e = None if lab == "20.3 B" else sqrt_msd(e, efull)
        self_c = None if lab == "20.3 B" else sqrt_msd(c, cfull)
        cross = sqrt_msd(c, efull)
        rows.append({"depth": lab, "contacts": contacts,
                     "eigen_vs_eigen_full": self_e,
                     "crush_vs_crush_full": self_c,
                     "crush_vs_eigen_full": cross})
        f = lambda v: "      NA" if v is None else f"{v:8.2f}"
        print(f"  {lab:>7}  eigen/eigen {f(self_e)}   CRUSH/CRUSH {f(self_c)}"
              f"   CRUSH/eigen {f(cross)}", flush=True)
    return rows


def check_workbook(rows):
    if not XLSX.exists():
        print("  workbook not found -- skipping the cross-check")
        return
    ws = openpyxl.load_workbook(XLSX, data_only=True)[TAB]
    want = {d[0] for d in DEPTHS}
    pub = {}
    for row in ws.iter_rows(values_only=True):
        lab = str(row[0]).strip() if row and row[0] is not None else ""
        if lab in want:
            pub[lab] = list(row[1:4])
    missing = want - set(pub)
    if missing:
        print(f"  workbook rows not found: {sorted(missing)}")
    keys = ["eigen_vs_eigen_full", "crush_vs_crush_full", "crush_vs_eigen_full"]
    bad = 0
    for r in rows:
        got = pub.get(r["depth"])
        if got is None:
            continue
        for k, v in zip(keys, got):
            mine = r[k]
            theirs = None if v in (None, "NA") else float(v)
            if mine is None and theirs is None:
                continue
            if (mine is None) != (theirs is None) or \
               abs(mine - theirs) > 0.005 * max(abs(theirs), 1.0):
                print(f"  DIFFERS  {r['depth']:>7} {k:<22} "
                      f"recomputed {mine!r} vs sheet {theirs!r}")
                bad += 1
    print(f"  workbook cross-check: {'all points agree' if not bad else f'{bad} disagree'}")


def figure(rows):
    x = np.array([r["contacts"] for r in rows]) / 1e6
    series = [("eigen_vs_eigen_full", BLUE, "-"),
              ("crush_vs_crush_full", RED, "-"),
              ("crush_vs_eigen_full", GREY, (0, (5, 2)))]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    lo, hi = np.inf, -np.inf
    for key, col, dash in series:
        y = np.array([np.nan if r[key] is None else r[key] for r in rows], float)
        if SIMILARITY:
            y = X100 / y
        ok = np.isfinite(y)
        lo, hi = min(lo, np.nanmin(y[ok])), max(hi, np.nanmax(y[ok]))
        ax.plot(x[ok], y[ok], color=col, linestyle=dash, linewidth=6.0, marker="o",
                markersize=11.0, markerfacecolor=col, markeredgecolor="white",
                markeredgewidth=2.0, zorder=3, clip_on=False,
                solid_capstyle="round", dash_capstyle="round")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
    ax.xaxis.set_major_locator(FixedLocator(x))
    ax.set_xticklabels([r["depth"] for r in rows])
    ax.xaxis.set_minor_formatter(NullFormatter())
    yt = [v for v in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
          if lo / 1.6 <= v <= hi * 1.6]
    ax.yaxis.set_major_locator(FixedLocator(yt))
    ax.set_yticklabels([f"{v:g}" for v in yt])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Intrachromosomal contacts", fontsize=10.5, color=INK, labelpad=8)
    ax.set_ylabel("100/sqrt(MSD)" if SIMILARITY else "sqrt(MSD)",
                  fontsize=10.5, color=INK, labelpad=8)
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
        fig.savefig(OUT / f"figS1C_line.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"\nwrote {OUT/'figS1C_line.png'} (+ .pdf)")


def main() -> None:
    print("measuring with MSD_tracks.py:")
    rows = measure()
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "figS1C_line_sourcedata.json"
    p.write_text(json.dumps(rows, indent=1))
    print(f"\nwrote {p}")
    check_workbook(rows)
    figure(rows)


if __name__ == "__main__":
    main()
