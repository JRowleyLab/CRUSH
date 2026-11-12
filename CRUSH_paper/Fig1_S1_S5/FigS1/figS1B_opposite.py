#!/usr/bin/env python3
from __future__ import annotations

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
import opposite_calls as oc

LCL = FIG / "LCL"
CRUSH_DIR, EIGEN_DIR = LCL / "CRUSHtracks", LCL / "EIGENtracks"
OUT = FIG / "results"
XLSX = REV / "SourceData_revision_UPDATED.xlsx"
TAB = "DATA S1_Fig. S1B"

RED, BLUE, GREY = "#C6303A", "#2E77B5", "#9aa0a8"
INK, INK2, GRID = "#16181c", "#5a5f68", "#e6e6e2"

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
        _cache[p] = oc.load(str(p), literal=False)
    return _cache[p]


def opposite(sub: Path, ctl: Path) -> float:
    """Percent of calling bins where the two tracks disagree, oc's own defaults.

    The full-depth track is the CONTROL and the subsampled one the sample, so the
    per-chromosome orientation flips the subsample onto the reference rather than
    the other way round. For this metric the choice is symmetric -- both orders
    give 48.211% on the worst pair -- but it is the direction the quantity means.
    """
    a, b = ctl, sub
    m = oc.align(track(a), track(b), False, f"{b.name} vs {a.name}")
    m = oc.orient_to_control(m, False, 100, 0.0, f"{b.name} vs {a.name}")
    called, opp, _one, _zero = oc.call_masks(m, 0.0)
    n = int(called.sum())
    if not n:
        raise SystemExit(f"no calling bins for {b.name} vs {a.name}")
    return 100.0 * int(opp.sum()) / n


def measure():
    cfull = CRUSH_DIR / f"{DEPTHS[0][2]}.bw"
    efull = EIGEN_DIR / f"{DEPTHS[0][3]}.bedgraph"
    rows = []
    for lab, contacts, cstem, estem in DEPTHS:
        c, e = CRUSH_DIR / f"{cstem}.bw", EIGEN_DIR / f"{estem}.bedgraph"
        self_e = None if lab == "20.3 B" else opposite(e, efull)
        self_c = None if lab == "20.3 B" else opposite(c, cfull)
        cross = opposite(c, efull)
        rows.append({"depth": lab, "contacts": contacts, "eigen_vs_eigen_full": self_e,
                     "crush_vs_crush_full": self_c, "crush_vs_eigen_full": cross})
        print(f"  {lab:>7}  eigen/eigen "
              f"{'   NA' if self_e is None else f'{self_e:5.2f}%'}   "
              f"CRUSH/CRUSH {'   NA' if self_c is None else f'{self_c:5.2f}%'}   "
              f"CRUSH/eigen {cross:5.2f}%", flush=True)
    return rows


def check_workbook(rows):
    """Compare against DATA S1, which holds the published fractions (not %)."""
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
    for row in rows:
        got = pub.get(row["depth"])
        if got is None:
            print(f"  {row['depth']}: no matching workbook row")
            continue
        for k, v in zip(keys, got):
            mine, theirs = row[k], (None if v in (None, "NA") else float(v) * 100)
            if mine is None and theirs is None:
                continue
            if mine is None or theirs is None or abs(mine - theirs) > 0.15:
                print(f"  DIFFERS  {row['depth']:>7} {k:<22} "
                      f"recomputed {mine!r} vs sheet {theirs!r}")
                bad += 1
    print(f"  workbook cross-check: {'all points agree' if not bad else f'{bad} disagree'}")


def figure(rows):
    x = np.array([r["contacts"] for r in rows]) / 1e6
    series = [("eigen_vs_eigen_full", "Eigenvector vs eigenvector @ 20.3 B", BLUE, "-"),
              ("crush_vs_crush_full", "CRUSH vs CRUSH @ 20.3 B", RED, "-"),
              ("crush_vs_eigen_full", "CRUSH vs eigenvector @ 20.3 B", GREY, (0, (5, 2)))]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    for key, lab, col, dash in series:
        y = np.array([np.nan if r[key] is None else r[key] for r in rows], float)
        ok = np.isfinite(y)
        ax.plot(x[ok], y[ok], color=col, linestyle=dash, linewidth=6.0, marker="o",
                markersize=11.0, markerfacecolor=col, markeredgecolor="white",
                markeredgewidth=2.0, label=lab, zorder=3, clip_on=False,
                solid_capstyle="round", dash_capstyle="round")
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.xaxis.set_major_locator(FixedLocator(x))
    ax.set_xticklabels([r["depth"] for r in rows])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(0, 55)
    ax.set_xlabel("Intrachromosomal contacts", fontsize=10.5, color=INK, labelpad=8)
    ax.set_ylabel("Bins with opposite calls (%)", fontsize=10.5, color=INK, labelpad=8)
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
        fig.savefig(OUT / f"figS1B_opposite.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"\nwrote {OUT/'figS1B_opposite.png'} (+ .pdf)")


def main() -> None:
    print("measuring with opposite_calls.py:")
    rows = measure()
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "figS1B_opposite_sourcedata.json"
    p.write_text(json.dumps(rows, indent=1))
    print(f"\nwrote {p}")
    check_workbook(rows)
    figure(rows)


if __name__ == "__main__":
    main()
