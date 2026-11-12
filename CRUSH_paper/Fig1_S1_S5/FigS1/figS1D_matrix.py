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
import pandas as pd

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
REV = FIG.parent
sys.path.insert(0, str(REV))
import MSD_tracks as msd

LCL = FIG / "LCL"
OUT = FIG / "results"
XLSX = REV / "SourceData_revision_UPDATED.xlsx"
TAB = "DATA_S3_Fig. S1D"

INK, INK2 = "#16181c", "#5a5f68"
VMIN, VMAX = 0.0, 60.0
DARKRED = np.array([139, 0, 0]) / 255.0

RES = [1000, 5000, 10000, 25000, 100000]
LABEL = {1000: "1 kb", 5000: "5 kb", 10000: "10 kb",
         25000: "25 kb", 100000: "100 kb"}
TRACKS = {
    "CRUSH": {r: LCL / "CRUSHtracks" / f"full_Crush_{r}.bw" for r in RES},
    "Eigen": {1000:   LCL / "EIGENtracks" / "e_full.bedgraph",
              5000:   LCL / "EIGENtracks" / "e_5kb.bedgraph",
              10000:  LCL / "EIGENtracks" / "e_10kb.bedgraph",
              25000:  LCL / "EIGENtracks" / "e_25kb.bedgraph",
              100000: LCL / "EIGENtracks" / "e_100kb.bedgraph"},
}

_raw: dict[Path, pd.DataFrame] = {}
_binned: dict[tuple[Path, int], pd.DataFrame] = {}


def raw(p: Path) -> pd.DataFrame:
    if p not in _raw:
        if not p.exists():
            raise SystemExit(f"missing track: {p}")
        _raw[p] = msd.load(str(p))
    return _raw[p]


def at(p: Path, res: int, native: int) -> pd.DataFrame:
    """The track on a `res` bp grid; averaged down when res is coarser."""
    if res == native:
        return raw(p)
    key = (p, res)
    if key not in _binned:
        d = raw(p).copy()
        d["Start"] = (d["Start"] // res) * res
        g = d.groupby(["Chromosome", "Start"], as_index=False)["Score"].mean()
        g["End"] = g["Start"] + res
        _binned[key] = g[["Chromosome", "Start", "End", "Score"]]
    return _binned[key]


def similarity(method: str, r_fine: int, r_coarse: int) -> float:
    """1/sqrt(MSD), bare IQR scale, both tracks on the coarser grid."""
    pf, pc = TRACKS[method][r_fine], TRACKS[method][r_coarse]
    a = at(pf, r_coarse, r_fine)
    b = at(pc, r_coarse, r_coarse)
    tag = f"{method} {LABEL[r_fine]} vs {LABEL[r_coarse]}"
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        m = msd.align(msd.keyed(b, False), msd.keyed(a, False), False, tag)
        m = msd.orient_to_control(m, False, 100, tag)
        xi, yi = msd.track_iqr(m)
    if xi == 0 or yi == 0:
        raise SystemExit(f"zero IQR for {tag}")
    d = ((m["Score_x"] / xi) - (m["Score_y"] / yi)) ** 2
    s = float(np.sqrt(d.mean()))
    return 1.0 / s if s else float("nan")


def measure():
    out = {}
    for method in ("CRUSH", "Eigen"):
        M = np.full((len(RES), len(RES)), np.nan)
        print(f"  {method}:")
        for i, rf in enumerate(RES):
            for j, rc in enumerate(RES):
                if j <= i:
                    continue
                v = similarity(method, rf, rc)
                M[i, j] = M[j, i] = v
                print(f"    {LABEL[rf]:>6} vs {LABEL[rc]:<6} {v:8.4f}", flush=True)
        out[method] = M
    return out


def check_workbook(mats):
    if not XLSX.exists():
        print("  workbook not found -- skipping the cross-check")
        return
    ws = openpyxl.load_workbook(XLSX, data_only=True)[TAB]
    rows = [[c for c in r] for r in ws.iter_rows(values_only=True)]
    bad = tot = 0
    for method in ("CRUSH", "Eigen"):
        hdr = next((k for k, r in enumerate(rows)
                    if r and str(r[0]).strip().upper() == method.upper()), None)
        if hdr is None:
            print(f"  no {method} block found in {TAB}")
            continue
        for i in range(4):
            r = rows[hdr + 1 + i]
            for j in range(i + 1, 5):
                v = r[j] if j < len(r) else None
                if v is None:
                    continue
                tot += 1
                mine = mats[method][i, j]
                if abs(mine - float(v)) > 0.005 * max(abs(float(v)), 1.0):
                    print(f"  DIFFERS  {method} {LABEL[RES[i]]} vs {LABEL[RES[j]]}: "
                          f"recomputed {mine:.4f} vs sheet {float(v):.4f}")
                    bad += 1
    print(f"  workbook cross-check: {tot - bad}/{tot} cells agree")


def figure(mats):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.9))
    fig.patch.set_facecolor("white")
    n = len(RES)
    for ax, method in zip(axes, ("CRUSH", "Eigen")):
        M = mats[method]
        ax.set_facecolor("white")
        for i in range(n):
            for j in range(n):
                if j <= i or not np.isfinite(M[i, j]):
                    continue
                t = np.clip((M[i, j] - VMIN) / (VMAX - VMIN), 0, 1)
                fill = (1 - t) * np.ones(3) + t * DARKRED
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=fill,
                                           edgecolor="white", linewidth=2.0))
                lum = 0.2126 * fill[0] + 0.7152 * fill[1] + 0.0722 * fill[2]
                ax.text(j + 0.5, i + 0.5, f"{M[i, j]:.2f}", ha="center",
                        va="center", fontsize=11,
                        color="white" if lum < 0.55 else INK)
        ax.set_xlim(0, n); ax.set_ylim(n, 0)
        ax.set_xticks(np.arange(n) + 0.5); ax.set_yticks(np.arange(n) + 0.5)
        ax.set_xticklabels([LABEL[r] for r in RES], fontsize=10, color=INK2)
        ax.set_yticklabels([LABEL[r] for r in RES], fontsize=10, color=INK2)
        ax.set_aspect("equal")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(method, fontsize=11.5, color=INK, pad=10)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"figS1D_matrix.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"\nwrote {OUT/'figS1D_matrix.png'} (+ .pdf)")


def main() -> None:
    print("measuring with MSD_tracks.py:")
    mats = measure()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [{"method": m, "res_a": LABEL[RES[i]], "res_b": LABEL[RES[j]],
             "inv_sqrt_MSD": mats[m][i, j]}
            for m in ("CRUSH", "Eigen")
            for i in range(len(RES)) for j in range(i + 1, len(RES))]
    p = OUT / "figS1D_matrix_sourcedata.json"
    p.write_text(json.dumps(rows, indent=1))
    print(f"\nwrote {p}")
    check_workbook(mats)
    figure(mats)


if __name__ == "__main__":
    main()
