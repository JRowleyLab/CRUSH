#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import openpyxl

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
REV = FIG.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REV))
import heat_panels as hp
import MSD_tracks as msd

OUT = FIG / "results"
XLSX = REV / "SourceData_revision_UPDATED.xlsx"
SHEET = "DATA_S4_Fig.1C"
CRUSHDEV = Path("CRUSHdev")
RES = 1000
X100 = 100.0

SAMPLES = ["H1_insitu", "H1_microC", "HFFc6_insitu", "HFFc6_microC"]
LABELS = ["H1 Hi-C", "H1 Micro-C", "HFFc6 Hi-C", "HFFc6 Micro-C"]


def crush_track(s: str) -> Path:
    return REV / "tracks" / "h1hff" / f"{s}_r{RES}" / f"{s}_mergedCrush_{RES}.bedgraph"


def eigen_track(s: str) -> Path:
    return CRUSHDEV / "h1hff" / "eigen" / f"{s}_eigen_{RES}.bedgraph"


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
    out = {}
    for method, get in (("CRUSH", crush_track), ("Eigenvector", eigen_track)):
        cells = {}
        print(f"  {method}:")
        for i in range(len(SAMPLES)):
            for j in range(i + 1, len(SAMPLES)):
                v = inv_sqrt_msd(get(SAMPLES[i]), get(SAMPLES[j]))
                cells[(i, j)] = v
                print(f"    {LABELS[i]:>14} vs {LABELS[j]:<14} {v:9.5f}", flush=True)
        out[method] = cells
    return out


def check_workbook(mats):
    """DATA S4 holds sqrt(MSD) on the x100 scale; bare 1/sqrt is 100/that."""
    if not XLSX.exists():
        print("  workbook not found -- skipping the cross-check")
        return
    ws = openpyxl.load_workbook(XLSX, data_only=True)[SHEET]
    bad = tot = 0
    for i in range(4):
        for j in range(4):
            v = ws.cell(row=5 + i, column=3 + j).value
            if not isinstance(v, (int, float)) or v <= 0:
                continue
            method = "CRUSH" if j > i else "Eigenvector"
            key = (min(i, j), max(i, j))
            theirs = X100 / float(v)
            mine = mats[method].get(key)
            tot += 1
            if mine is None or abs(mine - theirs) > 0.005 * max(abs(theirs), 1e-9):
                print(f"  DIFFERS  {method} {LABELS[key[0]]} vs {LABELS[key[1]]}: "
                      f"recomputed {mine!r} vs sheet {theirs:.5f}")
                bad += 1
    print(f"  workbook cross-check: {tot - bad}/{tot} cells agree")


def main() -> None:
    print("measuring with MSD_tracks.py:")
    mats = measure()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [{"method": m, "sample_a": LABELS[i], "sample_b": LABELS[j],
             "inv_sqrt_MSD_bare": v}
            for m, cells in mats.items() for (i, j), v in cells.items()]
    p = OUT / "fig1C_matrix_sourcedata.json"
    p.write_text(json.dumps(rows, indent=1))
    print(f"\nwrote {p}")
    check_workbook(mats)

    def top(cells, step):
        return float(np.ceil(max(cells.values()) / step) * step)
    ct, et = top(mats["CRUSH"], 1.0), top(mats["Eigenvector"], 0.01)
    hp.draw([("CRUSH", hp.DARKRED, ct, f"{ct:g}", mats["CRUSH"]),
             ("Eigenvector", hp.DARKBLUE, et, f"{et:g}", mats["Eigenvector"])],
            LABELS, OUT / "fig1C_matrix.png",
            fmt=lambda v: f"{v:.2f}" if v >= 0.5 else f"{v:.4f}", values=True)
    for nm, v in mats.items():
        print(f"  {nm:<12} {min(v.values()):.5f}-{max(v.values()):.5f}")
    print("  scales are independent -- darkness is NOT comparable between panels")


if __name__ == "__main__":
    main()
