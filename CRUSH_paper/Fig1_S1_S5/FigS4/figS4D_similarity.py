#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import csv
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
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
TAB = "DATA_S7_Fig.S2D"
FIGPREP = Path("chicken")
PC1DIR = Path("hicPCA")

EIGEN = "pc1"
SHARED_SCALE = True
SHARED_VMAX = 0.03
X100 = 100.0
SAMPLES = ["Fibroblast", "Erythroblast", "eRBC", "aRBC"]
LABELS = ["Fibro", "Erythro", "eRBC", "aRBC"]
SHORT = dict(zip(SAMPLES, LABELS))
ARMS = [("eigen_50kb", "Eigenvector 50 kb"),
        ("crush_50kb", "CRUSH 50 kb"),
        ("crush_1kb", "CRUSH 1 kb")]


def path_for(arm, sample):
    if arm == "eigen_50kb":
        if EIGEN == "pc1":
            return PC1DIR / f"{SHORT[sample]}1.bw"
        return FIGPREP / f"{SHORT[sample]}_Eigen_50kb_flipped.bedgraph"
    res = 50000 if arm == "crush_50kb" else 1000
    return REV / "tracks" / "chicken" / f"{sample}_r{res}" / \
        f"{sample}_mergedCrush_{res}.bedgraph"


_cache: dict[Path, object] = {}


def track(p: Path):
    if p not in _cache:
        if not p.exists():
            raise SystemExit(f"missing track: {p}")
        _cache[p] = msd.keyed(msd.load(str(p)), literal=False)
    return _cache[p]


def similarity(a: Path, b: Path) -> float:
    """1/sqrt(MSD) on the IQR x100 scale, as DATA S7 stores it."""
    tag = f"{a.name} vs {b.name}"
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        m = msd.align(track(a), track(b), False, tag)
        m = msd.orient_to_control(m, False, 100, tag)
        xi, yi = msd.track_iqr(m)
    if xi == 0 or yi == 0:
        raise SystemExit(f"zero IQR for {tag}")
    d = ((m["Score_x"] / xi) * X100 - (m["Score_y"] / yi) * X100) ** 2
    s = float(np.sqrt(d.mean()))
    return 1.0 / s if s else float("nan")


def measure():
    out = {}
    for arm, title in ARMS:
        cells = {}
        print(f"  {title}:")
        for i in range(len(SAMPLES)):
            for j in range(i + 1, len(SAMPLES)):
                v = similarity(path_for(arm, SAMPLES[i]), path_for(arm, SAMPLES[j]))
                cells[(i, j)] = v
                print(f"    {LABELS[i]:>8} vs {LABELS[j]:<8} {v:9.5f}", flush=True)
        out[arm] = cells
    return out


def check_workbook(mats):
    """DATA S7: three stacked staggered triangles, each under its own header."""
    if not XLSX.exists():
        print("  workbook not found -- skipping the cross-check")
        return
    ws = openpyxl.load_workbook(XLSX, data_only=True)[TAB]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    
    heads = {"eigen_50kb": ("eigen 50",), "crush_50kb": ("crush 50", "crush50"),
             "crush_1kb": ("crush 1", "crush1")}
    bad = tot = 0
    for arm, _title in ARMS:
        keys = heads[arm]
        hdr = next((k for k, r in enumerate(rows)
                    if r and r[0]
                    and any(str(r[0]).strip().lower().startswith(x) for x in keys)),
                   None)
        if hdr is None:
            print(f"  no block for {arm} in {TAB}")
            continue
        for i in range(3):
            r = rows[hdr + 2 + i]
            for j in range(i + 1, 4):
                v = r[j] if j < len(r) else None
                if not isinstance(v, (int, float)):
                    continue
                tot += 1
                mine = mats[arm][(i, j)]
                if abs(mine - float(v)) > 0.005 * max(abs(float(v)), 1e-9):
                    print(f"  DIFFERS  {arm} {LABELS[i]} vs {LABELS[j]}: "
                          f"recomputed {mine:.5f} vs sheet {float(v):.5f}")
                    bad += 1
    print(f"  workbook cross-check: {tot - bad}/{tot} cells agree")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"measuring with MSD_tracks.py (eigenvector = {EIGEN}):")
    mats = measure()

    q = OUT / "figS4D_similarity_sourcedata.csv"
    with q.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["arm", "cell_a", "cell_b", "inv_sqrt_MSD_x100"])
        for arm, title in ARMS:
            for (i, j), v in mats[arm].items():
                w.writerow([title, LABELS[i], LABELS[j], f"{v:.6f}"])
    print(f"\nwrote {q}")
    check_workbook(mats)

    def top(cells, step):
        return float(np.ceil(max(cells.values()) / step) * step)
    if SHARED_SCALE:
        allv = [v for arm, _ in ARMS for v in mats[arm].values()]
        vmax = (SHARED_VMAX if SHARED_VMAX is not None
                else float(np.ceil(max(allv) / 0.005) * 0.005))
        over = [v for v in allv if v > vmax]
        if over:
            print(f"  {len(over)} of {len(allv)} cells saturate at {vmax:g} "
                  f"(max {max(allv):.5f})")
        panels = [(title, hp.DARKRED, vmax, f"{vmax:g}", mats[arm])
                  for arm, title in ARMS]
        print(f"  shared scale 0-{vmax:g} across all three panels")
    else:
        panels = [(title, hp.DARKBLUE if arm == "eigen_50kb" else hp.DARKRED,
                   top(mats[arm], 0.005), f"{top(mats[arm], 0.005):g}", mats[arm])
                  for arm, title in ARMS]
    for labels in (True, False):
        stem = "figS4D_similarity" + ("_labeled" if labels else "")
        hp.draw(panels, LABELS, OUT / f"{stem}.png",
                fmt=lambda v: f"{v:.3f}", values=labels)
        print(f"wrote {OUT / (stem + '.png')}")
    for arm, title in ARMS:
        v = mats[arm]
        print(f"  {title:<20} {min(v.values()):.5f}-{max(v.values()):.5f}")


if __name__ == "__main__":
    main()
