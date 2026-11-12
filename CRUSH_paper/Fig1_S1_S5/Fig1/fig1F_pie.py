#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
TRACKS = Path("tnbc")
LINES = {"hmec": "HMEC", "t47D": "T47D", "zr7530": "ZR7530",
         "hcc1954": "HCC1954", "hcc70": "HCC70", "bt549": "BT549"}

THRESH = 50.0
INK, INK2 = "#16181c", "#5a5f68"
SLICE = ["#8B1F24", "#C6303A", "#E08A87", "#F3D3CF"]
NAMES = ["Invariant", "Differs in 1", "Differs in 2", "Differs in 3"]


def load(line: str) -> pd.DataFrame:
    p = TRACKS / f"{line}_r1000" / f"{line}_mergedCrush_1000.bedgraph"
    if not p.exists():
        raise SystemExit(f"missing {p}")
    df = pd.read_csv(p, sep="\t", header=None, skiprows=1,
                     names=["chrom", "start", "end", "score"],
                     dtype={"chrom": str})
    if df["score"].isna().any():
        raise SystemExit(f"unparsed rows in {p}")
    df = df[df["score"] != 0]
    df = df[df["score"].abs() >= THRESH]
    return pd.DataFrame({"key": df["chrom"] + "_" + df["start"].astype(str),
                         line: (df["score"].to_numpy() > 0).astype(np.int8)})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    calls = None
    for line, lab in LINES.items():
        s = load(line)
        print(f"  {lab:<9}{len(s):>10,} bins", flush=True)
        calls = s if calls is None else calls.merge(s, on="key")
    print(f"  {'shared':<9}{len(calls):>10,} bins\n", flush=True)

    mat = calls[list(LINES)].to_numpy()
    n_a = mat.sum(axis=1)
    minority = np.minimum(n_a, 6 - n_a)
    counts = Counter(minority.tolist())
    total = len(minority)
    vals = [counts.get(k, 0) for k in range(4)]
    frac = [v / total for v in vals]

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    fig.patch.set_facecolor("white")
    wedges, _, autotexts = ax.pie(
        vals, colors=SLICE, startangle=90, counterclock=False,
        autopct=lambda p: f"{p:.1f}%", pctdistance=0.72,
        wedgeprops={"edgecolor": "white", "linewidth": 2.5},
        textprops={"fontsize": 12.5, "fontweight": "bold"})
    for t, c, w, f in zip(autotexts, SLICE, wedges, frac):
        r, g, b = (int(c[i:i + 2], 16) / 255 for i in (1, 3, 5))
        lum = .2126 * r + .7152 * g + .0722 * b
        t.set_color("white" if lum < 0.55 else INK)
        if f < 0.08:
            ang = np.deg2rad((w.theta1 + w.theta2) / 2)
            t.set_position((1.16 * np.cos(ang), 1.16 * np.sin(ang)))
            t.set_color(INK)
            t.set_fontsize(11)
    ax.set_aspect("equal")

    fig.tight_layout()
    fig.savefig(OUT / "fig1F_pie.png", dpi=300, facecolor="white")
    (OUT / "fig1F_pie.json").write_text(json.dumps(
        {"n_bins": total, "lines": list(LINES.values()),
         "counts": dict(zip(NAMES, vals)),
         "fraction": dict(zip(NAMES, frac)),
         "threshold": THRESH,
         "source": "revision/tracks/tnbc, crush-rs v0.2.0, 1 kb, |score|>=50"}, indent=1))

    print(f"{'category':<16}{'bins':>12}{'share':>9}")
    for n, v, f in zip(NAMES, vals, frac):
        print(f"{n:<16}{v:>12,}{100*f:>8.2f}%")
    print(f"\nwrote {OUT / 'fig1F_pie.png'}")


if __name__ == "__main__":
    main()
