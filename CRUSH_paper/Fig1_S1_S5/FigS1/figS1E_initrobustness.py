#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
DATA = FIG / "init_robustness"
STEM = "figS1E_initrobustness"
RED, GREY, INK, INK2 = "#C6303A", "#b9bec6", "#16181c", "#5a5f68"


def sourcedata(res, rand, raw, ali, gc):
    """Every point the panel draws, one row each, for the SourceData workbook."""
    rows = []
    if gc is not None:
        rows.append({"group": "GC-only seeds", "arm": "gconly", "measure": "raw",
                     "agreement_pct": 100 * gc["raw"], "n_bins": gc["n_bins"],
                     "disagreeing_bins": gc["disagreeing_bins"],
                     "chromosomes_flipped": len(gc["chromosomes_flipped"])})
        rows.append({"group": "GC-only seeds", "arm": "gconly", "measure": "aligned",
                     "agreement_pct": 100 * gc["aligned"], "n_bins": gc["n_bins"],
                     "disagreeing_bins": gc["disagreeing_bins"],
                     "chromosomes_flipped": len(gc["chromosomes_flipped"])})
    for k in rand:
        r = res[k]
        for measure, val in (("raw", 100 * r["raw"]), ("aligned", 100 * r["aligned"])):
            rows.append({"group": f"Random seeds (n={len(rand)})",
                         "arm": k.split()[-1], "measure": measure,
                         "agreement_pct": val, "n_bins": r["n_bins"],
                         "disagreeing_bins": r["disagreeing_bins"],
                         "chromosomes_flipped": len(r["chromosomes_flipped"])})
    summary = [("Random seeds", "raw", raw), ("Random seeds", "aligned", ali)]
    for grp, measure, vals in summary:
        rows.append({"group": grp, "arm": "MEAN", "measure": measure,
                     "agreement_pct": vals.mean(), "n_bins": "",
                     "disagreeing_bins": "", "chromosomes_flipped": ""})
        rows.append({"group": grp, "arm": "SD", "measure": measure,
                     "agreement_pct": vals.std(ddof=1), "n_bins": "",
                     "disagreeing_bins": "", "chromosomes_flipped": ""})
    out = OUT / f"{STEM}_sourcedata.csv"
    cols = ["group", "arm", "measure", "agreement_pct", "n_bins",
            "disagreeing_bins", "chromosomes_flipped"]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")


def main() -> None:
    agree = OUT / f"{STEM}_agreement.json"
    if not agree.exists():
        print(f"{agree.name} not found -- running the analysis first\n")
        import subprocess
        subprocess.run([sys.executable,
                        str(HERE / "figS1E_initrobustness_analyze.py")], check=True)
        print()
    res = json.loads(agree.read_text())
    rand = sorted(k for k in res if k.startswith("std vs rand"))
    raw = np.array([100 * res[k]["raw"] for k in rand])
    ali = np.array([100 * res[k]["aligned"] for k in rand])
    gc = res.get("std vs gconly")
    n = len(rand)

    groups = []
    if gc is not None:
        groups.append(("GC-only seeds",
                       np.array([100 * gc["raw"]]), np.array([100 * gc["aligned"]])))
    groups.append((f"Random seeds (n={n})", raw, ali))

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    fig.patch.set_facecolor("white")
    rng = np.random.default_rng(0)
    w = 0.34
    for g, (glab, graw, gali) in enumerate(groups):
        for k, (vals, col) in enumerate(((graw, GREY), (gali, RED))):
            x = g + (k - 0.5) * (w + 0.04)
            ax.bar(x, vals.mean(), w, facecolor=col, edgecolor="#5a5f68",
                   linewidth=1.1, zorder=2)
            if len(vals) > 1:
                sd = vals.std(ddof=1)
                ax.errorbar(x, vals.mean(), yerr=sd, fmt="none", ecolor=INK2,
                            elinewidth=2.0, capsize=7, capthick=2.0, zorder=4)
                ax.scatter(x + rng.uniform(-0.11, 0.11, len(vals)), vals, s=50,
                           facecolor="white", edgecolor=INK, linewidth=1.6,
                           zorder=5, clip_on=False)
            m = vals.mean()
            ax.text(x, 103.0, f"{m:.1f}" if m < 99.9 else f"{m:.2f}",
                    ha="center", fontsize=10, color=INK2 if col == GREY else RED,
                    fontweight="bold")
            ax.text(x, -4.0, "raw" if col == GREY else "aligned", ha="center",
                    va="top", fontsize=9.5, color=INK2, clip_on=False)
        ax.text(g, -12.5, glab, ha="center", va="top", fontsize=10.5,
                color=INK, clip_on=False)

    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_ylim(0, 108)
    ax.set_xticks([])

    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("A/B call agreement with standard\n(per 1 kb bin, %)",
                  fontsize=10.5, color=INK, labelpad=8)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#d7d7d3")
    ax.tick_params(colors=INK2, labelsize=10, length=0)
    ax.grid(True, axis="y", color="#ececE8", linewidth=0.8)
    ax.set_axisbelow(True)
    nb = res[rand[0]]["n_bins"]
    ax.text(0.0, 1.04, f"{n} independent random initialisations, {nb:,} bins",
            transform=ax.transAxes, fontsize=9.5, color=INK2)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(bottom=0.18, top=0.88, left=0.16, right=0.97)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{STEM}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT / f'{STEM}.png'} (+ .pdf)\n")
    sourcedata(res, rand, raw, ali, gc)
    print(f"  {'arm':<10}{'raw':>9}{'aligned':>10}{'disagree':>11}{'chr flipped':>13}")
    for k in rand:
        r = res[k]
        print(f"  {k.split()[-1]:<10}{100*r['raw']:>8.1f}%{100*r['aligned']:>9.3f}%"
              f"{r['disagreeing_bins']:>11,}{len(r['chromosomes_flipped']):>13}")
    print(f"\n  raw     mean {raw.mean():.2f}%  SD {raw.std(ddof=1) if n>1 else 0:.2f}")
    print(f"  aligned mean {ali.mean():.3f}%  SD {ali.std(ddof=1) if n>1 else 0:.3f}")


if __name__ == "__main__":
    main()
