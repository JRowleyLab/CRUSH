#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, NullFormatter

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
PAC = FIG.parent / "reviewer" / "pachytene"
RED, BLUE, GREY = "#C6303A", "#2E77B5", "#8a8f98"
INK, INK2 = "#16181c", "#5a5f68"
LCL = json.loads((FIG.parent / "reviewer" / "validation" / "mindist"
                  / "lcl_gw_agreement.json").read_text())
LCL = {int(k): 100 * v["aligned"] for k, v in LCL.items()}


def main() -> None:
    d = json.loads((PAC / "mindist_agreement.json").read_text())
    fig, ax = plt.subplots(figsize=(6.8, 4.7))
    fig.patch.set_facecolor("white")

    series = [("pachytene", RED, "-", d.get("pachytene", {})),
              ("fibroblast (somatic)", BLUE, "-", d.get("fibroblast (somatic)", {}))]
    for lab, col, ls, arm in series:
        if not arm:
            continue
        x = np.array(sorted(int(k) for k in arm))
        y = np.array([100 * arm[str(k)]["aligned"] for k in x])
        ax.plot(x / 1e6, y, color=col, linestyle=ls, linewidth=5.0, marker="o",
                markersize=10, markerfacecolor=col, markeredgecolor="white",
                markeredgewidth=2.0, solid_capstyle="round", clip_on=False,
                zorder=3, label=lab)
    lx = np.array(sorted(LCL))
    ax.plot(lx / 1e6, [LCL[k] for k in lx], color=GREY, linestyle=(0, (5, 2.5)),
            linewidth=3.0, marker="s", markersize=7, markerfacecolor=GREY,
            markeredgecolor="white", markeredgewidth=1.5, zorder=2,
            label="LCL, human (somatic, genome-wide)")

    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([0.1, 0.5, 1, 2, 5, 10]))
    ax.set_xticklabels(["0.1", "0.5", "1", "2", "5", "10"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(70, 101)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d7d7d3")
    ax.tick_params(colors=INK2, labelsize=10, length=0)
    ax.grid(True, color="#ececE8", linewidth=0.8)
    ax.set_axisbelow(True)


    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "figS2E_pachytene_mindist.png"
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, facecolor="white")
    print(f"wrote {out} (+ .pdf)\n")

    import csv
    q = OUT / "figS2E_pachytene_mindist_sourcedata.csv"
    with q.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["series", "min_distance_Mb", "percent_agreement_aligned",
                    "percent_agreement_raw", "n_bins", "disagreeing_bins"])
        for lab, _, _, arm in series:
            for k in sorted(arm, key=int):
                v = arm[k]
                w.writerow([lab, f"{int(k)/1e6:g}", f"{100*v['aligned']:.4f}",
                            f"{100*v['raw']:.4f}", v["n"], v["disagree"]])
        for k in sorted(LCL):
            w.writerow(["LCL, human (somatic, genome-wide)", f"{k/1e6:g}",
                        f"{LCL[k]:.4f}", "", "", ""])
    print(f"wrote {q}")
    for lab, _, _, arm in series:
        if arm:
            print(f"  {lab:<24}" + "  ".join(
                f"{int(k)/1e6:g}M={100*arm[k]['aligned']:.1f}%" for k in sorted(arm, key=int)))


if __name__ == "__main__":
    main()
