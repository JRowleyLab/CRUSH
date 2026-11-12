#!/usr/bin/env python3
from __future__ import annotations

import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
DATA = FIG.parent / "reviewer" / "validation"
CMAP = LinearSegmentedColormap.from_list("oe", ["#2E77B5", "#ffffff", "#C6303A"])
RED, BLUE, INK, INK2 = "#C6303A", "#2E77B5", "#16181c", "#5a5f68"
NDISP = 1000


def block_mean(M, target=NDISP):
    n = M.shape[0]; k = max(1, n // target); m = (n // k) * k
    return np.nanmean(M[:m, :m].reshape(m // k, k, m // k, k), axis=(1, 3))


def down(v, target=NDISP):
    n = len(v); k = max(1, n // target); m = (n // k) * k
    return v[:m].reshape(-1, k).mean(1)


def strength(OE, sc):
    o = np.argsort(sc); M = OE[np.ix_(o, o)]
    q = len(sc) // 4
    aa = np.nanmean(M[-q:, -q:]); bb = np.nanmean(M[:q, :q]); ab = np.nanmean(M[:q, -q:])
    return M, o, (aa + bb) / (2 * ab)


def main() -> None:
    f = sorted(glob.glob(str(DATA / "win_chr*.npz")))[0]
    z = np.load(f)
    OE, cr, ei = z["OE"], z["crush"], z["eigen"]
    chrom, start, res = str(z["chrom"]), int(z["start"]), int(z["res"])
    Mc, oc, sc_ = strength(OE, cr)
    Me, oe_, se_ = strength(OE, ei)

    disp = [block_mean(OE), block_mean(Mc), block_mean(Me)]
    VMIN, VCENTER, VMAX = 0.0, 1.0, 1.5
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=VCENTER, vmax=VMAX)
    clipped = float(np.mean([np.nanmean(d > VMAX) for d in disp]))

    fig = plt.figure(figsize=(12.6, 5.4))
    fig.patch.set_facecolor("white")
    gs = GridSpec(3, 3, height_ratios=[1, 0.075, 0.075], hspace=0.10, wspace=0.07,
                  left=0.035, right=0.90, top=0.90, bottom=0.10)
    end = start + len(cr) * res
    titles = [f"genomic order\nchr{chrom}:{start/1e6:.1f}-{end/1e6:.1f} Mb",
              f"sorted by CRUSH\n(AA+BB)/2AB = {sc_:.2f}",
              f"sorted by eigenvector\n(AA+BB)/2AB = {se_:.2f}"]
    strips = [[(cr, RED), (ei, BLUE)],
              [(cr[oc], RED), (ei[oc], BLUE)],
              [(cr[oe_], RED), (ei[oe_], BLUE)]]

    for j in range(3):
        ax = fig.add_subplot(gs[0, j])
        ax.imshow(disp[j], cmap=CMAP, norm=norm, interpolation="nearest")
        ax.set_title(titles[j], fontsize=10.5, color=INK)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color("#d7d7d3")
        for r, (vals, col) in enumerate(strips[j]):
            tx = fig.add_subplot(gs[1 + r, j])
            y = down(vals)
            y = y / (np.nanmax(np.abs(y)) or 1)
            tx.fill_between(np.arange(len(y)), 0, y, where=y > 0, color=RED, lw=0)
            tx.fill_between(np.arange(len(y)), 0, y, where=y <= 0, color=BLUE, lw=0)
            tx.set_xlim(0, len(y)); tx.set_ylim(-1.05, 1.05)
            tx.set_xticks([]); tx.set_yticks([])
            for s in tx.spines.values():
                s.set_visible(False)
            if j == 0:
                tx.text(-0.012, 0.5, "CRUSH" if r == 0 else "Eigen",
                        transform=tx.transAxes, ha="right", va="center",
                        fontsize=8.5, color=col)

    cax = fig.add_axes([0.915, 0.42, 0.012, 0.40])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                      ticks=[VMIN, 0.5, VCENTER, 1.25, VMAX])
    cb.ax.tick_params(labelsize=9, colors=INK2, length=0)
    cb.outline.set_edgecolor("#8a8f98")
    cax.set_title("O/E", fontsize=9.5, color=INK, pad=8)

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "figS2A_sorted.png"
    for ext in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, facecolor="white")
    import json
    (OUT / "figS2A_sorted_sourcedata.json").write_text(json.dumps(
        {"chrom": chrom, "start": start, "res": res, "n_bins": int(len(cr)),
         "crush_strength": float(sc_), "eigen_strength": float(se_)}, indent=1))
    print(f"wrote {OUT/'figS2A_sorted_sourcedata.json'}")
    print(f"wrote {out}")
    print(f"  {len(cr):,} bins called by both;  CRUSH {sc_:.3f}  "
          f"eigenvector {se_:.3f}  ({100*(sc_/se_-1):+.1f}%)")
    print(f"  scale 0-1.5 centred on 1.0; {100*clipped:.1f}% of displayed cells "
          f"exceed 1.5 and saturate")


if __name__ == "__main__":
    main()
