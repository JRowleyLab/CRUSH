#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Ellipse, Rectangle

HERE = Path(__file__).resolve().parent
FIG = HERE.parent
OUT = FIG / "results"
CACHE = OUT / "tcga_expression.json"

RAW = "--raw" in sys.argv
COLS = ["Normal", "Luminal", "HER2+", "TNBC"]
ROWS = [("CRUSH 1 kb", "crush"), ("Eigen 100 kb", "eigen")]
VLIM, BOXED = 1.0, "TNBC"
BLUE, RED = "#2E77B5", "#C6303A"
INK, INK2 = "#16181c", "#5a5f68"
CMAP = LinearSegmentedColormap.from_list("crush_div", [BLUE, "#ffffff", RED])
PAM50 = {"BRCA_LumA": "Luminal", "BRCA_LumB": "Luminal",
         "BRCA_Her2": "HER2+", "BRCA_Basal": "TNBC"}


def fetch_expression(genes) -> dict:
    if CACHE.exists():
        d = json.loads(CACHE.read_text())
        if all(g in d for g in genes):
            return d
    api = "https://www.cbioportal.org/api"
    req = urllib.request.Request(
        f"{api}/studies/brca_tcga_pan_can_atlas_2018/clinical-data/fetch"
        "?clinicalDataType=PATIENT",
        data=json.dumps({"attributeIds": ["SUBTYPE"]}).encode(),
        headers={"Content-Type": "application/json"})
    pam = {c["patientId"]: PAM50.get(c["value"])
           for c in json.load(urllib.request.urlopen(req, timeout=180))
           if c.get("clinicalAttributeId") == "SUBTYPE"}
    print(f"  PAM50 assignments: {len(pam)} patients", flush=True)

    out = {}
    for g in genes:
        vals = {c: [] for c in COLS}
        for st in ("NT", "TP"):
            url = ("https://firebrowse.org/api/v1/Samples/mRNASeq?format=json"
                   f"&gene={g}&cohort=BRCA&sample_type={st}&protocol=RSEM"
                   "&page_size=3000")
            recs = json.load(urllib.request.urlopen(url, timeout=240)).get("mRNASeq", [])
            for r in recs:
                v = r.get("expression_log2")
                if v is None:
                    continue
                if st == "NT":
                    vals["Normal"].append(float(v))
                else:
                    grp = pam.get(r["tcga_participant_barcode"])
                    if grp:
                        vals[grp].append(float(v))
        out[g] = vals
        print(f"  {g:<8}" + "  ".join(f"{c}={len(vals[c])}" for c in COLS), flush=True)
    CACHE.write_text(json.dumps(out, indent=1))
    return out


def panel(gene: str, expr: dict) -> None:
    tag = "_promoter_raw" if RAW else "_promoter"
    d = json.loads((OUT / f"fig1G_{gene}{tag}_subtype.json").read_text())
    fig = plt.figure(figsize=(5.4, 5.0))
    fig.patch.set_facecolor("white")
    L, W = 0.245, 0.575
    ax = fig.add_axes([L, 0.665, W, 0.245])
    bx = fig.add_axes([L, 0.215, W, 0.395])

    ax.set_xlim(-0.5, len(COLS) - 0.5)
    ax.set_ylim(len(ROWS) - 0.5, -0.5)
    ax.set_axis_off()
    norm = Normalize(-VLIM, VLIM)
    for r, (rlab, key) in enumerate(ROWS):
        for c, col in enumerate(COLS):
            ax.add_patch(Ellipse((c, r), 0.62, 0.52,
                                 facecolor=CMAP(norm(np.clip(d[key][col], -VLIM, VLIM))),
                                 edgecolor="#8a8f98", linewidth=1.0, zorder=3))
        ax.text(-0.60, r, rlab, ha="right", va="center", fontsize=10, color=INK)

    vals = [expr[gene][c] for c in COLS]
    bp = bx.boxplot(vals, positions=range(len(COLS)), widths=0.52,
                    patch_artist=True, showfliers=False, whis=(5, 95))
    for i, box in enumerate(bp["boxes"]):
        box.set(facecolor=RED if COLS[i] == BOXED else "#c8cdd4",
                edgecolor="#5a5f68", linewidth=1.1, alpha=0.95)
    for part in ("whiskers", "caps"):
        for a in bp[part]:
            a.set(color="#5a5f68", linewidth=1.1)
    for m in bp["medians"]:
        m.set(color=INK, linewidth=1.8)
    bx.set_xlim(-0.5, len(COLS) - 0.5)
    bx.set_ylabel("Expression\n(log2 RSEM)", fontsize=9.5, color=INK, labelpad=6)
    bx.set_xticks(range(len(COLS)))
    bx.set_xticklabels(COLS, rotation=32, ha="right", fontsize=10, color=INK)
    for s in ("top", "right"):
        bx.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        bx.spines[s].set_color("#d7d7d3")
    bx.tick_params(colors=INK2, labelsize=9, length=0)
    bx.grid(True, axis="y", color="#ececE8", linewidth=0.8)
    bx.set_axisbelow(True)

    i = COLS.index(BOXED)
    x0 = L + W * (i - 0.42 + 0.5) / len(COLS)
    fig.patches.append(Rectangle((x0, 0.215), W * 0.84 / len(COLS), 0.695,
                                 transform=fig.transFigure, fill=False,
                                 edgecolor=INK, linewidth=2.0,
                                 linestyle=(0, (4, 2.4)), zorder=10))

    cax = fig.add_axes([0.855, 0.675, 0.030, 0.225])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                      ticks=[-VLIM, VLIM])
    cb.ax.set_yticklabels([f"-{VLIM:g}", f"{VLIM:g}"])
    cb.outline.set_edgecolor("#8a8f98")
    cb.ax.tick_params(labelsize=9, colors=INK2, length=0)
    cax.text(2.5, 0.5, "score / IQR" if RAW else "ΔIQR\n(v.s. HMEC)",
             transform=cax.transAxes, ha="center",
             va="center", rotation=90, fontsize=9, color=INK)
    fig.text(L, 0.955, gene, fontsize=13, fontweight="bold", color=INK)

    out = OUT / f"fig1G_{gene}_full{'_raw' if RAW else ''}.png"
    fig.savefig(out, dpi=300, facecolor="white")
    plt.close(fig)
    med = {c: float(np.median(expr[gene][c])) for c in COLS}
    print(f"\n{gene}  {d['locus']}  [{d.get('window','')}]")
    print(f"  {'subtype':<9}{'CRUSH':>8}{'Eigen':>8}{'log2 RSEM':>11}{'n':>6}")
    for c in COLS:
        print(f"  {c:<9}{d['crush'][c]:>8.2f}{d['eigen'][c]:>8.3f}"
              f"{med[c]:>11.2f}{len(expr[gene][c]):>6}")
    print(f"  -> {out.name}")


if __name__ == "__main__":
    genes = [a for a in sys.argv[1:] if not a.startswith("-")] or ["RSU1", "STAM"]
    print("fetching TCGA expression ...", flush=True)
    expr = fetch_expression(genes)
    for g in genes:
        panel(g, expr)
