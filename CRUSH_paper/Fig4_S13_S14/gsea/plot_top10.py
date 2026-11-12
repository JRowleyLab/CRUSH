#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
E_COL, M_COL = "#b5307a", "#2f6fe0"
N = 5
METHODS = {"gene": "Normalized ΔCRUSH (mammoth − elephant)",
           "state": "State score (logistic, scale 50)",
           "stateA": "State score, scale 50, A in either species",
           "stateAnarrow": "State score, scale 15, A in either species"}
LIBS = {"GO_Biological_Process_2023": "GO BP 2023", "GO_BP_full": "GO BP (full)",
        "GO_Molecular_Function_2023": "GO MF 2023", "GO_Cellular_Component_2023": "GO CC 2023",
        "MSigDB_Hallmark_2020": "MSigDB Hallmark", "KEGG_2021_Human": "KEGG 2021",
        "Reactome_Pathways_2024": "Reactome 2024", "WikiPathways_2024_Human": "WikiPathways 2024"}


def top(f):
    r = pd.read_csv(f, sep="\t").sort_values("NES")
    r = pd.concat([r[r.NES < 0].drop_duplicates("leading_edge", keep="first").head(N),
                   r[r.NES > 0].drop_duplicates("leading_edge", keep="last").tail(N)])
    return r.sort_values("NES")


def label(t, lib):
    t = t.split(" (GO:")[0] if "GO" in lib else t.split(" WP")[0] if "Wiki" in lib else t.split(" R-HSA")[0]
    return t if len(t) <= 55 else t[:53] + "…"


def main():
    for v, title in METHODS.items():
        for notext in (False, True):
            fig, axs = plt.subplots(4, 2, figsize=(7.5, 9.0))
            for ax, (lib, name) in zip(axs.flat, LIBS.items()):
                t = top(HERE / f"gsea_{v}_{lib}.tsv")
                y = np.arange(len(t))
                ax.barh(y, t.NES, color=np.where(t.NES > 0, M_COL, E_COL), height=0.7, lw=0)
                ax.axvline(0, color="black", lw=0.6)
                ax.set_yticks(y)
                ax.set_xlim(-2.4, 2.4)
                ax.tick_params(labelsize=6)
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
                if notext:
                    ax.set_yticklabels([])
                    ax.set_xticklabels([])
                    continue
                ax.set_yticklabels([f"{'† ' if c == 'yes' else ''}{label(s, lib)} ({q:.2g})" for s, q, c in
                                    zip(t.term, t.FDR, t.clustered.fillna(""))], fontsize=5.5)
                ax.set_title(name, fontsize=7, loc="left")
                ax.set_xlabel("NES", fontsize=6.5)
            if not notext:
                fig.suptitle(f"{title}   (< 0 A in elephant, > 0 A in mammoth; FDR in brackets; † clustered)",
                             fontsize=7.5)
            fig.tight_layout()
            for ext in ("png", "pdf", "svg"):
                fig.savefig(HERE / f"gsea_top10_{v}{'_notext' if notext else ''}.{ext}", dpi=300,
                            bbox_inches="tight")
            plt.close(fig)


if __name__ == "__main__":
    main()
