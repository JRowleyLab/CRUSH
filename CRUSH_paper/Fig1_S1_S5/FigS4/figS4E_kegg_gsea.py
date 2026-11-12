#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results"
WINDOW = sys.argv[1] if len(sys.argv) > 1 else "body"
LIB = OUT / "KEGG_2021_Human.gmt"


def main() -> None:
    sc = pd.read_csv(OUT / f"figS4_genescores_{WINDOW}.tsv", sep="\t", index_col=0)
    sc.index = [g.upper() for g in sc.index]
    sc = sc[~sc.index.duplicated()]
    delta = (sc[["Erythro", "eRBC", "aRBC"]].mean(axis=1)
             - sc["Fibro"]).sort_values(ascending=False)
    rnk = OUT / f"figS4E_ranks_{WINDOW}.rnk"
    delta.to_csv(rnk, sep="\t", header=False)
    print(f"  {len(delta):,} genes, delta {delta.min():+.1f} to {delta.max():+.1f}")

    gsea = OUT / f"figS4E_kegg_gsea_{WINDOW}.tsv"
    subprocess.run(["Rscript", str(HERE / "figS4_gsea.R"), str(rnk), str(LIB),
                    str(gsea)], check=True)
    rot = OUT / f"figS4E_kegg_rotation_{WINDOW}.tsv"
    subprocess.run([sys.executable, str(HERE / "figS4_rotation.py"), str(rnk),
                    str(gsea), str(LIB), str(rot)], check=True)

    d = pd.read_csv(rot, sep="\t")
    print(f"\n  {'term':<48}{'NES':>8}{'n':>5}{'p(rot)':>9}")
    for _, r in d.iterrows():
        print(f"  {r['term'][:48]:<48}{r['NES']:>8.2f}"
              f"{r['n_genes']:>5}{r['p_rotation']:>9.3f}")


if __name__ == "__main__":
    main()
