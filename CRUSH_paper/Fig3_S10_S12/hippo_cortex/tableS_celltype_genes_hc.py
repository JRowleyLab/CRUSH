#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

import fig3F_signatures_hc as f

HERE = Path(__file__).resolve().parent
FDR = 0.05
SHOWN = {"Kcnh1": "Fig. 3G", "St18": "Fig. 3G", "Dock10": "Fig. 3G", "Hhat": "Fig. S11C",
         "Brinp3": "Fig. S11C", "Slc40a1": "Fig. S11C", "Prex2": "Fig. S11C", "Gli2": "Fig. S11C",
         "Sntg1": "Fig. S11C"}
SHEET = {"Neuron": "Neuron", "Mature Oligodendrocyte": "Mature oligodendrocyte",
         "Microglia Etc.": "Microglia", "Adult Astrocyte": "Adult astrocyte",
         "Neonatal Astrocyte": "Neonatal astrocyte", "Oligodendrocyte Progenitor": "OPC"}


def main():
    U = pd.read_csv(f.UMAP)
    G = f.gene_matrix()
    typed = (U.Group != "Unknown").to_numpy()
    Gt, grp = G.loc[:, typed], U.Group.to_numpy()[typed]
    Z = G.sub(G.mean(axis=1), axis=0).div(G.std(axis=1, ddof=0), axis=0).loc[:, typed]
    bed = pd.read_csv(f.GENES, sep="\t", header=None, usecols=[0, 1, 2, 3, 5], dtype={0: str})
    bed = bed[bed[0].isin(["1", "chr1"])].drop_duplicates(3).set_index(3)
    sig = pd.read_csv(HERE / "fig3F_hc_signatures.tsv", sep="\t")

    out = []
    for t in f.TYPES:
        y = (grp == t).astype(int)
        rows = []
        for gene, v in zip(Gt.index, Gt.to_numpy()):
            if v.std() == 0:
                continue
            z = Z.loc[gene].to_numpy()
            rows.append({"gene": gene, "AUC": roc_auc_score(y, v),
                         "median_z_type": np.median(z[y == 1]), "median_z_other": np.median(z[y == 0]),
                         "p": mannwhitneyu(v[y == 1], v[y == 0]).pvalue})
        r = pd.DataFrame(rows)
        r["spec"] = 2 * r.AUC - 1
        r["FDR"] = f.bh(r.p)
        r = r[r.FDR < FDR].copy()
        r["direction"] = np.where(r.spec > 0, "A in type", "B in type")
        r["top150_signature"] = r.gene.isin(sig.loc[sig.cell_type == t, "gene"]).map({True: "yes", False: ""})
        r["shown_in"] = [SHOWN[g] if g in SHOWN and GENES_TYPE.get(g) == t else "" for g in r.gene]
        r.insert(0, "cell_type", SHEET[t])
        r.insert(1, "n_cells_type", int(y.sum()))
        r.insert(2, "n_cells_other", int((1 - y).sum()))
        b = bed.loc[r.gene]
        r.insert(4, "chrom", "chr1")
        r.insert(5, "start", b[1].to_numpy())
        r.insert(6, "end", b[2].to_numpy())
        r.insert(7, "strand", b[5].to_numpy())
        out.append(r.sort_values(["direction", "spec"], key=lambda s: s.abs() if s.name == "spec" else s,
                                 ascending=[True, False]))
        print(t, (r.direction == "A in type").sum(), "A,", (r.direction == "B in type").sum(), "B")
    T = pd.concat(out)[["cell_type", "n_cells_type", "n_cells_other", "gene", "chrom", "start", "end",
                         "strand", "direction", "AUC", "spec", "median_z_type", "median_z_other", "p",
                         "FDR", "top150_signature", "shown_in"]]
    T.to_csv(HERE / "tableS_celltype_genes_hc.tsv", sep="\t", index=False, float_format="%.4g")

    xl = HERE / "tableS_celltype_genes_hc.xlsx"
    readme = [
        ["Genes with cell-type-specific compartment patterns, mouse hippocampus + cortex Dip-C (Tan et al. 2021)"],
        [""],
        ["Data", "CRUSH v1 single-cell scores, 25 kb bins, chr1; gene value = mean over the 25 kb bins the gene "
                 "overlaps (mm10 annotation); missing bins filled with the gene's mean across cells"],
        ["Comparison", "Cells of one type vs all other typed cells (neuronal subtypes merged; Unknown cells excluded)"],
        ["Kept", f"Genes with FDR < {FDR}"],
        [""],
        ["Column", "Description"],
        ["cell_type", "Cell type tested"],
        ["n_cells_type / n_cells_other", "Number of cells of that type / of all other typed cells"],
        ["gene, chrom, start, end, strand", "Gene and mm10 coordinates"],
        ["direction", "A in type: higher CRUSH (more A) in that cell type; B in type: lower CRUSH (more B)"],
        ["AUC", "Area under the ROC curve for separating that cell type from the rest by the gene's CRUSH score"],
        ["spec", "2 x AUC - 1 (-1 to 1)"],
        ["median_z_type / median_z_other", "Median CRUSH z-score (z-scored across all cells) in that type / the rest"],
        ["p", "Two-sided Wilcoxon rank-sum (Mann-Whitney U) test"],
        ["FDR", "Benjamini-Hochberg, within each cell type over all tested chr1 genes"],
        ["top150_signature", "Gene is in the top-150 A signature used for the Fig. 3F ROC"],
        ["shown_in", "Figure panel showing this gene"],
    ]
    with pd.ExcelWriter(xl, engine="openpyxl") as w:
        pd.DataFrame(readme).to_excel(w, sheet_name="README", index=False, header=False)
        for t in f.TYPES:
            T[T.cell_type == SHEET[t]].to_excel(w, sheet_name=SHEET[t], index=False)
    wb = load_workbook(xl)
    for ws in wb:
        for row in ws.iter_rows():
            for c in row:
                c.font = Font(name="Arial", size=10, bold=c.font.bold)
                if isinstance(c.value, float):
                    c.number_format = "0.00E+00" if ws.cell(1, c.column).value in ("p", "FDR") else "0.000"
        ws.freeze_panes = "A2" if ws.title != "README" else None
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = (
                max(len(str(c.value or "")) for c in col) + 2 if ws.title != "README" or col[0].column > 1 else 32)
    ws = wb["README"]
    ws["A1"].font = Font(name="Arial", size=11, bold=True)
    ws["A7"].font = ws["B7"].font = Font(name="Arial", size=10, bold=True)
    ws.column_dimensions["B"].width = 100
    wb.save(xl)


GENES_TYPE = {"Kcnh1": "Neuron", "St18": "Mature Oligodendrocyte", "Dock10": "Microglia Etc.",
              "Hhat": "Neuron", "Brinp3": "Mature Oligodendrocyte", "Slc40a1": "Microglia Etc.",
              "Prex2": "Adult Astrocyte", "Gli2": "Neonatal Astrocyte", "Sntg1": "Oligodendrocyte Progenitor"}


if __name__ == "__main__":
    main()
