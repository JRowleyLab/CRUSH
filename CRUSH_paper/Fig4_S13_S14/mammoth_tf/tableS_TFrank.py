#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font

from fig4E_targets_thr import THR, adjusted

HERE = Path(__file__).resolve().parent
NEW4E = ["GMEB2", "BHLHE40", "ADNP", "OTX2", "FOXR2", "POU6F2"]


def main():
    r = pd.read_csv(HERE / "fig4E_rank.tsv", sep="\t")
    d, reg = adjusted()
    T = pd.DataFrame({
        "rank": r["rank"], "TF": r.TF, "direction": np.where(r.combined < 0, "A in elephant", "A in mammoth"),
        "elephant_CRUSH": r.elephant_CRUSH, "mammoth_CRUSH": r.mammoth_CRUSH,
        "TF_dCRUSH_norm": r.TF_dCRUSH_norm,
        "compartment_switch": np.where(np.sign(r.elephant_CRUSH) != np.sign(r.mammoth_CRUSH), "yes", ""),
        "n_targets": r.n_targets, "target_r": r.rank_biserial, "target_p": r.p_rotation,
        "target_FDR": r.FDR_rotation, "TF_score": r.TF_score, "target_score": r.target_score,
        "combined_score": r.combined,
        "shown_in": np.where(r.TF.isin(NEW4E), "Fig. 4E", "")})
    for T_ in THR:
        conc, disc = [], []
        for tf, c in zip(r.TF, r.combined):
            x = np.sign(c) * d[d.index.isin(reg[tf] - {tf})]
            conc.append(int((x > 0).sum() if T_ == 0 else (x >= T_).sum()))
            disc.append(int((x < 0).sum() if T_ == 0 else (x <= -T_).sum()))
        T[f"concordant_T{T_}"], T[f"discordant_T{T_}"] = conc, disc
    assert (T.concordant_T0 + T.discordant_T0 <= T.n_targets).all()
    T.to_csv(HERE / "tableS_TFrank.tsv", sep="\t", index=False, float_format="%.4g")
    print(T[T.shown_in != ""].to_string(index=False))

    xl = HERE / "tableS_TFrank.xlsx"
    readme = [
        ["Transcription factors ranked by the compartment shift of the TF gene and its ChEA3 targets, "
         "woolly mammoth vs Asian elephant skin"],
        [""],
        ["Gene CRUSH", "Median CRUSH v1 score (5 kb) over each gene; elephant and mammoth quantile-normalized to "
                       "their mean distribution; dCRUSH = mammoth - elephant (> 0 = more A in mammoth)"],
        ["Targets", "Union of the TF's gene sets in the six ChEA3 libraries (ARCHS4 and GTEx co-expression, "
                    "Enrichr queries, ENCODE, ReMap and literature ChIP-seq), TF itself removed; genes limited to "
                    "those in any ChEA3 library with a CRUSH value (14,484)"],
        ["Adjustment", "Target dCRUSH minus the median dCRUSH of genes listed by a similar number of TFs "
                       "(20 quantile bins), removing the larger elephant A bias of broadly targeted genes"],
        ["Ranking", "TF_score and target_score are signed percentiles (-1 to 1) of TF_dCRUSH_norm and target_r "
                    "over the 1,173 TFs with a CRUSH value; combined_score = their mean; ranked low (elephant) "
                    "to high (mammoth)"],
        [""],
        ["Column", "Description"],
        ["rank", "Rank by combined_score (1 = most elephant A)"],
        ["TF", "Transcription factor"],
        ["direction", "Sign of combined_score"],
        ["elephant_CRUSH / mammoth_CRUSH", "Median CRUSH of the TF gene (not normalized)"],
        ["TF_dCRUSH_norm", "Normalized dCRUSH of the TF gene"],
        ["compartment_switch", "TF gene has opposite CRUSH signs in elephant and mammoth"],
        ["n_targets", "Number of ChEA3 targets tested"],
        ["target_r", "Rank-biserial correlation (2 x AUC - 1) of target vs non-target adjusted dCRUSH "
                     "(> 0 = targets more A in mammoth)"],
        ["target_p", "Two-sided rotation p: dCRUSH ranks circularly shifted along the genome 2,000 times "
                     "(minimum 1/2,001)"],
        ["target_FDR", "Benjamini-Hochberg FDR of target_p over all 1,632 ChEA3 TFs tested"],
        ["TF_score / target_score / combined_score", "See Ranking"],
        ["shown_in", "Figure panel showing this TF"],
        ["concordant_T / discordant_T", "Targets with adjusted dCRUSH >= T in the TF's direction / in the "
                                         "opposite direction (T = 0: split by sign)"],
    ]
    with pd.ExcelWriter(xl, engine="openpyxl") as w:
        pd.DataFrame(readme).to_excel(w, sheet_name="README", index=False, header=False)
        T.to_excel(w, sheet_name="TF ranking", index=False)
    wb = load_workbook(xl)
    for ws in wb:
        for row in ws.iter_rows():
            for c in row:
                c.font = Font(name="Arial", size=10, bold=c.font.bold)
                if isinstance(c.value, float):
                    c.number_format = "0.00E+00" if ws.cell(1, c.column).value in ("target_p", "target_FDR") \
                        else "0.000"
        ws.freeze_panes = "C2" if ws.title != "README" else None
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = (
                max(len(str(c.value or "")) for c in col) + 2 if ws.title != "README" or col[0].column > 1 else 36)
    ws = wb["README"]
    ws["A1"].font = Font(name="Arial", size=11, bold=True)
    ws["A8"].font = ws["B8"].font = Font(name="Arial", size=10, bold=True)
    ws.column_dimensions["B"].width = 110
    wb.save(xl)


if __name__ == "__main__":
    main()
