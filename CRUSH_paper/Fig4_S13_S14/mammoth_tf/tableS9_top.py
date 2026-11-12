#!/usr/bin/env python3
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font

HERE = Path(__file__).resolve().parent
NTOP = 60
COLS = {"TF": "TF", "elephant_CRUSH": "elephant_CRUSH", "mammoth_CRUSH": "mammoth_CRUSH",
        "TF_dCRUSH_norm": "TF_dCRUSH_norm", "compartment_switch": "compartment_switch",
        "n_targets": "n_targets", "target_r": "target_r", "target_FDR": "target_FDR",
        "combined_score": "combined_score", "shown_in": "shown_in"}


def main():
    r = pd.read_csv(HERE / "tableS_TFrank.tsv", sep="\t", keep_default_na=False)
    e = r.head(NTOP).assign(direction="A in elephant")
    m = r.tail(NTOP).iloc[::-1].assign(direction="A in mammoth")
    T = pd.concat([e, m])
    T.insert(0, "rank_in_direction", list(range(1, NTOP + 1)) * 2)
    T = T[["direction", "rank_in_direction", *COLS]]
    assert set(r.TF[r.shown_in != ""]) <= set(T.TF)
    T.to_csv(HERE / "tableS9_top.tsv", sep="\t", index=False, float_format="%.4g")
    xl = HERE / "tableS9_top.xlsx"
    T.to_excel(xl, sheet_name="Table S9", index=False)
    wb = load_workbook(xl)
    ws = wb.active
    for row in ws.iter_rows():
        for c in row:
            c.font = Font(name="Arial", size=10, bold=c.row == 1)
            if isinstance(c.value, float):
                c.number_format = "0.00E+00" if ws.cell(1, c.column).value == "target_FDR" else "0.00"
    for col in ws.columns:
        ws.column_dimensions[col[0].column_letter].width = max(len(str(c.value or "")) for c in col) + 2
    ws.freeze_panes = "C2"
    wb.save(xl)
    print(T.groupby("direction").agg(n=("TF", "size"), FDR05=("target_FDR", lambda x: (x < 0.05).sum()),
                                      switch=("compartment_switch", lambda x: (x == "yes").sum())))
    print(T[T.shown_in != ""][["direction", "rank_in_direction", "TF"]].to_string(index=False))


if __name__ == "__main__":
    main()
