#!/usr/bin/env python3
import pyBigWig
import pandas as pd
import argparse
import os
import numpy as np
from pathlib import Path


def build_chrom_map(chrom_sizes: dict) -> dict:
    """
    Build a case-insensitive mapping from any chromosome name variant -> canonical name.
    Handles UCSC (chr1), Ensembl (1), mixed case, chrM/chrMT/MT/M aliases, scaffolds.
    """
    mapping = {}
    for canonical in chrom_sizes:
        mapping[canonical.lower()] = canonical
        stripped = canonical.lower()
        if stripped.startswith("chr"):
            mapping.setdefault(stripped[3:], canonical)
        else:
            mapping.setdefault("chr" + stripped, canonical)

    mt_aliases = ["mt", "chrmt", "m", "chrm", "mitochondria", "mitochondrion",
                  "dmel_mitochondrion_genome"]
    for a in mt_aliases:
        if a in mapping:
            canonical = mapping[a]
            for b in mt_aliases:
                mapping.setdefault(b, canonical)

    return mapping


def resolve_chrom(chrom: str, chrom_map: dict):
    return chrom_map.get(chrom.lower())


def quantify_bigwig_over_genes(bw_path, genes_df):
    bw = pyBigWig.open(bw_path)
    chrom_sizes = bw.chroms()
    chrom_map = build_chrom_map(chrom_sizes)

    rows = []
    unresolved = set()

    for _, row in genes_df.iterrows():
        chrom_raw = str(row["Chr"])
        start = int(row["Start"])
        end = int(row["Stop"])
        gname = row["Genename"]

        canonical = resolve_chrom(chrom_raw, chrom_map)
        if canonical is None:
            unresolved.add(chrom_raw)
            rows.append((gname, end - start, 0.0))
            continue

        clen = chrom_sizes[canonical]
        start = max(0, start)
        end = min(clen, end)
        if end <= start:
            rows.append((gname, 0, 0.0))
            continue

        arr = bw.values(canonical, start, end, numpy=True)
        total_signal = 0.0 if arr is None else float(np.nansum(arr))
        rows.append((gname, end - start, total_signal))

    bw.close()

    if unresolved:
        print(f"  WARNING: could not resolve {len(unresolved)} chrom name(s): "
              f"{sorted(unresolved)[:10]}{'...' if len(unresolved) > 10 else ''}")
        print(f"  BigWig chroms (first 10): {sorted(chrom_sizes.keys())[:10]}")

    tmp = pd.DataFrame(rows, columns=["Genename", "length_bp", "raw_sum_coverage"])
    tmp = tmp.groupby("Genename", as_index=False).agg(
        {"length_bp": "max", "raw_sum_coverage": "sum"}
    )
    return tmp


def compute_pseudo_tpm(df):
    df = df.copy()
    df["length_kb"] = df["length_bp"] / 1000.0
    df.loc[df["length_kb"] <= 0, "length_kb"] = np.nan
    df["RPK"] = df["raw_sum_coverage"] / df["length_kb"]

    total_RPK = np.nansum(df["RPK"].values)
    if total_RPK <= 0 or not np.isfinite(total_RPK):
        df["pseudo_TPM"] = 0.0
    else:
        df["pseudo_TPM"] = (df["RPK"] / total_RPK) * 1e6

    df["pseudo_TPM"] = df["pseudo_TPM"].fillna(0.0).astype(float)
    return df[["Genename", "pseudo_TPM"]]


def main():
    parser = argparse.ArgumentParser(
        description="Quantify gene-level pseudo-TPM from bigWig files"
    )
    parser.add_argument("--genes_bed", required=True,
                        help="BED file (with or without header)")
    parser.add_argument("--bigwigs", required=True,
                        help="Comma-separated list of bigWig files")
    parser.add_argument("--labels", required=False,
                        help="Optional sample labels, comma-separated")
    parser.add_argument("--genes_list", required=False,
                        help="Optional file with one Genename per line to subset")
    parser.add_argument("--name_col", required=False, type=int, default=None,
                        help="0-based column index for gene name in a headerless BED "
                             "(default: 3). Use this if auto-detection picks the wrong column.")
    parser.add_argument("--out_tsv", required=True,
                        help="Output long-format TSV file")
    args = parser.parse_args()

    with open(args.genes_bed) as fh:
        first_line = fh.readline().rstrip("\n").split("\t")
        second_line = fh.readline().rstrip("\n").split("\t")

    def _looks_numeric(s):
        try:
            int(s); return True
        except ValueError:
            return False

    has_header = not _looks_numeric(first_line[1]) if len(first_line) > 1 else True

    print("\n=== BED diagnostic ===")
    print(f"  Header detected : {has_header}")
    print(f"  Row 0 (first 8) : {first_line[:8]}")
    if second_line:
        print(f"  Row 1 (first 8) : {second_line[:8]}")

    if has_header:
        genes_df = pd.read_csv(args.genes_bed, sep="\t", header=0)
        genes_df = genes_df.rename(columns={
            "chr": "Chr", "chrom": "Chr",
            "start": "Start", "end": "Stop",
            "gene_name": "Genename", "name": "Genename",
        })
    else:
        genes_df = pd.read_csv(args.genes_bed, sep="\t", header=None)

        name_idx = args.name_col if args.name_col is not None else 3

        col_map = {0: "Chr", 1: "Start", 2: "Stop", name_idx: "Genename"}
        genes_df = genes_df.rename(columns=col_map)

        sample_names = genes_df["Genename"].unique()[:8].tolist() if "Genename" in genes_df.columns else "N/A"
        print(f"\n  Headerless BED -- Genename from column index {name_idx}.")
        print(f"  Sample Genename values : {sample_names}")
        print(f"  If these look wrong (e.g. '+'/'-'), rerun with --name_col <correct_index>.\n")

    required = {"Chr", "Start", "Stop", "Genename"}
    missing = required - set(genes_df.columns)
    if missing:
        raise ValueError(
            f"BED file is missing required columns: {missing}.\n"
            f"Found columns: {list(genes_df.columns)}\n"
            f"If your BED has no header, try --name_col <0-based index of gene name column>."
        )

    print(f"  BED loaded: {len(genes_df)} rows, {genes_df['Genename'].nunique()} unique gene names.")

    if args.genes_list:
        keep = {line.strip() for line in open(args.genes_list)}
        genes_df = genes_df[genes_df["Genename"].isin(keep)].copy()
        print(f"  After genes_list filter: {len(genes_df)} rows.")

    bw_paths = [p.strip() for p in args.bigwigs.split(",") if p.strip()]
    labels = (
        [os.path.splitext(os.path.basename(p))[0] for p in bw_paths]
        if not args.labels
        else [l.strip() for l in args.labels.split(",")]
    )

    if len(bw_paths) != len(labels):
        raise ValueError(f"Mismatch: {len(bw_paths)} bigWigs vs {len(labels)} labels")

    all_rows = []
    for bw_path, label in zip(bw_paths, labels):
        print(f"\nQuantifying [{label}]  <-  {bw_path}")
        raw = quantify_bigwig_over_genes(bw_path, genes_df)
        print(raw.head())
        tpm = compute_pseudo_tpm(raw)
        tpm["Sample"] = label
        all_rows.append(tpm)

    final_df = pd.concat(all_rows, axis=0)
    final_df = final_df[["Genename", "Sample", "pseudo_TPM"]]
    final_df.to_csv(args.out_tsv, sep="\t", index=False)

    print(f"\nWrote : {args.out_tsv}")
    print(f"Rows  : {len(final_df)}  |  Unique genes : {final_df['Genename'].nunique()}")


if __name__ == "__main__":
    main()
