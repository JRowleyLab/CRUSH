#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from itertools import product

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import gaussian_kde, spearmanr


def build_chrom_map(chrom_sizes: dict) -> dict:
    mapping = {}
    for canonical in chrom_sizes:
        mapping[canonical.lower()] = canonical
        stripped = canonical.lower()
        if stripped.startswith("chr"):
            mapping.setdefault(stripped[3:], canonical)
        else:
            mapping.setdefault("chr" + stripped, canonical)
    for a in ["mt", "chrmt", "m", "chrm", "mitochondria", "mitochondrion"]:
        if a in mapping:
            canon = mapping[a]
            for b in ["mt", "chrmt", "m", "chrm", "mitochondria", "mitochondrion"]:
                mapping.setdefault(b, canon)
    return mapping


def resolve_chrom(chrom: str, chrom_map: dict):
    return chrom_map.get(chrom.lower())


def load_genes_bed(path: str, name_col: int = None) -> pd.DataFrame:
    with open(path, encoding="latin-1") as fh:
        first  = fh.readline().rstrip("\n").split("\t")
        second = fh.readline().rstrip("\n").split("\t")

    def _numeric(s):
        try: int(s); return True
        except ValueError: return False

    has_header = not _numeric(first[1]) if len(first) > 1 else True

    print("\n=== BED diagnostic ===")
    print(f"  Header          : {has_header}")
    print(f"  Row 0 (8 cols)  : {first[:8]}")
    print(f"  Row 1 (8 cols)  : {second[:8]}")

    if has_header:
        df = pd.read_csv(path, sep="\t", header=0,
                         dtype={"Chr": str, "Strand": str}, encoding="latin-1")
        df = df.rename(columns={
            "chr": "Chr", "chrom": "Chr", "start": "Start", "end": "Stop",
            "strand": "Strand", "gene_name": "Genename", "name": "Genename",
        })
    else:
        df = pd.read_csv(path, sep="\t", header=None,
                         dtype={0: str}, encoding="latin-1")
        idx = name_col if name_col is not None else 3
        col_map = {0: "Chr", 1: "Start", 2: "Stop", idx: "Genename"}
        if 5 not in col_map and df.shape[1] > 5:
            col_map[5] = "Strand"
        df = df.rename(columns=col_map)
        print(f"  Headerless BED  : Genename col {idx}, "
              f"samples: {df['Genename'].unique()[:5].tolist()}")

    if "Strand" not in df.columns:
        print("  WARNING: no Strand column, defaulting all to '+'")
        df["Strand"] = "+"

    for col in ("Chr", "Genename", "Strand"):
        df[col] = df[col].astype(str).str.strip()

    required = {"Chr", "Start", "Stop", "Genename"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"BED missing columns: {missing}\n"
                         f"Found: {list(df.columns)}\n"
                         f"Tip: use --name_col to set 0-based gene name column.")

    print(f"  Loaded          : {len(df)} rows, {df['Genename'].nunique()} unique names\n")
    return df


def extract_compartment_scores(comp_bw: str, genes_df: pd.DataFrame,
                               window_bp: int) -> pd.DataFrame:
    bw = pyBigWig.open(comp_bw)
    chrom_sizes = bw.chroms()
    chrom_map   = build_chrom_map(chrom_sizes)
    unresolved  = set()
    rows = []

    for _, row in genes_df.iterrows():
        chrn    = str(row["Chr"])
        strand  = str(row.get("Strand", "+"))
        g_start = int(row["Start"])
        g_stop  = int(row["Stop"])
        gname   = row["Genename"]

        canonical = resolve_chrom(chrn, chrom_map)
        if canonical is None:
            unresolved.add(chrn); rows.append((gname, np.nan)); continue

        tss = g_start if strand == "+" else g_stop
        raw_s, raw_e = (tss, tss + window_bp) if strand == "+" \
                       else (tss - window_bp, tss)

        clen = chrom_sizes[canonical]
        s = max(0, min(raw_s, raw_e))
        e = min(clen, max(raw_s, raw_e))
        if e <= s:
            rows.append((gname, np.nan)); continue

        vals     = bw.stats(canonical, s, e, type="mean")
        mean_val = np.nan if not vals or vals[0] is None else vals[0]
        rows.append((gname, mean_val))

    bw.close()
    if unresolved:
        print(f"  WARNING: unresolved chroms ({len(unresolved)}): "
              f"{sorted(unresolved)[:8]}")

    df = pd.DataFrame(rows, columns=["Genename", "CompartmentScore"])
    df = df.groupby("Genename", as_index=False)["CompartmentScore"].mean()
    n_valid = df["CompartmentScore"].notna().sum()
    print(f"  Comp scores: {n_valid} genes with values, {len(df)-n_valid} NaN")
    return df


def quantify_signal_bigwig(bw_path: str, genes_df: pd.DataFrame) -> pd.DataFrame:
    """
    Works with any signal bigWig: RNA-seq, GRO-seq, BRU-seq,
    Histone ChIP, ATAC-seq, etc.
    Computes gene-body sum coverage then normalises to pseudo-TPM.
    """
    bw = pyBigWig.open(bw_path)
    chrom_sizes = bw.chroms()
    chrom_map   = build_chrom_map(chrom_sizes)
    unresolved  = set()
    rows = []

    for _, row in genes_df.iterrows():
        chrom_raw = str(row["Chr"])
        start     = int(row["Start"])
        end       = int(row["Stop"])
        gname     = row["Genename"]

        canonical = resolve_chrom(chrom_raw, chrom_map)
        if canonical is None:
            unresolved.add(chrom_raw); rows.append((gname, end-start, 0.0)); continue

        clen  = chrom_sizes[canonical]
        start = max(0, start)
        end   = min(clen, end)
        if end <= start:
            rows.append((gname, 0, 0.0)); continue

        arr          = bw.values(canonical, start, end, numpy=True)
        total_signal = 0.0 if arr is None else float(np.nansum(arr))
        rows.append((gname, end - start, total_signal))

    bw.close()
    if unresolved:
        print(f"  WARNING: unresolved chroms ({len(unresolved)}): "
              f"{sorted(unresolved)[:8]}")

    tmp = pd.DataFrame(rows, columns=["Genename", "length_bp", "raw_sum"])
    tmp = tmp.groupby("Genename", as_index=False).agg(
        {"length_bp": "max", "raw_sum": "sum"})

    tmp["length_kb"] = tmp["length_bp"] / 1000.0
    tmp.loc[tmp["length_kb"] <= 0, "length_kb"] = np.nan
    tmp["RPK"]   = tmp["raw_sum"] / tmp["length_kb"]
    total        = np.nansum(tmp["RPK"].values)
    tmp["pseudo_TPM"] = 0.0 if (total <= 0 or not np.isfinite(total)) \
                        else (tmp["RPK"] / total) * 1e6
    tmp["pseudo_TPM"] = tmp["pseudo_TPM"].fillna(0.0)
    return tmp[["Genename", "pseudo_TPM"]]


def load_signal_tsv(tsv_path: str, label: str) -> pd.DataFrame:
    """
    Accepts two formats:
      A) Genename <TAB> pseudo_TPM                  (single-sample)
      B) Genename <TAB> Sample <TAB> pseudo_TPM     (multi-sample, filter by label)
    Returns DataFrame with columns: Genename, pseudo_TPM
    """
    df = pd.read_csv(tsv_path, sep="\t", encoding="latin-1")
    df["Genename"] = df["Genename"].astype(str).str.strip()

    if "Sample" in df.columns:
        df["Sample"] = df["Sample"].astype(str).str.strip()
        available = df["Sample"].unique().tolist()
        filtered = df[df["Sample"] == label][["Genename", "pseudo_TPM"]].copy()
        if filtered.empty:
            if len(available) == 1:
                print(f"  NOTE: label '{label}' not in {tsv_path}; "
                      f"using only available sample '{available[0]}'")
                filtered = df[df["Sample"] == available[0]][["Genename", "pseudo_TPM"]].copy()
            else:
                raise ValueError(f"Label '{label}' not found in {tsv_path}.\n"
                                 f"Available: {available}")
        return filtered
    else:
        return df[["Genename", "pseudo_TPM"]].copy()


def rank_transform(s: pd.Series) -> pd.Series:
    if len(s) <= 1:
        return pd.Series([0.0] * len(s), index=s.index)
    return (s.rank(method="average") - 1) / (len(s) - 1)


KDE_QUANTILES = [0.05, 0.10, 0.20, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.99]


def draw_contour_panel(ax, x: np.ndarray, y: np.ndarray,
                       title: str, cmap, signal_type: str = "Signal",
                       min_genes: int = 50) -> tuple:
    """Draw one KDE contour panel onto ax. Returns (rho, pval, n)."""
    n = len(x)
    if n < min_genes:
        ax.text(0.5, 0.5, f"n={n} < {min_genes}", ha="center", va="center",
                transform=ax.transAxes, fontsize=8)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_title(title, fontsize=8)
        return np.nan, np.nan, n

    try:
        kde  = gaussian_kde(np.vstack([x, y]))
        gx   = np.linspace(0, 1, 120)
        gy   = np.linspace(0, 1, 120)
        X, Y = np.meshgrid(gx, gy)
        Z    = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
        Z    = (Z - Z.min()) / (Z.max() - Z.min() + 1e-12)
        levs = np.unique(np.quantile(Z, KDE_QUANTILES))
        ax.contourf(X, Y, Z, levels=levs, extend="max", cmap=cmap)
        ax.contour(X, Y, Z,  levels=levs, linewidths=0.5, colors="k", alpha=0.3)
    except Exception as exc:
        ax.scatter(x, y, s=3, alpha=0.3, color="steelblue")
        print(f"  KDE fallback for '{title}': {exc}")

    rho, pval = spearmanr(x, y)
    ax.set_title(f"{title}\nρ={rho:.3f} n={n:,}", fontsize=8)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.tick_params(labelsize=7)
    return rho, pval, n


def build_merged(comp_df: pd.DataFrame, signal_df: pd.DataFrame) -> pd.DataFrame:
    merged = comp_df.merge(signal_df, on="Genename", how="inner") \
                    .dropna(subset=["CompartmentScore", "pseudo_TPM"])
    if merged.empty:
        return merged
    merged["CompPct"]   = rank_transform(merged["CompartmentScore"])
    merged["SignalPct"] = rank_transform(merged["pseudo_TPM"])
    return merged


def main():
    parser = argparse.ArgumentParser(
        description="Compartment × Any-Signal KDE contour pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument("--mode", choices=["single", "matrix"], default="single",
                        help="'single': one-to-one comp×signal pairs. "
                             "'matrix': every comp × every signal (N×M grid).")

    parser.add_argument("--genes_bed",    required=True,
                        help="Gene BED file (headerless or with header).")
    parser.add_argument("--comp_bws",     required=True,
                        help="Comma-separated compartment bigWig paths.")
    parser.add_argument("--comp_labels",  required=True,
                        help="Comma-separated labels for compartment bigWigs.")

    sig_group = parser.add_mutually_exclusive_group(required=True)
    sig_group.add_argument("--signal_tsvs",
                        help="Comma-separated pre-computed signal TSV paths "
                             "(Genename + pseudo_TPM columns). Use when signal is "
                             "already quantified (pseudo-TPM, normalised counts, etc.).")
    sig_group.add_argument("--signal_bws",
                        help="Comma-separated signal bigWig paths to quantify on the fly. "
                             "Works with ANY bigWig type: RNA-seq, GRO-seq, BRU-seq, "
                             "Histone ChIP-seq, ATAC-seq, etc.")

    parser.add_argument("--signal_labels", required=True,
                        help="Comma-separated labels for signal files.")
    parser.add_argument("--signal_type",  default="Signal",
                        help="Y-axis label describing the signal type, e.g. "
                             "'BRU-seq', 'GRO-seq', 'H3K27ac ChIP', 'ATAC-seq'. "
                             "Default: 'Signal'.")

    parser.add_argument("--out_dir", required=True)

    parser.add_argument("--window_bp", type=int, default=10_000,
                        help="TSS window in bp, stranded (default: 10000).")
    parser.add_argument("--min_genes", type=int, default=50,
                        help="Min genes per panel to attempt KDE (default: 50).")
    parser.add_argument("--name_col",  type=int, default=None,
                        help="0-based gene name col in headerless BED (default: 3).")
    parser.add_argument("--dpi",       type=int, default=300)
    parser.add_argument("--cmap",      default="magma")

    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmap = plt.get_cmap(args.cmap)

    comp_bws    = [p.strip() for p in args.comp_bws.split(",")    if p.strip()]
    comp_labels = [l.strip() for l in args.comp_labels.split(",") if l.strip()]
    sig_labels  = [l.strip() for l in args.signal_labels.split(",") if l.strip()]

    if args.signal_bws:
        sig_paths = [p.strip() for p in args.signal_bws.split(",") if p.strip()]
        use_bw    = True
    else:
        sig_paths = [p.strip() for p in args.signal_tsvs.split(",") if p.strip()]
        if len(sig_paths) == 1 and len(sig_labels) > 1:
            print(f"  Single signal TSV — will filter by label for each sample.")
            sig_paths = sig_paths * len(sig_labels)
        use_bw = False

    if len(comp_bws) != len(comp_labels):
        sys.exit(f"ERROR: --comp_bws ({len(comp_bws)}) and "
                 f"--comp_labels ({len(comp_labels)}) must match.")
    if len(sig_paths) != len(sig_labels):
        sys.exit(f"ERROR: signal paths ({len(sig_paths)}) and "
                 f"--signal_labels ({len(sig_labels)}) must match.")
    if args.mode == "single" and len(comp_bws) != len(sig_paths):
        sys.exit(f"ERROR: in 'single' mode, --comp_bws ({len(comp_bws)}) and "
                 f"signal files ({len(sig_paths)}) must have the same count.\n"
                 f"Use --mode matrix for all-vs-all.")

    genes = load_genes_bed(args.genes_bed, name_col=args.name_col)

    print("\n=== Extracting compartment scores ===")
    comp_dfs = {}
    for label, bw in zip(comp_labels, comp_bws):
        print(f"  [{label}]  {bw}")
        df = extract_compartment_scores(bw, genes, args.window_bp)
        df.to_csv(out_dir / f"{label}_compartment_scores.tsv", sep="\t", index=False)
        comp_dfs[label] = df

    print(f"\n=== Loading signal ({args.signal_type}) ===")
    sig_dfs = {}
    for label, path in zip(sig_labels, sig_paths):
        print(f"  [{label}]  {path}")
        if use_bw:
            df = quantify_signal_bigwig(path, genes)
            df["Sample"] = label
            df[["Genename","Sample","pseudo_TPM"]].to_csv(
                out_dir / f"{label}_signal_pseudoTPM.tsv", sep="\t", index=False)
        else:
            df = load_signal_tsv(path, label)
        print(f"    {len(df)} genes loaded")
        sig_dfs[label] = df

    if args.mode == "single":
        print("\n=== MODE: single ===")
        all_stats = []
        for c_label, s_label in zip(comp_labels, sig_labels):
            merged = build_merged(comp_dfs[c_label], sig_dfs[s_label])
            tag    = f"{c_label}_vs_{s_label}"
            merged.to_csv(out_dir / f"{tag}_merged.tsv", sep="\t", index=False)

            fig, ax = plt.subplots(figsize=(4.5, 4.2), facecolor="white")
            x = merged["CompPct"].to_numpy()
            y = merged["SignalPct"].to_numpy()
            rho, pval, n = draw_contour_panel(ax, x, y, tag, cmap,
                                              args.signal_type, args.min_genes)

            ax.set_xlabel("Rank Compartment Score", fontsize=9)
            ax.set_ylabel(f"Rank {args.signal_type}", fontsize=9)

            norm = mpl.colors.Normalize(vmin=0, vmax=1)
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                         ax=ax, fraction=0.046, pad=0.04).set_label(
                "Relative density (KDE)", fontsize=7)
            fig.tight_layout()
            out_png = out_dir / f"{tag}_contour.png"
            fig.savefig(out_png, dpi=args.dpi, facecolor="white", bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved → {out_png}  ρ={rho:.4f}  n={n}")
            all_stats.append({"Comp": c_label, "Signal": s_label,
                               "N_genes": n, "Spearman_rho": rho,
                               "Spearman_pval": f"{pval:.3e}" if not np.isnan(pval) else "NA"})

        pd.DataFrame(all_stats).to_csv(out_dir / "stats.tsv", sep="\t", index=False)

    else:
        print(f"\n=== MODE: matrix  "
              f"({len(comp_labels)} comps × {len(sig_labels)} signals) ===")

        nrows = len(comp_labels)
        ncols = len(sig_labels)

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(3.2 * ncols, 3.0 * nrows),
            squeeze=False,
            facecolor="white"
        )

        all_stats = []
        for i, c_label in enumerate(comp_labels):
            for j, s_label in enumerate(sig_labels):
                ax     = axes[i][j]
                merged = build_merged(comp_dfs[c_label], sig_dfs[s_label])
                tag    = f"{c_label}_vs_{s_label}"
                merged.to_csv(out_dir / f"{tag}_merged.tsv", sep="\t", index=False)

                x = merged["CompPct"].to_numpy()   if not merged.empty else np.array([])
                y = merged["SignalPct"].to_numpy()  if not merged.empty else np.array([])

                rho, pval, n = draw_contour_panel(ax, x, y, tag, cmap,
                                                  args.signal_type, args.min_genes)

                if j == 0:
                    ax.set_ylabel(f"Rank {args.signal_type}\n({c_label} comp)",
                                  fontsize=7)
                if i == nrows - 1:
                    ax.set_xlabel(f"Rank Comp\n({s_label} signal)", fontsize=7)

                fig_single, ax_single = plt.subplots(figsize=(4.5, 4.2), facecolor="white")
                draw_contour_panel(ax_single, x, y, tag, cmap,
                                   args.signal_type, args.min_genes)
                ax_single.set_xlabel("Rank Compartment Score", fontsize=9)
                ax_single.set_ylabel(f"Rank {args.signal_type}", fontsize=9)
                norm = mpl.colors.Normalize(vmin=0, vmax=1)
                fig_single.colorbar(
                    mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                    ax=ax_single, fraction=0.046, pad=0.04
                ).set_label("Relative density (KDE)", fontsize=7)
                fig_single.tight_layout()
                fig_single.savefig(out_dir / f"{tag}_contour.png",
                                   dpi=args.dpi, facecolor="white", bbox_inches="tight")
                plt.close(fig_single)

                all_stats.append({
                    "Comp": c_label, "Signal": s_label, "N_genes": n,
                    "Spearman_rho": round(rho, 4) if not np.isnan(rho) else np.nan,
                    "Spearman_pval": f"{pval:.3e}" if not np.isnan(pval) else "NA"
                })
                print(f"  {c_label} × {s_label}  ρ={rho:.4f}  n={n}")

        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        fig.subplots_adjust(right=0.88, wspace=0.35, hspace=0.50)
        cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.70])
        mpl.colorbar.ColorbarBase(
            cbar_ax, cmap=cmap, norm=norm, orientation="vertical"
        ).set_label("Relative density (KDE, 0–1)", fontsize=8)

        for i, c_label in enumerate(comp_labels):
            axes[i][0].annotate(
                f"{c_label} comp", xy=(0, 0.5),
                xytext=(-axes[i][0].yaxis.labelpad - 30, 0),
                xycoords=axes[i][0].yaxis.label,
                textcoords="offset points",
                ha="right", va="center", fontsize=8, fontweight="bold",
                rotation=90
            )

        grid_path = out_dir / "matrix_contours.png"
        fig.savefig(grid_path, dpi=args.dpi, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        print(f"\n  Grid figure → {grid_path}")

        stats_df = pd.DataFrame(all_stats)
        stats_df.to_csv(out_dir / "matrix_stats.tsv", sep="\t", index=False)

        pivot = stats_df.pivot(index="Comp", columns="Signal", values="Spearman_rho")
        pivot = pivot.reindex(index=comp_labels, columns=sig_labels)
        print("\n  Spearman ρ matrix:")
        print(pivot.to_string())
        pivot.to_csv(out_dir / "matrix_rho_pivot.tsv", sep="\t")

    print(f"\nDone. All outputs in: {out_dir}")


if __name__ == "__main__":
    main()
