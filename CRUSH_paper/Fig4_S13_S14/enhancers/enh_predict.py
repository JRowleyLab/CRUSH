#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve

HERE = Path(__file__).resolve().parent
RES = 1000
ENH = Path("enhancers")
TRK = Path("compBedgraphs")
GENES = "hg19genes.bed"
SIZES = "hg19.sizes"
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
LINES = {"A549": ("a549", "A549", "a549"), "HCT116": ("hct", "HCT", "hct"),
         "HeLa-S3": ("helaS3", "HELA", "hela"), "HepG2": ("hepg2", "HEPG2", "hepg2"),
         "HUVEC": ("huvec", "HUVEC", "huvec"), "IMR90": ("imr90", "IMR90", "imr90"),
         "K562": ("k562", "K562", "k562")}
COL = dict(zip(LINES, plt.get_cmap("tab10").colors))


def bed(path, cols=(0, 1, 2)):
    gz = open(path, "rb").read(2) == b"\x1f\x8b"
    t = pd.read_csv(path, sep="\t", header=None, usecols=list(cols), comment="#", dtype={cols[0]: str},
                    compression="gzip" if gz else None)
    t.columns = ["c", "s", "e"]
    t["c"] = t["c"].str.replace("chr", "", regex=False)
    return t[t.c.isin(CHROMS)]


def mark(sizes, t):
    """Boolean per-bin arrays: bins touched by any interval in t."""
    out = {c: np.zeros(sizes[c] // RES + 1, bool) for c in CHROMS}
    for c, g in t.groupby("c"):
        a = out[c]
        for s, e in zip(g.s.to_numpy() // RES, (g.e.to_numpy() - 1) // RES):
            a[s:min(e, len(a) - 1) + 1] = True
    return out


def flat(d):
    return np.concatenate([d[c] for c in CHROMS])


def main():
    sizes = {l.split()[0]: int(l.split()[1]) for l in open(SIZES) if l.split()[0] in CHROMS}
    g = pd.read_csv(GENES, sep="\t", header=None, usecols=[0, 1, 2, 5], dtype={0: str})
    g.columns = ["c", "s", "e", "strand"]
    tss = np.where(g.strand == "+", g.s, g.e - 1)
    prom = flat(mark(sizes, pd.DataFrame({"c": g.c, "s": np.maximum(tss - 1000, 0), "e": tss + 1001})
                     .query("c in @CHROMS")))
    body = flat(mark(sizes, g[g.c.isin(CHROMS)][["c", "s", "e"]]))

    rows, comp, curves = [], [], {}
    for lab, (trk, edb, pk) in LINES.items():
        t = pd.read_csv(TRK / f"{trk}_noY_1kb.bedgraph", sep="\t", header=None, comment="t",
                        names=["c", "s", "e", "v"], dtype={"c": str})
        t["c"] = t["c"].str.replace("chr", "", regex=False)
        cr = {c: np.full(sizes[c] // RES + 1, np.nan) for c in CHROMS}
        for c, gg in t[t.c.isin(CHROMS)].groupby("c"):
            cr[c][(gg.s.to_numpy() // RES).clip(0, len(cr[c]) - 1)] = gg.v.to_numpy()
        score = flat(cr)
        ok = np.isfinite(score) & (score != 0)
        enh = flat(mark(sizes, bed(ENH / f"enhdb_{edb}.bed")))
        k27 = flat(mark(sizes, bed(next((ENH / "h3k27ac").glob(f"{pk}_*.bed")))))
        k4 = flat(mark(sizes, bed(next((ENH / "h3k4me1").glob(f"{pk}_*.bed")))))
        act = enh & k27 & k4
        A = score > 0

        Aok = A & ok
        nA = Aok.sum()
        peak = (k27 | k4) & ~enh & ~prom
        comp.append({"line": lab, "n_A_bins": int(nA), "frac_genome_A": nA / ok.sum(),
                     "promoter": (Aok & prom).sum() / nA,
                     "enhancer (not promoter)": (Aok & enh & ~prom).sum() / nA,
                     "other H3K27ac/H3K4me1 peak": (Aok & peak).sum() / nA,
                     "gene body": (Aok & body & ~prom & ~enh & ~peak).sum() / nA,
                     "intergenic": (Aok & ~body & ~prom & ~enh & ~peak).sum() / nA,
                     "active enhancer, of A": (Aok & act).sum() / nA,
                     "B: promoter": (~A & ok & prom).sum() / (~A & ok).sum(),
                     "B: enhancer (not promoter)": (~A & ok & enh & ~prom).sum() / (~A & ok).sum(),
                     "B: other H3K27ac/H3K4me1 peak": (~A & ok & peak).sum() / (~A & ok).sum(),
                     "B: gene body": (~A & ok & body & ~prom & ~enh & ~peak).sum() / (~A & ok).sum()})
        for name, y, sub in (("enhancer", enh, ok), ("active enhancer", act, ok),
                             ("enhancer, distal", enh, ok & ~prom), ("active enhancer, distal", act, ok & ~prom)):
            yy, ss = y[sub], score[sub]
            prev = yy.mean()
            r = {"line": lab, "label": name, "n_bins": int(sub.sum()), "n_pos": int(yy.sum()),
                 "prevalence": prev,
                 "precision_A": yy[ss > 0].mean(), "recall_A": (ss[yy] > 0).mean(),
                 "enrichment_A": yy[ss > 0].mean() / prev,
                 "AUROC": roc_auc_score(yy, ss), "AUPRC": average_precision_score(yy, ss)}
            order = np.argsort(-ss)
            for q in (1, 5, 10):
                top = order[: int(len(ss) * q / 100)]
                r[f"precision_top{q}pct"] = yy[top].mean()
            rows.append(r)
            if name in ("enhancer", "active enhancer"):
                curves[(lab, name)] = roc_curve(yy, ss)[:2]
            print(lab, name, {k: round(v, 3) if isinstance(v, float) else v for k, v in r.items() if k not in ("line", "label")}, flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(HERE / "enh_predict.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(comp).to_csv(HERE / "a_composition.tsv", sep="\t", index=False, float_format="%.4g")

    cd = pd.DataFrame(comp).set_index("line")
    for notext in (False, True):
        fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), gridspec_kw={"width_ratios": [1, 1, 1.1]})
        for k, name in enumerate(("enhancer", "active enhancer")):
            for lab in LINES:
                fpr, tpr = curves[(lab, name)]
                step = max(1, len(fpr) // 2000)
                ax[k].plot(fpr[::step], tpr[::step], color=COL[lab], lw=1.1, label=lab)
            ax[k].plot([0, 1], [0, 1], ls=":", color="grey", lw=0.8)
            ax[k].set_aspect("equal")
        bottom = np.zeros(len(cd))
        for part, c in (("promoter", "#7570b3"), ("enhancer (not promoter)", "#d62728"),
                        ("other H3K27ac/H3K4me1 peak", "#fc9272"), ("gene body", "#969696"), ("intergenic", "#e0e0e0")):
            ax[2].bar(range(len(cd)), cd[part] * 100, bottom=bottom, color=c, edgecolor="black", lw=0.4, label=part)
            bottom += cd[part].to_numpy() * 100
        ax[2].set_xticks(range(len(cd)))
        for a in ax:
            for s in ("top", "right"):
                a.spines[s].set_visible(False)
            a.tick_params(labelsize=8)
        if notext:
            ax[2].set_xticklabels([])
        else:
            for k, name in enumerate(("enhancerDB enhancers", "active enhancers (+H3K27ac, +H3K4me1)")):
                ax[k].set_title(name, fontsize=9)
                ax[k].set_xlabel("False positive rate", fontsize=9)
                ax[k].set_ylabel("True positive rate", fontsize=9)
            ax[0].legend(frameon=False, fontsize=7, loc="lower right")
            ax[2].set_xticklabels(cd.index, rotation=45, ha="right", fontsize=8)
            ax[2].set_ylabel("% of A compartment bins", fontsize=9)
            ax[2].legend(frameon=False, fontsize=7, bbox_to_anchor=(1.0, 1.0), loc="upper left")
        fig.tight_layout()
        stem = "enh_predict" + ("_notext" if notext else "")
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"{stem}.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
