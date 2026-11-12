#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
from concurrent.futures import ThreadPoolExecutor

import multistate_crush as ms

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
FIGS = HERE.parent / "figures"
R = 100_000
MARKS = Path("ultraMap_Histones_hg19")
TRACKS = {"K36": "H3K36me3.bigwig", "K27ac": "H3K27ac.bigwig", "GRO": "GMgro.bigwig",
          "K27me3": "H3K27me3.bigwig", "K9me3": "H3K9me3.bigwig", "HP1": "CBX5_HP1a.bigWig"}
GENES = "hg19genes.bed"
BBINS = "Bbins.bed"
ROUNDS = 5


def bin_signal(path: Path, chrom: str, n: int) -> np.ndarray:
    b = pyBigWig.open(str(path))
    name = next((k for k in (chrom, "chr" + chrom) if k in b.chroms()), None)
    out = np.full(n, np.nan)
    if name:
        m = min(n, b.chroms()[name] // R)
        out[:m] = [np.nan if x is None else x for x in b.stats(name, 0, m * R, nBins=m)]
    b.close()
    return out


def density(bed: str, chrom: str, n: int) -> np.ndarray:
    t = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"], dtype={"c": str})
    t = t[t["c"].str.replace("chr", "", regex=False) == chrom]
    d = np.zeros(n)
    np.add.at(d, np.clip(((t["s"] + t["e"]) // 2 // R).to_numpy(), 0, n - 1), 1)
    return d


def main():
    sizes = {l.split()[0]: int(l.split()[1]) for l in open(ms.SIZES)}
    lab = ms.labels_at(R, sizes)
    rows = []
    for c in ms.CHROMS:
        n = sizes[c] // R + 1
        df = pd.DataFrame({k: bin_signal(MARKS / f, c, n) for k, f in TRACKS.items()})
        df["genes"] = density(GENES, c, n)
        df["lowgc"] = density(BBINS, c, n)
        df["c"], df["i"], df["rao"] = c, np.arange(n), lab[c]
        rows.append(df)
    d = pd.concat(rows, ignore_index=True)
    for k in TRACKS:
        d[k] = (d[k] - d[k].mean()) / d[k].std()
    d["act"] = d[["K36", "K27ac", "GRO"]].mean(axis=1)
    q = lambda col, p: d[col].quantile(p)

    rules = {
        "A1": (d.act >= q("act", .85)) & (d.genes >= q("genes", .75)),
        "A2": (d.act >= q("act", .50)) & (d.act < q("act", .85)) & (d.K9me3 >= q("K9me3", .50)),
        "B1": (d.K27me3 >= q("K27me3", .85)) & (d.act < q("act", .50)),
        "B2": (np.maximum(d.K9me3, d.HP1) >= np.maximum(q("K9me3", .85), q("HP1", .85)))
              & (d.act < q("act", .50)) & (d.K27me3 < q("K27me3", .50)),
        "B3": (d.act < q("act", .30)) & (d.genes <= q("genes", .25)) & (d.lowgc >= q("lowgc", .85)),
    }
    hits = sum(m.astype(int) for m in rules.values())
    seed_state = pd.Series("", index=d.index, dtype=object)
    for st, m in rules.items():
        seed_state[m & (hits == 1)] = st
    d["seed0"] = seed_state

    res = {"seed_counts": d["seed0"].value_counts().drop("", errors="ignore").to_dict(),
           "seed_purity_vs_rao": {st: float((d.loc[d.seed0 == st, "rao"] == st).mean()) for st in ms.STATES},
           "mark_means_by_rao_state": d[d.rao != ""].groupby("rao")[["act", "K27me3", "K9me3", "HP1", "genes", "lowgc"]]
           .mean().round(3).to_dict()}

    with ThreadPoolExecutor(8) as ex:
        mats = dict(zip(ms.CHROMS, ex.map(lambda c: ms.oe_matrix(c, R, sizes[c] // R + 1), ms.CHROMS)))

    test = (d.rao != "") & (d.seed0 == "")
    seeds = d["seed0"].copy()
    per_round = []
    for rnd in range(ROUNDS):
        sc = []
        for c in ms.CHROMS:
            sub = d.index[d.c == c]
            s = {st: d.loc[sub, "i"].to_numpy()[(seeds[sub] == st).to_numpy()] for st in ms.STATES}
            sc.append(ms.state_scores(mats[c], s))
        z = ms.standardise(np.vstack(sc))
        z = np.where(np.isfinite(z), z, -np.inf)
        ok = np.isfinite(z).sum(1) >= 2
        pred = np.full(len(d), "", dtype=object)
        pred[ok] = np.array(ms.STATES)[np.argmax(z[ok], axis=1)]
        srt = np.sort(np.where(ok[:, None], z, -np.inf), axis=1)
        margin = srt[:, -1] - srt[:, -2]
        m = test.to_numpy() & ok
        acc = float((pred[m] == d.rao.to_numpy()[m]).mean())
        ab = lambda x: np.where(np.isin(x, ["A1", "A2"]), "A", "B")
        acc2 = float((ab(pred[m]) == ab(d.rao.to_numpy()[m])).mean())
        conf = pd.crosstab(pd.Series(d.rao.to_numpy()[m], name="rao"), pd.Series(pred[m], name="pred"), normalize="index")
        per_round.append({"round": rnd, "five_state_acc": acc, "two_state_acc": acc2,
                          "confusion": conf.reindex(index=ms.STATES, columns=ms.STATES).fillna(0).round(3).to_dict(orient="index")})
        print(f"round {rnd}: 5-state {acc:.3f}  A/B {acc2:.3f}", flush=True)
        new = pd.Series("", index=d.index, dtype=object)
        for st in ms.STATES:
            idx = np.flatnonzero(ok & (pred == st))
            if idx.size:
                keep = idx[margin[idx] >= np.median(margin[idx])]
                new.iloc[keep] = st
        seeds = new
    freq = d.loc[test, "rao"].value_counts(normalize=True)
    res["chance_five_state"] = float((freq ** 2).sum())
    res["rounds"] = per_round
    res["n_test_bins"] = int(test.sum())
    (RESULTS / "multistate_annot.json").write_text(json.dumps(res, indent=1, default=float))

    best = per_round[-1]
    cm = pd.DataFrame(best["confusion"]).T.reindex(index=ms.STATES, columns=ms.STATES).to_numpy()
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.2))
    im = ax[0].imshow(cm, cmap="Reds", vmin=0, vmax=1)
    for i in range(5):
        for j in range(5):
            ax[0].text(j, i, f"{cm[i, j]:.2f}", ha="center", va="center", fontsize=7,
                       color="white" if cm[i, j] > 0.5 else "black")
    ax[0].set_xticks(range(5))
    ax[0].set_yticks(range(5))
    ax[0].set_xticklabels(ms.STATES)
    ax[0].set_yticklabels(ms.STATES)
    ax[0].set_xlabel("annotation-seeded multi-state CRUSH")
    ax[0].set_ylabel("Rao 2014 subcompartment")
    cb = fig.colorbar(im, ax=ax[0], fraction=0.046)
    cb.set_label("Fraction of bins per\nRao subcompartment", fontsize=7)
    ax[1].plot([r["round"] for r in per_round], [r["five_state_acc"] for r in per_round], "o-", color="#b2182b", label="5-state")
    ax[1].plot([r["round"] for r in per_round], [r["two_state_acc"] for r in per_round], "o-", color="#2166ac", label="A/B")
    ax[1].axhline(res["chance_five_state"], ls="--", color="grey", lw=0.8, label="5-state chance")
    ax[1].set_xlabel("re-seeding round")
    ax[1].set_ylabel("accuracy vs Rao (non-seed bins)")
    ax[1].set_ylim(0, 1)
    ax[1].legend(frameon=False, fontsize=7)
    for s in ("top", "right"):
        ax[1].spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIGS / f"multistate_annot.{ext}", dpi=300, bbox_inches="tight")
    print(json.dumps({k: res[k] for k in ("seed_counts", "seed_purity_vs_rao", "chance_five_state", "n_test_bins")}, indent=1))


if __name__ == "__main__":
    main()
