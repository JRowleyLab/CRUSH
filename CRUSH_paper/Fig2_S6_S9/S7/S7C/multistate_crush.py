#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
FIGS = HERE.parent / "figures"
HIC = "full.hic"
STRAW = "straw"
SUBC = Path("GSE63525_GM12878_subcompartments.bed")
SIZES = "hg19.sizes"
CHROMS = [str(i) for i in range(1, 23)]
STATES = ["A1", "A2", "B1", "B2", "B3"]
NSPLIT = 10
MIN_SEP = 2


def oe_matrix(chrom: str, res: int, n: int) -> np.ndarray:
    out = subprocess.run([STRAW, "oe", "KR", HIC, chrom, chrom, "BP", str(res)],
                         capture_output=True, text=True, check=True).stdout
    a = np.loadtxt(out.splitlines(), dtype=float, ndmin=2) if out else np.empty((0, 3))
    m = np.full((n, n), np.nan)
    i, j, v = (a[:, 0] // res).astype(int), (a[:, 1] // res).astype(int), a[:, 2]
    ok = (i < n) & (j < n) & np.isfinite(v) & (v > 0)
    lv = np.log2(v[ok])
    m[i[ok], j[ok]] = lv
    m[j[ok], i[ok]] = lv
    covered = np.isfinite(m).sum(1) > 0
    floor = np.nanpercentile(lv, 1) if lv.size else 0.0
    sub = m[np.ix_(covered, covered)]
    sub[np.isnan(sub)] = floor
    m[np.ix_(covered, covered)] = sub
    idx = np.arange(n)
    m[np.abs(idx[:, None] - idx[None, :]) < MIN_SEP] = np.nan
    return m


def labels_at(res: int, sizes: dict) -> dict[str, np.ndarray]:
    sc = pd.read_csv(SUBC, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["c", "s", "e", "st"])
    sc["c"] = sc["c"].str.replace("chr", "", regex=False)
    out = {}
    for c in CHROMS:
        n = sizes[c] // res + 1
        lab = np.full(n, "", dtype=object)
        for s, e, st in sc.loc[sc["c"] == c, ["s", "e", "st"]].itertuples(index=False):
            if st in STATES:
                lab[s // res:min(n, e // res)] = st
        out[c] = lab
    return out


def state_scores(m: np.ndarray, seeds: dict[str, np.ndarray]) -> np.ndarray:
    """bins x states: mean log2 O/E to each state's seed bins (NaN if none)."""
    sc = np.full((m.shape[0], len(seeds)), np.nan)
    for k, idx in enumerate(seeds.values()):
        if idx.size:
            sc[:, k] = np.nanmean(m[:, idx], axis=1)
    return sc


def standardise(x: np.ndarray) -> np.ndarray:
    return (x - np.nanmean(x, axis=0)) / np.nanstd(x, axis=0)


def main() -> None:
    sizes = {l.split()[0]: int(l.split()[1]) for l in open(SIZES)}
    rng = np.random.default_rng(0)
    res = {}

    R = 100_000
    lab = labels_at(R, sizes)
    with ThreadPoolExecutor(8) as ex:
        mats = dict(zip(CHROMS, ex.map(lambda c: oe_matrix(c, R, sizes[c] // R + 1), CHROMS)))
    acc5, acc2, chance5, conf = [], [], [], np.zeros((5, 5))
    for split in range(NSPLIT):
        truth5, pred5, truth2, pred2 = [], [], [], []
        allsc5, allsc2, alltrue = [], [], []
        for c in CHROMS:
            m, L = mats[c], lab[c]
            labelled = np.flatnonzero(L != "")
            seed_mask = np.zeros(len(L), bool)
            seed_mask[rng.permutation(labelled)[: len(labelled) // 2]] = True
            test = np.flatnonzero((L != "") & ~seed_mask)
            s5 = {st: np.flatnonzero(seed_mask & (L == st)) for st in STATES}
            s2 = {"A": np.flatnonzero(seed_mask & np.isin(L, ["A1", "A2"])),
                  "B": np.flatnonzero(seed_mask & np.isin(L, ["B1", "B2", "B3"]))}
            allsc5.append(state_scores(m, s5)[test])
            allsc2.append(state_scores(m, s2)[test])
            alltrue.append(L[test])
        sc5, sc2, tr = standardise(np.vstack(allsc5)), standardise(np.vstack(allsc2)), np.concatenate(alltrue)
        ok = np.isfinite(sc5).all(1) & np.isfinite(sc2).all(1)
        p5 = np.array(STATES)[np.argmax(sc5[ok], axis=1)]
        p2 = np.where(sc2[ok, 0] > sc2[ok, 1], "A", "B")
        t5 = tr[ok]
        t2 = np.where(np.isin(t5, ["A1", "A2"]), "A", "B")
        acc5.append(float((p5 == t5).mean()))
        acc2.append(float((p2 == t2).mean()))
        freq = pd.Series(t5).value_counts(normalize=True)
        chance5.append(float((freq ** 2).sum()))
        for a, b in zip(t5, p5):
            conf[STATES.index(a), STATES.index(b)] += 1
    conf_row = conf / conf.sum(1, keepdims=True)
    res["heldout_100kb"] = {
        "five_state_accuracy_mean": float(np.mean(acc5)), "five_state_accuracy_sd": float(np.std(acc5)),
        "chance_five_state": float(np.mean(chance5)),
        "two_state_accuracy_mean": float(np.mean(acc2)),
        "confusion_rows_truth_cols_pred": {s: dict(zip(STATES, map(float, r))) for s, r in zip(STATES, conf_row)},
        "n_splits": NSPLIT,
    }

    r_f = 25_000
    lab_f = labels_at(r_f, sizes)
    diff, total, per_state = 0, 0, {s: [0, 0] for s in STATES}
    fine_sc, fine_true = [], []
    for c in CHROMS:
        mf = oe_matrix(c, r_f, sizes[c] // r_f + 1)
        Lf = lab_f[c]
        seeds = {st: np.flatnonzero(Lf == st) for st in STATES}
        fine_sc.append(state_scores(mf, seeds))
        fine_true.append(Lf)
        del mf
    sc = standardise(np.vstack(fine_sc))
    tr = np.concatenate(fine_true)
    ok = np.isfinite(sc).all(1) & (tr != "")
    pred = np.array(STATES)[np.argmax(sc[ok], axis=1)]
    t = tr[ok]
    for s in STATES:
        m = t == s
        per_state[s] = float((pred[m] != s).mean())
    res["refine_25kb"] = {
        "frac_25kb_bins_differing_from_100kb_label": float((pred != t).mean()),
        "per_state_frac_differing": per_state,
        "note": "seeds are all 25 kb bins inside each 100 kb label (in-sample); "
                "the held-out accuracy above is the generalisation test",
    }
    (RESULTS / "multistate_crush.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))

    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.2))
    im = ax[0].imshow(conf_row, cmap="Reds", vmin=0, vmax=1)
    ax[0].set_xticks(range(5))
    ax[0].set_yticks(range(5))
    ax[0].set_xticklabels(STATES)
    ax[0].set_yticklabels(STATES)
    ax[0].set_xlabel("multi-state CRUSH call")
    ax[0].set_ylabel("Rao 2014 subcompartment")
    for i in range(5):
        for j in range(5):
            ax[0].text(j, i, f"{conf_row[i, j]:.2f}", ha="center", va="center", fontsize=7,
                       color="white" if conf_row[i, j] > 0.5 else "black")
    cb = fig.colorbar(im, ax=ax[0], fraction=0.046)
    cb.set_label("Fraction of bins per\nRao subcompartment", fontsize=7)
    h = res["heldout_100kb"]
    ax[1].bar([0, 1, 2], [h["five_state_accuracy_mean"], h["chance_five_state"], h["two_state_accuracy_mean"]],
              yerr=[h["five_state_accuracy_sd"], 0, 0], color=["#b2182b", "#bbbbbb", "#2166ac"])
    ax[1].set_xticks([0, 1, 2])
    ax[1].set_xticklabels(["5-state", "chance\n(5-state)", "A/B only"], fontsize=8)
    ax[1].set_ylabel("held-out accuracy (100 kb)")
    ax[1].set_ylim(0, 1)
    for a in ax[1:]:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIGS / f"multistate_crush.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
