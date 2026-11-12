#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = "validated_unique_enhancers.tsv"
POOL = "iqr_normalized.tsv"
LINES = ["A549", "HCT", "HELA", "HEPG2", "HUVEC", "IMR90", "K562"]
LABEL = dict(zip(LINES, ["A549", "HCT-116", "HeLa-S3", "HepG2", "HUVEC", "IMR90", "K562"]))
COL = dict(zip(LINES, ["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC", "#CA9161", "#FBAFE4"]))
N = 10_000


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def main():
    d = pd.read_csv(SRC, sep="\t")
    S = d[[f"{l}_score" for l in LINES]].to_numpy(float)
    med = np.nanmedian(S, axis=1)
    lab = d.enhancer_called_in.map({l: i for i, l in enumerate(LINES)}).to_numpy()
    idx = np.arange(len(S))
    rng = np.random.default_rng(1)

    def stats(labels):
        v = S[idx, labels] - med
        return np.nanmean(v), np.array([np.nanmean(v[labels == i]) for i in range(7)])

    obs_all, obs_line = stats(lab)
    shuf_all, shuf_line = np.empty(N), np.empty((N, 7))
    for k in range(N):
        shuf_all[k], shuf_line[k] = stats(rng.permutation(lab))

    P = pd.read_csv(POOL, sep="\t", header=None, usecols=range(3, 10)).to_numpy(float)
    Pd = P - np.nanmedian(P, axis=1, keepdims=True)
    nl = np.bincount(lab, minlength=7)
    null_all, null_line = np.empty(N), np.empty((N, 7))
    for k in range(N):
        sums = []
        for i in range(7):
            x = Pd[rng.choice(len(Pd), nl[i], replace=False), i]
            null_line[k, i] = np.nanmean(x)
            sums.append(np.nansum(x))
        null_all[k] = np.sum(sums) / nl.sum()

    v = S[idx, lab] - med
    boot = np.array([np.nanmean(v[rng.integers(0, len(v), len(v))]) for _ in range(N)])
    boot_line = np.array([[np.nanmean(w[rng.integers(0, len(w), len(w))]) for _ in range(N)]
                          for w in (v[lab == i] for i in range(7))])
    p_all = ((null_all >= obs_all).sum() + 1) / (N + 1)
    p_line = ((null_line >= obs_line).sum(0) + 1) / (N + 1)
    ps_all = ((shuf_all >= obs_all).sum() + 1) / (N + 1)
    ps_line = ((shuf_line >= obs_line).sum(0) + 1) / (N + 1)
    null_ln = null_line.mean(0)
    obs_cl, null_cl = obs_line.mean(), null_ln.mean()
    p_cl = ((null_line.mean(1) >= obs_cl).sum() + 1) / (N + 1)

    rows = [{"set": "mean of seven lines (4B)", "n": 7, "observed": obs_cl, "observed_sd": obs_line.std(ddof=1),
             "null_mean": null_cl, "null_sd": null_ln.std(ddof=1), "p_empirical": p_cl},
            {"set": "all unique enhancers, pooled", "n": len(v), "observed": obs_all,
             "ci95_lo": np.percentile(boot, 2.5), "ci95_hi": np.percentile(boot, 97.5),
             "null_mean": null_all.mean(), "null_q025": np.percentile(null_all, 2.5),
             "null_q975": np.percentile(null_all, 97.5), "p_empirical": p_all, "p_label_shuffle": ps_all}]
    for i, l in enumerate(LINES):
        rows.append({"set": LABEL[l], "n": int((lab == i).sum()), "observed": obs_line[i],
                     "ci95_lo": np.percentile(boot_line[i], 2.5), "ci95_hi": np.percentile(boot_line[i], 97.5),
                     "null_mean": null_line[:, i].mean(), "null_q025": np.percentile(null_line[:, i], 2.5),
                     "null_q975": np.percentile(null_line[:, i], 97.5), "p_empirical": p_line[i],
                     "p_label_shuffle": ps_line[i]})
    st = pd.DataFrame(rows)
    st.to_csv(HERE / "fig4B_S9A_stats.tsv", sep="\t", index=False, float_format="%.4g")
    print(st.round(4).to_string(index=False))

    for notext in (False, True):
        fig, ax = plt.subplots(figsize=(2.7, 3.0))
        for x, m, vals, c in ((0, obs_cl, obs_line, "#00e5ee"), (1, null_cl, null_ln, "#808080")):
            ax.bar(x, m, width=0.6, color=c, edgecolor="black", lw=0.8)
            ax.errorbar(x, m, yerr=vals.std(ddof=1), color="black", capsize=4, lw=0.8)
            jit = np.random.default_rng(3).uniform(-0.15, 0.15, 7)
            for i, l in enumerate(LINES):
                ax.scatter(x + jit[i], vals[i], s=14, color=COL[l], edgecolor="black", lw=0.3, zorder=3)
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks([0, 1])
        ax.set_xlim(-0.6, 1.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=8)
        if notext:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(["Unique\nenhancers", "Random\nenhancers"], fontsize=8)
            ax.set_ylabel("ΔCRUSH [IQR]\n(unique line - median)", fontsize=9)
            ax.text(0.5, max(obs_line.max(), obs_cl + obs_line.std(ddof=1)) * 1.03, stars(p_cl),
                    ha="center", fontsize=10)
        fig.tight_layout()
        stem = "fig4B" + ("_notext" if notext else "")
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)

        fig, axs = plt.subplots(1, 7, figsize=(8.4, 2.6), sharey=True)
        for i, (ax, l) in enumerate(zip(axs, LINES)):
            ax.bar(0, obs_line[i], width=0.6, color="#00e5ee", edgecolor="black", lw=0.8)
            nm = null_line[:, i].mean()
            nlo, nhi = np.percentile(null_line[:, i], [2.5, 97.5])
            ax.bar(1, nm, width=0.6, color="#808080", edgecolor="black", lw=0.8)
            ax.errorbar(1, nm, yerr=[[nm - nlo], [nhi - nm]], color="black", capsize=3, lw=0.8)
            ax.axhline(0, color="black", lw=0.6)
            ax.set_xticks([0, 1])
            ax.set_xlim(-0.6, 1.6)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.tick_params(labelsize=8)
            if notext:
                ax.set_xticklabels([])
            else:
                ax.set_xticklabels(["Unique", "Random"], fontsize=7, rotation=45, ha="right")
                ax.set_title(f"{LABEL[l]}\n(n={nl[i]:,})", fontsize=8)
                ax.text(0.5, max(obs_line[i], nhi) + 0.03, stars(p_line[i]), ha="center", fontsize=9)
        if not notext:
            axs[0].set_ylabel("ΔCRUSH [IQR]\n(unique line - median)", fontsize=9)
        fig.tight_layout()
        stem = "figS9A" + ("_notext" if notext else "")
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
