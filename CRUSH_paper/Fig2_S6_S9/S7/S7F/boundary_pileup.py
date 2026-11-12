#!/usr/bin/env python3
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse

import domain_ada as ada
import domain_classes as dc
import loop_vs_compartment_domains as lv
import pandas as pd

HERE = Path(__file__).resolve().parent
RES = 5000
W = 20


def load_band(chrom, n):
    cmd = f"{dc.STRAW} oe KR {dc.HIC} {chrom} {chrom} BP {RES} | awk '($2-$1) <= {2 * W * RES}'"
    txt = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True).stdout
    a = np.loadtxt(txt.splitlines(), ndmin=2) if txt else np.empty((0, 3))
    a = a[np.isfinite(a[:, 2]) & (a[:, 2] > 0)] if len(a) else a
    i, j = (a[:, 0] // RES).astype(int), (a[:, 1] // RES).astype(int)
    return sparse.coo_matrix((a[:, 2], (i, j)), shape=(n, n)).tocsr()


def one_chrom(chrom, n, site_sets):
    m = load_band(chrom, n)
    out = {}
    for key, sites in site_sets.items():
        acc, cnt, k = np.zeros((2 * W + 1,) * 2), np.zeros((2 * W + 1,) * 2), 0
        for p, flip in sites.get(chrom, []):
            if p - W < 0 or p + W >= n:
                continue
            sub = m[p - W:p + W + 1, p - W:p + W + 1].toarray()
            sub = np.triu(sub) + np.triu(sub, 1).T
            sub[sub == 0] = np.nan
            if flip:
                sub = sub[::-1, ::-1]
            ok = np.isfinite(sub)
            acc[ok] += sub[ok]
            cnt[ok] += 1
            k += 1
        out[key] = (acc, cnt, k)
    return out


def main():
    sizes = {l.split()[0]: int(l.split()[1]) for l in open(dc.SIZES)}
    lp = pd.read_csv(dc.LOOPS, sep="\t", usecols=["chr1", "x1", "x2", "y1", "y2"], dtype={"chr1": str})
    lp["chr1"] = lp["chr1"].str.replace("chr", "", regex=False)
    lp = lp[lp.chr1.isin(dc.CHROMS)]
    sets = {"loop domains": {}}
    for c, g in lp.groupby("chr1"):
        s, e = ((g.x1 + g.x2) // 2).to_numpy(), ((g.y1 + g.y2) // 2).to_numpy()
        keep = (e - s >= ada.MIN_L) & (e - s <= ada.MAX_L)
        sets["loop domains"][c] = [(x // RES, False) for x in s[keep]] + [(y // RES, True) for y in e[keep]]
    for name in ["CRUSH v2 1 kb", "eigenvector 1 kb", "eigenvector 25 kb"]:
        tr = dc.load_track(*dc.TRACKS[name], sizes)
        key = f"A compartment domains, {name}"
        sets[key] = {}
        for c, a in tr.items():
            d = [(s * dc.RES, e * dc.RES) for s, e, sg in lv.comp_domains(a) if sg > 0]
            d = [(s, e) for s, e in d if ada.MIN_L <= e - s <= ada.MAX_L]
            sets[key][c] = [(s // RES, False) for s, _ in d] + [(e // RES, True) for _, e in d]

    tot = {k: [np.zeros((2 * W + 1,) * 2), np.zeros((2 * W + 1,) * 2), 0] for k in sets}
    with ThreadPoolExecutor(6) as ex:
        for f in [ex.submit(one_chrom, c, sizes[c] // RES + 1, sets) for c in dc.CHROMS]:
            for k, (acc, cnt, n) in f.result().items():
                tot[k][0] += acc
                tot[k][1] += cnt
                tot[k][2] += n
    mats = {k: np.log2(v[0] / np.where(v[1] > 0, v[1], np.nan)) for k, v in tot.items()}
    idx = np.arange(2 * W + 1)
    for mtx in mats.values():
        mtx[np.abs(idx[:, None] - idx[None, :]) < 1] = np.nan
    res = {}
    for k, mtx in mats.items():
        up = mtx[:W, :W][np.triu_indices(W, 1)]
        down = mtx[W + 1:, W + 1:][np.triu_indices(W, 1)]
        cross = mtx[:W, W + 1:].ravel()
        res[k] = {"n_sites": int(tot[k][2]), "upstream": float(np.nanmean(up)),
                  "downstream_domain": float(np.nanmean(down)), "across": float(np.nanmean(cross))}
        print(k, res[k], flush=True)
    (HERE / "boundary_pileup.json").write_text(json.dumps(res, indent=1))

    keys = list(sets)
    lim = np.nanpercentile(np.abs(np.concatenate([m.ravel() for m in mats.values()])), 98)
    for notext in (False, True):
        fig, ax = plt.subplots(1, len(keys), figsize=(3.0 * len(keys), 3.2))
        for a, k in zip(ax, keys):
            im = a.imshow(mats[k], cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
            a.axhline(W, color="black", lw=0.4, ls=":")
            a.axvline(W, color="black", lw=0.4, ls=":")
            a.set_xticks([])
            a.set_yticks([])
            if not notext:
                a.set_title(f"{k.replace('A compartment domains, ', 'A comp. domain edges\n')}\n(n = {res[k]['n_sites']:,})", fontsize=7)
        cb = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01)
        cb.ax.tick_params(labelsize=7)
        if not notext:
            cb.set_label("log2(O/E)", fontsize=8)
        stem = "boundary_pileup" + ("_notext" if notext else "")
        for ext in ("png", "pdf", "svg"):
            fig.savefig(HERE / f"{stem}.{ext}", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
