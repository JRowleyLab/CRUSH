#!/usr/bin/env python3
import json
import numpy as np
import pandas as pd
import multistate_annot as ma
import multistate_crush as ms

sizes = {l.split()[0]: int(l.split()[1]) for l in open(ms.SIZES)}
lab = ms.labels_at(ma.R, sizes)
rows = []
for c in ms.CHROMS:
    n = sizes[c] // ma.R + 1
    df = pd.DataFrame({k: ma.bin_signal(ma.MARKS / f, c, n) for k, f in ma.TRACKS.items()})
    df["genes"], df["lowgc"] = ma.density(ma.GENES, c, n), ma.density(ma.BBINS, c, n)
    df["rao"] = lab[c]
    rows.append(df)
d = pd.concat(rows, ignore_index=True)
for k in ma.TRACKS:
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
seed = pd.Series("", index=d.index, dtype=object)
for st, m in rules.items():
    seed[m & (hits == 1)] = st
feats = ["act", "K27me3", "K9me3", "HP1", "genes", "lowgc"]
X = d[feats].copy()
X = (X - X.mean()) / X.std()
cent = np.vstack([X[seed == st].mean().to_numpy() for st in ms.STATES])
test = (d.rao != "") & (seed == "") & X.notna().all(1)
dist = ((X[test].to_numpy()[:, None, :] - cent[None]) ** 2).sum(-1)
pred = np.array(ms.STATES)[dist.argmin(1)]
truth = d.loc[test, "rao"].to_numpy()
ab = lambda x: np.where(np.isin(x, ["A1", "A2"]), "A", "B")
out = {"n_test": int(test.sum()), "five_state_acc": float((pred == truth).mean()),
       "two_state_acc": float((ab(pred) == ab(truth)).mean()),
       "per_state_recall": {s: float((pred[truth == s] == s).mean()) for s in ms.STATES}}
print(json.dumps(out, indent=1))
open(ma.RESULTS / "marks_only_control.json", "w").write(json.dumps(out, indent=1))
