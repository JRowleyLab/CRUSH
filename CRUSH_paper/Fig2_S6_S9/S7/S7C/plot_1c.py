#!/usr/bin/env python3
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"
FIGS = Path(__file__).resolve().parent.parent / "figures"

S = ["A1", "A2", "B1", "B2", "B3"]
an = json.load(open(RESULTS / "multistate_annot.json"))
lb = json.load(open(RESULTS / "multistate_crush.json"))
mo = json.load(open(RESULTS / "marks_only_control.json"))
cm = pd.DataFrame(an["rounds"][0]["confusion"]).T.reindex(index=S, columns=S).to_numpy()
fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.3), gridspec_kw={"width_ratios": [1, 1.15]})
im = ax[0].imshow(cm, cmap="Reds", vmin=0, vmax=1)
for i in range(5):
    for j in range(5):
        ax[0].text(j, i, f"{cm[i, j]:.2f}", ha="center", va="center", fontsize=7,
                   color="white" if cm[i, j] > 0.5 else "black")
ax[0].set_xticks(range(5)); ax[0].set_yticks(range(5))
ax[0].set_xticklabels(S); ax[0].set_yticklabels(S)
ax[0].set_xlabel("multi-state CRUSH (annotation seeds)", fontsize=8)
ax[0].set_ylabel("Rao 2014 subcompartment", fontsize=8)
cb = fig.colorbar(im, ax=ax[0], fraction=0.046)
cb.set_label("Fraction of bins per\nRao subcompartment", fontsize=7)
cb.ax.tick_params(labelsize=7)
names = ["label seeds\n(upper bound)", "annotation\nseeds", "marks only\n(no contacts)", "chance"]
five = [lb["heldout_100kb"]["five_state_accuracy_mean"], an["rounds"][0]["five_state_acc"],
        mo["five_state_acc"], an["chance_five_state"]]
two = [lb["heldout_100kb"]["two_state_accuracy_mean"], an["rounds"][0]["two_state_acc"], mo["two_state_acc"], np.nan]
x = np.arange(4)
ax[1].bar(x - 0.18, five, 0.36, color="#b2182b", label="5-state")
ax[1].bar(x + 0.18, two, 0.36, color="#2166ac", label="A/B")
ax[1].set_xticks(x); ax[1].set_xticklabels(names, fontsize=7)
ax[1].set_ylabel("accuracy vs Rao 2014 (held-out bins)", fontsize=8)
ax[1].set_ylim(0, 1); ax[1].legend(frameon=False, fontsize=7)
for s in ("top", "right"):
    ax[1].spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(FIGS / f"panel_1c.{ext}", dpi=300, bbox_inches="tight")
