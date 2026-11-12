#!/usr/bin/env bash
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA=$(cd "$HERE/.." && pwd)/init_robustness
CRUSH=${CRUSH:-crush-rs}
HIC=${HIC:-full_remake.hic}
A=hg19_nochr
SIZES=$A/hg19.sizes
CPU=${CPU:-14}
JOBS=${JOBS:-2}

python3 - <<'PY'
import numpy as np, pandas as pd
from pathlib import Path
S=Path("seeds")
pool=pd.concat([pd.read_csv(S/f,sep="\t",header=None,names=["c","s","e"],dtype={"c":str})
                for f in ("gc_high_A.bed","gc_low_B.bed")],ignore_index=True)
pool=pool.drop_duplicates().sort_values(["c","s"]).reset_index(drop=True)
nA,nB=56779,56863
for k in range(2,6):
    rng=np.random.default_rng(1000+k)
    idx=rng.permutation(len(pool))
    a=pool.iloc[np.sort(idx[:nA])]; b=pool.iloc[np.sort(idx[nA:nA+nB])]
    a.to_csv(S/f"rand{k}_A.bed",sep="\t",header=False,index=False)
    b.to_csv(S/f"rand{k}_B.bed",sep="\t",header=False,index=False)
    print(f"  rand{k}: {len(a):,} A / {len(b):,} B windows (disjoint)")
PY

run() {
    local arm=$1 out="$DATA/$1"
    if [ -s "$out/${arm}_mergedCrush_1000.bedgraph" ]; then echo "skip $arm"; return; fi
    mkdir -p "$out"
    echo "[$(date +%H:%M:%S)] start $arm"
    "$CRUSH" -i "$HIC" -r 1000 -g "$SIZES" \
        -a "$DATA/seeds/${arm}_A.bed" -b "$DATA/seeds/${arm}_B.bed" \
        -o "${arm}_" -O "$out" -c "$CPU" > "$out/run.log" 2>&1 \
        && echo "[$(date +%H:%M:%S)] done  $arm" || echo "[$(date +%H:%M:%S)] FAILED $arm" >&2
}
export -f run; export CRUSH HIC SIZES CPU DATA
printf "%s\n" rand2 rand3 rand4 rand5 | xargs -P "$JOBS" -I{} bash -c 'run "$@"' _ {}
echo "all random arms finished"
