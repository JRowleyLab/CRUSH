#!/usr/bin/env bash
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
DATA=$(cd "$HERE/.." && pwd)/init_robustness
CRUSH=${CRUSH:-crush-rs}
HIC=${HIC:-full_remake.hic}
A=hg19_nochr
SIZES=$A/hg19.sizes
GENES=$A/hg19genes.bed
FASTA=$A/hg19_nochr.fa
S=$DATA/seeds
CPU=${CPU:-22}

for f in "$CRUSH" "$HIC" "$SIZES" "$GENES" "$FASTA" \
         "$S/gc_high_A.bed" "$S/gc_low_B.bed" "$S/rand_A.bed" "$S/rand_B.bed"; do
    [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

run() {
    local arm=$1 a=$2 b=$3 out="$DATA/$1"
    if [ -s "$out/${arm}_mergedCrush_1000.bedgraph" ]; then
        echo "skip (done): $arm"; return
    fi
    mkdir -p "$out"
    echo "[$(date +%H:%M:%S)] start $arm"
    "$CRUSH" -i "$HIC" -r 1000 -g "$SIZES" -a "$a" -b "$b" \
        -o "${arm}_" -O "$out" -c "$CPU" > "$out/run.log" 2>&1 \
        && echo "[$(date +%H:%M:%S)] done  $arm" \
        || echo "[$(date +%H:%M:%S)] FAILED $arm" >&2
}

run std    "$GENES"          "$FASTA"          &
run gconly "$S/gc_high_A.bed" "$S/gc_low_B.bed" &
run rand   "$S/rand_A.bed"    "$S/rand_B.bed"   &
wait
echo "all arms finished"
