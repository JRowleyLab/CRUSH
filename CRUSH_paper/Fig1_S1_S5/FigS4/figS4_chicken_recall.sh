#!/usr/bin/env bash
set -uo pipefail

CRUSH=${CRUSH:-crush-rs}
OUTROOT=${OUTROOT:-chicken}
G=ucsc
CHIC=hicfiles
CPU=${CPU:-24}

[ -x "$CRUSH" ] || { echo "no crush-rs at $CRUSH" >&2; exit 1; }

for S in Fibroblast Erythroblast eRBC aRBC; do
  for RES in 50000 1000; do
    OUT="$OUTROOT/${S}_r${RES}"
    if [ -s "$OUT/${S}_mergedCrush_${RES}.bedgraph" ]; then
      echo "  skip (done): ${S} @ ${RES}"; continue
    fi
    [ -s "$CHIC/GSE206192_${S}.hic" ] || { echo "  MISSING $CHIC/GSE206192_${S}.hic" >&2; continue; }
    mkdir -p "$OUT"
    echo "[$(date +%H:%M:%S)] start ${S} @ ${RES}"
    "$CRUSH" -i "$CHIC/GSE206192_${S}.hic" -r "$RES" \
        -g "$G/chicken_nochr.sizes" -a "$G/chicken_nochr_genes.bed" \
        -b "$G/chicken_nochr.fa" -o "${S}_" -O "$OUT" -c "$CPU" \
        > "$OUT/run.log" 2>&1
    if [ $? -eq 0 ] && [ -s "$OUT/${S}_mergedCrush_${RES}.bw" ]; then
      echo "[$(date +%H:%M:%S)] done  ${S} @ ${RES}"
    else
      echo "[$(date +%H:%M:%S)] FAILED ${S} @ ${RES} (see $OUT/run.log)" >&2
    fi
  done
done
echo "all chicken recalls finished"
