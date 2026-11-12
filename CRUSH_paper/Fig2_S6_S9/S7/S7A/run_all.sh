#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$HERE/results/full}"

compartment-benchmark \
  --samples  "$HERE/config/samples.tsv" \
  --outdir   "$OUT" \
  --binsize  5000 \
  --baseline crush_v1 \
  --chromhmm-mode bp

echo
echo "Headline comparison:"
column -t -s $'\t' "$OUT/summary/active_separation_by_callset.tsv" 2>/dev/null || true
