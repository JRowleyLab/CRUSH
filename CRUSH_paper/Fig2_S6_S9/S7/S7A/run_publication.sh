#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$HERE/results/publication}"

compartment-benchmark \
  --samples  "$HERE/config/samples_publication.tsv" \
  --outdir   "$OUT" \
  --binsize  5000 \
  --analysis vs_calder \
  --analysis chromhmm \
  --chromhmm-mode bp \
  --heatmap-scale sample

compartment-figures --results "$OUT"

echo
echo "Main figures:"
ls -1 "$OUT/figures/"
echo
echo "Headline comparison (active separation, percentage points):"
column -t -s $'\t' "$OUT/summary/active_separation_by_callset.tsv" 2>/dev/null || true
