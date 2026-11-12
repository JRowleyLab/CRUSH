#!/bin/bash
# CRUSH demo: chr20 of the 100M-contact subsample of the LCL Hi-C map, called at 1 kb.
# Run from the demo/ folder. Uses crush-rs on your PATH, or the repository's
# prebuilt binary if there is none.
set -euo pipefail
cd "$(dirname "$0")"
CRUSH=$(command -v crush-rs || echo ../crush-rs)

$CRUSH -i LCL_100M_chr20.hic \
       -g hg19_chr20.sizes \
       -a hg19_chr20_genes.bed \
       -b hg19_chr20_lowGC.bed \
       -r 1000 -c 4 -o demo_

# Check the result against the expected track
if gunzip -c expected/demo_mergedCrush_1000.bedgraph.gz | cmp -s - CRUSHtracks/demo_mergedCrush_1000.bedgraph; then
    echo "Demo output matches expected/demo_mergedCrush_1000.bedgraph.gz"
else
    echo "Demo output differs from expected/demo_mergedCrush_1000.bedgraph.gz" >&2
    exit 1
fi
