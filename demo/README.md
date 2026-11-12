# CRUSH demo: LCL chr20, 100 million contacts

This folder holds a small dataset for checking a CRUSH installation and seeing what a run looks like. It runs in about 10 seconds on 4 threads and uses about 100 MB of memory.

## Files

| File | Contents |
|---|---|
| `LCL_100M_chr20.hic` | Chromosome 20 (2.33 M contacts) of a 100 M-contact random subsample of the 20.3 B-contact lymphoblastoid (LCL) Hi-C map used in Fig. 1 of the CRUSH paper. It has resolutions from 1 Mb to 1 kb, no normalisation, and hg19 coordinates without a `chr` prefix. |
| `hg19_chr20.sizes` | Chromosome sizes (chr20 only) |
| `hg19_chr20_genes.bed` | A seeds: hg19 genes on chr20 |
| `hg19_chr20_lowGC.bed` | B seeds: low-GC 500 bp bins on chr20, the same bins CRUSH picks when `-b` is the hg19 FASTA |
| `make_lowGC_bed.py` | The script that made `hg19_chr20_lowGC.bed` from the hg19 FASTA |
| `run_demo.sh` | Runs the demo and checks the result |
| `expected/demo_mergedCrush_1000.bedgraph.gz` | The output you should get |

## Run

```
cd demo
./run_demo.sh
```

This is equivalent to:

```
crush-rs -i LCL_100M_chr20.hic -g hg19_chr20.sizes -a hg19_chr20_genes.bed -b hg19_chr20_lowGC.bed -r 1000 -c 4 -o demo_
```

The script uses `crush-rs` from your `PATH`, or the prebuilt `../crush-rs` if there is none. It ends with:

```
  compartment check
      (AA+BB)/2AB over quartiles of the merged track, scored at 100 kb,
      contacts >= 1 Mb apart; 1.0 = the track shifted along itself
        20      2.141   null 0.982   z    56.9       599 bins
      median over 1 chromosomes: 2.141  [supported]

Finished.
Demo output matches expected/demo_mergedCrush_1000.bedgraph.gz
```

## Output

`CRUSHtracks/demo_mergedCrush_1000.bedgraph` holds the 1 kb compartment track for chr20 (A positive, B negative). If UCSC `bedGraphToBigWig` is installed, you also get `demo_mergedCrush_1000.bw`. The bigWig uses `chr20` as the chromosome name so that it loads directly against IGV's hg19 genome.

The track is identical to chr20 of a CRUSH run on the full 100 M-contact map when that run is given the same chr20 size file.

To see the per-resolution tracks, add `--keepres all`:

```
crush-rs -i LCL_100M_chr20.hic -g hg19_chr20.sizes -a hg19_chr20_genes.bed -b hg19_chr20_lowGC.bed -r 1000 -c 4 -o demo_ --keepres all
```

## How the `.hic` file was made

```
java -jar juicer_tools_1.22.01.jar dump observed NONE 100M.hic 20 20 BP 1000 chr20_1kb.txt
awk '{print 0,"20",$1,0,0,"20",$2,1,$3}' chr20_1kb.txt > chr20.short
java -jar juicer_tools_1.22.01.jar pre -n -r 1000000,500000,250000,100000,50000,25000,10000,5000,1000 \
     chr20.short LCL_100M_chr20.hic hg19_chr20.sizes
```
