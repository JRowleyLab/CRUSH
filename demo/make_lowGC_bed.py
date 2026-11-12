#!/usr/bin/env python3
"""Low-GC B-seed BED from a FASTA, as crush-rs derives it for -b FASTA
(src/gc.rs): 500 bp bins, GC = (C+G)/(A+C+G+T), threshold mean - 1 SD over
bins with callable bases on the chromosomes in the sizes file; keeps bins
below it. Usage: lowgc_bed.py genome.fa sizes > lowgc.bed
"""
import sys
import numpy as np

BIN = 500
fa, sizes = sys.argv[1], sys.argv[2]
keep = {l.split()[0]: int(l.split()[1]) for l in open(sizes) if l.strip()}
seqs, name, buf = {}, None, []
for line in open(fa):
    if line.startswith(">"):
        if name in keep:
            seqs[name] = "".join(buf)
        name, buf = line[1:].split()[0], []
    elif name in keep:
        buf.append(line.strip())
if name in keep:
    seqs[name] = "".join(buf)
out, fr = {}, []
for c, s in seqs.items():
    a = np.frombuffer(s.encode(), dtype=np.uint8) | 0x20
    n = keep[c] // BIN + 1
    idx = np.arange(len(a)) // BIN
    ok = idx < n
    gc = np.bincount(idx[ok], weights=np.isin(a[ok], (ord("c"), ord("g"))), minlength=n)
    at = np.bincount(idx[ok], weights=np.isin(a[ok], (ord("a"), ord("t"))), minlength=n)
    tot = gc + at
    out[c] = (gc, tot)
    fr.append(gc[tot > 0] / tot[tot > 0])
fr = np.concatenate(fr)
thr = fr.mean() - fr.std()
for c in keep:
    if c not in out:
        continue
    gc, tot = out[c]
    for i in np.flatnonzero((tot > 0) & (gc / np.where(tot > 0, tot, 1) < thr)):
        print(f"{c}\t{i * BIN}\t{i * BIN + BIN}")
print(f"threshold {thr:.4f}", file=sys.stderr)
