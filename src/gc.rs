//! GC-content B-compartment seeding from a genome FASTA.
//!
//! Mirrors the bash pipeline's `-b <FASTA>` path, which shells out to
//! `bedtools nuc` over fixed-width bins and then keeps the low-GC tail:
//!
//! ```text
//!   bins      : 500 bp (or the finest resolution, if that is smaller)
//!   gc(bin)   : (C + G) / (A + C + G + T)      -- N excluded from the divisor
//!   threshold : mean(gc) - 1 * sd(gc)
//!   Bbins     : bins with gc < threshold
//! ```
//!
//! Bins that are entirely N (assembly gaps, telomeres, centromere) have no
//! callable bases and are dropped rather than counted as GC-poor, matching the
//! `$6 > 0` guard in the bash version.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::bed::{BedSet, ChromSizes, Interval};

/// True if the file starts with a FASTA header line.
pub fn looks_like_fasta<P: AsRef<Path>>(path: P) -> bool {
    let Ok(f) = File::open(path) else {
        return false;
    };
    let mut first = String::new();
    if BufReader::new(f).read_line(&mut first).is_err() {
        return false;
    }
    first.starts_with('>')
}

/// Compute low-GC bins from a FASTA, restricted to chromosomes in `sizes`.
pub fn bbins_from_fasta<P: AsRef<Path>>(
    path: P,
    sizes: &ChromSizes,
    bin: u64,
) -> io::Result<BedSet> {
    let f = BufReader::new(File::open(&path)?);

    // (chrom, bin_index) -> (gc_count, at_count)
    let mut counts: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut chrom = String::new();
    let mut pos: u64 = 0;

    for line in f.lines() {
        let line = line?;
        let l = line.trim_end();
        if l.is_empty() {
            continue;
        }
        if let Some(rest) = l.strip_prefix('>') {
            // FASTA headers may carry a description after the name.
            chrom = rest.split_whitespace().next().unwrap_or("").to_string();
            pos = 0;
            if let Some(len) = sizes.length_of(&chrom) {
                let n = (len / bin) as usize + 1;
                counts.insert(chrom.clone(), vec![(0, 0); n]);
            }
            continue;
        }
        let Some(v) = counts.get_mut(&chrom) else {
            // Chromosome not in the sizes file: skip its sequence entirely.
            continue;
        };
        for b in l.bytes() {
            let idx = (pos / bin) as usize;
            if idx < v.len() {
                match b | 0x20 {
                    b'c' | b'g' => v[idx].0 += 1,
                    b'a' | b't' => v[idx].1 += 1,
                    _ => {} // N and anything else: not callable
                }
            }
            pos += 1;
        }
    }

    // Genome-wide mean and SD of the per-bin GC fraction.
    let mut fracs: Vec<f64> = Vec::new();
    for v in counts.values() {
        for &(gc, at) in v {
            let tot = gc as f64 + at as f64;
            if tot > 0.0 {
                fracs.push(gc as f64 / tot);
            }
        }
    }
    if fracs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: no callable bases found; are the FASTA sequence names the \
                 same as in the chromosome sizes file?",
                path.as_ref().display()
            ),
        ));
    }

    let n = fracs.len() as f64;
    let mean = fracs.iter().sum::<f64>() / n;
    let var = fracs.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / n;
    let threshold = mean - var.sqrt();

    let mut by_chrom: HashMap<String, Vec<Interval>> = HashMap::new();
    for (c, v) in &counts {
        let mut out = Vec::new();
        for (i, &(gc, at)) in v.iter().enumerate() {
            let tot = gc as f64 + at as f64;
            if tot > 0.0 && (gc as f64 / tot) < threshold {
                let s = i as u64 * bin;
                out.push(Interval {
                    start: s,
                    end: s + bin,
                });
            }
        }
        if !out.is_empty() {
            by_chrom.insert(c.clone(), out);
        }
    }

    eprintln!(
        "      GC seeding: {} bins scored, mean={:.4} sd={:.4}, threshold={:.4}, {} B bins",
        fracs.len(),
        mean,
        var.sqrt(),
        threshold,
        by_chrom.values().map(|v| v.len()).sum::<usize>()
    );

    Ok(BedSet { by_chrom })
}
