//! Chromosome sizes, BED interval loading, and bin-occupancy counting.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct ChromSizes {
    pub names: Vec<String>,
    pub lengths: Vec<u64>,
    index: HashMap<String, usize>,
}

impl ChromSizes {
    pub fn load<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let f = BufReader::new(File::open(&path)?);
        let mut names = Vec::new();
        let mut lengths = Vec::new();
        let mut index = HashMap::new();

        for (lineno, line) in f.lines().enumerate() {
            let line = line?;
            let t = line.trim_end();
            if t.is_empty() {
                continue;
            }
            let cols: Vec<&str> = t.split('\t').collect();
            if cols.len() != 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "{}:{}: chromosome sizes file must have exactly two \
                         tab-separated columns (name, length); found {} column(s)",
                        path.as_ref().display(),
                        lineno + 1,
                        cols.len()
                    ),
                ));
            }
            let len: u64 = cols[1].trim().parse().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "{}:{}: chromosome length '{}' is not an integer",
                        path.as_ref().display(),
                        lineno + 1,
                        cols[1]
                    ),
                )
            })?;
            index.insert(cols[0].to_string(), names.len());
            names.push(cols[0].to_string());
            lengths.push(len);
        }

        if names.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{}: no chromosomes found", path.as_ref().display()),
            ));
        }

        Ok(ChromSizes {
            names,
            lengths,
            index,
        })
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn length_of(&self, name: &str) -> Option<u64> {
        self.index.get(name).map(|&i| self.lengths[i])
    }
}

/// A half-open genomic interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

/// BED intervals grouped by chromosome. Columns beyond the first three are
/// discarded, matching the `cut -f 1-3` the bash pipeline applies.
#[derive(Debug, Clone, Default)]
pub struct BedSet {
    pub by_chrom: HashMap<String, Vec<Interval>>,
}

impl BedSet {
    pub fn load<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let f = BufReader::new(File::open(&path)?);
        let mut by_chrom: HashMap<String, Vec<Interval>> = HashMap::new();

        for line in f.lines() {
            let line = line?;
            let t = line.trim_end();
            if t.is_empty() || t.starts_with('#') || t.starts_with("track") {
                continue;
            }
            let mut it = t.split('\t');
            let (Some(chrom), Some(s), Some(e)) = (it.next(), it.next(), it.next()) else {
                continue;
            };
            let (Ok(start), Ok(end)) = (s.trim().parse::<u64>(), e.trim().parse::<u64>()) else {
                continue;
            };
            by_chrom
                .entry(chrom.to_string())
                .or_default()
                .push(Interval { start, end });
        }

        Ok(BedSet { by_chrom })
    }

    pub fn get(&self, chrom: &str) -> &[Interval] {
        self.by_chrom.get(chrom).map(|v| &v[..]).unwrap_or(&[])
    }

    pub fn total(&self) -> usize {
        self.by_chrom.values().map(|v| v.len()).sum()
    }
}

/// Set of bin start coordinates (in bp) covered by any interval of `bed` on
/// `chrom`, at resolution `res`.
///
/// This reproduces the bash idiom
/// ```text
/// mawk '{if ($1==chr) print $1, int($2/r)*r, int($3/r)*r}' \
///   | awk '{for (i=$2; i<=$3; i+=r) b[i]+=1} END {for (j in b) print j, b[j]}'
/// ```
/// Note the `<=`: a feature ending exactly on a bin boundary also marks that
/// bin. Kept deliberately so the Rust and bash implementations agree.
pub fn occupied_bins(
    bed: &BedSet,
    chrom: &str,
    res: u64,
) -> HashMap<u64, u32> {
    let mut counts: HashMap<u64, u32> = HashMap::new();
    for iv in bed.get(chrom) {
        let s = (iv.start / res) * res;
        let e = (iv.end / res) * res;
        let mut i = s;
        while i <= e {
            *counts.entry(i).or_insert(0) += 1;
            i += res;
        }
    }
    counts
}
