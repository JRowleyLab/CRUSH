//! `--v1`: a line-for-line port of the bash pipeline that produced the
//! published ENCODE, TNBC and HCT-116 tracks
//! (Crush_FV_OnePunch_2_linkadded_1_Use.bash, "CRUSH_v1.0").
//!
//! This is deliberately NOT the corrected algorithm in main.rs. It reproduces
//! v1 exactly, bugs included, so published numbers can be regenerated from
//! Rust. Each function names the bash function it mirrors. Behaviours that
//! differ from the corrected pipeline, all kept on purpose:
//!
//!   * score per bin = (ΣzA·nA − ΣzB·nB) / (total seed count); a bin seen
//!     only by B seeds gets +ΣzB·nB, not −ΣzB·nB          (oppocheck mawk)
//!   * z is normalised by the PARTNER bin's mean and SD    (RowZs)
//!   * A and B seed sets may overlap; a seed interval ending on a bin edge
//!     also seeds the next bin                              (process_genes_Bbins)
//!   * re-initialised seeds replace the original set per chromosome
//!                                                          (separating_ABbins)
//!   * the shifter never reads its coarse input: a pure high-pass filter
//!     over every bin, including empty ones                 (run_midshifter, B4)
//!   * the reprocess pass uses the sliding window of the CURRENT resolution
//!     for every coarser rung                               ($newwindow is global)
//!   * the merge divides each rung's score by its bin width before averaging
//!   * awk's %.6g number formatting is applied wherever bash writes a value
//!     through awk and reads it back
//!
//! `--eigen-init 0` reproduces the older v1 runs (HCT-116, Oct 2024) that had
//! no eigenvector step and seeded round 0 from the -a/-b annotations.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use nalgebra::DMatrix;
use rayon::prelude::*;

use crate::bed::ChromSizes;
use crate::hic::{HicFile, ValueKind};

// ---------------------------------------------------------------------------
// awk number formatting
// ---------------------------------------------------------------------------

/// The value awk would read back after printing `x` with the default
/// OFMT/CONVFMT "%.6g": integral values print exactly, everything else is
/// rounded to six significant digits.
#[inline]
pub fn awk(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    if x == x.trunc() && x.abs() < 9.0e15 {
        return x;
    }
    format!("{:.5e}", x).parse().unwrap_or(x)
}

// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------

/// Intervals per chromosome, in file order.
pub type Ivs = HashMap<String, Vec<(i64, i64)>>;

pub fn load_ivs(path: &Path) -> Result<Ivs, String> {
    let f = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut out: Ivs = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line.map_err(|e| e.to_string())?;
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 3 || line.starts_with("track") || line.starts_with('#') {
            continue;
        }
        let (Ok(s), Ok(e)) = (c[1].trim().parse::<f64>(), c[2].trim().parse::<f64>()) else {
            continue;
        };
        out.entry(c[0].to_string()).or_default().push((s as i64, e as i64));
    }
    Ok(out)
}

/// Low-GC B seeds from a FASTA, as the bash v1 block does it: 500 bp bins from
/// 0 to the chromosome size (bedtools nuc skips the partial bin past the end),
/// GC/(GC+AT) from bedtools' pct columns read back through awk's %.6g, and a
/// genome-wide threshold of mean - 1 SD, itself printed through %.6g.
pub fn bbins_from_fasta(path: &Path, sizes: &ChromSizes) -> Result<(Ivs, f64), String> {
    const BIN: usize = 500;
    let f = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    // chrom -> per-bin (gc, at)
    let mut counts: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut cur: Option<String> = None;
    let mut pos = 0usize;
    for line in BufReader::new(f).lines() {
        let line = line.map_err(|e| e.to_string())?;
        if let Some(h) = line.strip_prefix('>') {
            let name = h.split_whitespace().next().unwrap_or("").to_string();
            cur = sizes.length_of(&name).map(|len| {
                counts.insert(name.clone(), vec![(0, 0); len as usize / BIN + 1]);
                name
            });
            pos = 0;
            continue;
        }
        let Some(c) = &cur else { continue };
        let v = counts.get_mut(c).unwrap();
        for b in line.trim_end().bytes() {
            let i = pos / BIN;
            if i < v.len() {
                match b | 0x20 {
                    b'g' | b'c' => v[i].0 += 1,
                    b'a' | b't' => v[i].1 += 1,
                    _ => {}
                }
            }
            pos += 1;
        }
    }
    let mut bins: Vec<(String, i64, f64)> = Vec::new();
    let (mut s, mut ss, mut n) = (0.0f64, 0.0f64, 0.0f64);
    for (ci, chrom) in sizes.names.iter().enumerate() {
        let size = sizes.lengths[ci] as usize;
        let Some(v) = counts.get(chrom) else { continue };
        let mut i = 0usize;
        while i <= size {
            if i + BIN <= size {
                let (gc, at) = v[i / BIN];
                let (pgc, pat) = (gc as f64 / BIN as f64, at as f64 / BIN as f64);
                let tot = pat + pgc;
                if tot > 0.0 {
                    let x = awk(pgc / tot);
                    s += x;
                    ss += x * x;
                    n += 1.0;
                    bins.push((chrom.clone(), i as i64, x));
                }
            }
            i += BIN;
        }
    }
    if n == 0.0 {
        return Err(format!("{}: no callable bases; do the FASTA names match the sizes file?", path.display()));
    }
    let m = s / n;
    let thresh = awk(awk(m) - awk((ss / n - m * m).sqrt()));
    let mut out: Ivs = HashMap::new();
    for (c, st, x) in bins {
        if x < thresh {
            out.entry(c).or_default().push((st, st + BIN as i64));
        }
    }
    Ok((out, thresh))
}

pub struct V1Config {
    pub hic: PathBuf,
    pub sizes: ChromSizes,
    pub genes: Ivs,
    pub bbins: Ivs,
    pub res: u64,
    pub coarsest_res: u64,
    pub eigen_res: u64,
    pub eigen_init: bool,
    pub window_kb: f64,
    pub threshold: f64,
    pub reshift: bool,
    pub smoothing: bool,
    pub outdir: PathBuf,
    pub outpre: String,
    pub verbose: bool,
}

// ---------------------------------------------------------------------------
// Seeds: process_genes_Bbins
// ---------------------------------------------------------------------------

/// Bins a set of intervals marks at `res`: every bin from int(s/res) to
/// int(e/res) inclusive, so an interval ending exactly on a bin edge also
/// marks the next bin.
fn bin_keys(ivs: Option<&Vec<(i64, i64)>>, res: i64) -> HashSet<i64> {
    let mut k = HashSet::new();
    if let Some(v) = ivs {
        for &(s, e) in v {
            let a = (s / res) * res;
            let b = (e / res) * res;
            let mut i = a;
            while i <= b {
                k.insert(i);
                i += res;
            }
        }
    }
    k
}

// ---------------------------------------------------------------------------
// One chromosome at one resolution: dump, O/E, row statistics, RowZs
// ---------------------------------------------------------------------------

/// The RowZs table for one chromosome and resolution. Each off-diagonal
/// contact appears twice, as (row, col, z) with z normalised by `col`.
struct Prep {
    rows: Vec<i64>,
    cols: Vec<i64>,
    z: Vec<f64>,
}

fn prepare(hic_path: &Path, chrom: &str, size: i64, res: i64, threshold: f64) -> Result<Prep, String> {
    let mut f = HicFile::open(hic_path).map_err(|e| e.to_string())?;
    // bash's dumper sends hicstraw errors to a side file and carries on with an
    // empty dump, so a chromosome the map lacks simply gets no scores.
    let contacts = match f.straw_cis(ValueKind::Observed, "NONE", chrom, res as u32) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("      warning: {chrom} @ {res}: {e}; treated as empty, as bash does");
            Vec::new()
        }
    };

    // expecteds: per-diagonal sum over (size/res - d), awk-formatted
    let mut diag: BTreeMap<i64, f64> = BTreeMap::new();
    for c in &contacts {
        let d = (c.y as i64 - c.x as i64) / res;
        *diag.entry(d).or_insert(0.0) += c.value as f64;
    }
    let nb = size as f64 / res as f64;
    let expd: HashMap<i64, f64> = diag
        .into_iter()
        .map(|(d, s)| (d, awk(s / (nb - d as f64))))
        .collect();

    // symmonfile: off-diagonal O/E, clamped at the threshold, both orientations
    let n = contacts.len();
    let mut rows = Vec::with_capacity(2 * n);
    let mut cols = Vec::with_capacity(2 * n);
    let mut vals = Vec::with_capacity(2 * n);
    for c in &contacts {
        let (x, y) = (c.x as i64, c.y as i64);
        if x == y {
            continue;
        }
        let d = (y - x) / res;
        let mut v = awk((c.value as f64 + 1.0) / (expd[&d] + 1.0));
        if !(v < threshold) {
            v = threshold;
        }
        rows.push(x);
        cols.push(y);
        vals.push(v);
        rows.push(y);
        cols.push(x);
        vals.push(v);
    }

    // tmpmean: row mean and SD over rowbins = int(size/res) + 1, zeros implied
    let rowbins = (size / res + 1) as f64;
    let mut m: HashMap<i64, (f64, f64)> = HashMap::new();
    for (&r, &v) in rows.iter().zip(&vals) {
        let e = m.entry(r).or_insert((0.0, 0.0));
        e.0 += v;
        e.1 += v * v;
    }
    let stats: HashMap<i64, (f64, f64)> = m
        .into_iter()
        .map(|(k, (s, ss))| {
            let mean = s / rowbins;
            (k, (awk(mean), awk((ss / rowbins - mean * mean).sqrt())))
        })
        .collect();

    // RowZs: (v - m[col]) / sd[col], kept only where sd[col] > 0
    let mut out = Prep { rows: Vec::with_capacity(rows.len()), cols: Vec::with_capacity(rows.len()), z: Vec::with_capacity(rows.len()) };
    for i in 0..rows.len() {
        let Some(&(mc, dc)) = stats.get(&cols[i]) else { continue };
        if dc > 0.0 {
            out.rows.push(rows[i]);
            out.cols.push(cols[i]);
            out.z.push(awk((vals[i] - mc) / dc));
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// process_oppocheck_statement
// ---------------------------------------------------------------------------

/// One output line of a per-chromosome CRUSH file: (start, end, value).
type Line = (i64, i64, f64);

fn oppocheck(prep: &Prep, a: &HashSet<i64>, b: &HashSet<i64>, size: i64, res: i64, window: i64) -> Vec<Line> {
    // A2removedfile / B2removedfile: per partner bin, Σz and count over
    // contacts whose row is a seed, plus 0.1 once if the bin is itself a seed.
    let side = |seeds: &HashSet<i64>| -> HashMap<i64, (f64, u64)> {
        let mut acc: HashMap<i64, (f64, u64)> = HashMap::new();
        for i in 0..prep.z.len() {
            if seeds.contains(&prep.rows[i]) {
                let e = acc.entry(prep.cols[i]).or_insert((0.0, 0));
                e.0 += prep.z[i];
                e.1 += 1;
            }
        }
        let mut i = 0;
        while i <= size {
            if seeds.contains(&i) {
                let e = acc.entry(i).or_insert((0.0, 0));
                e.0 += 0.1;
                e.1 += 1;
            }
            i += res;
        }
        acc
    };
    let am = side(a);
    let bm = side(b);
    let ab_sum: f64 = am.values().chain(bm.values()).map(|&(_, c)| c as f64).sum();

    // A minus B, where a bin seen only by B keeps +B
    let mut sc: BTreeMap<i64, f64> = BTreeMap::new();
    for (&k, &(s, c)) in &am {
        sc.insert(k, awk(s * c as f64));
    }
    for (&k, &(s, c)) in &bm {
        let v = awk(s * c as f64);
        match sc.get_mut(&k) {
            Some(x) => *x -= v,
            None => {
                sc.insert(k, v);
            }
        }
    }

    // sliding window, numerically sorted by bin
    let mut raw: Vec<(i64, f64)> = Vec::new();
    let w = window.max(1);
    let mut ring = vec![0.0; w as usize];
    let (mut sum, mut count) = (0.0f64, 0.0f64);
    let mut nr: i64 = 0;
    for (&k, &x) in &sc {
        if x == 0.0 {
            continue;
        }
        let v = awk(x / ab_sum);
        nr += 1;
        let md = (nr % w) as usize;
        if nr <= w {
            raw.push((k, v));
            count += 1.0;
        } else {
            sum -= ring[md];
        }
        sum += v;
        ring[md] = v;
        let shift = ((res * w) as f64 / 2.0) as i64;
        raw.push((k - shift, awk(sum / count)));
    }

    let half = res / 2;
    let mut out: Vec<Line> = raw
        .into_iter()
        .map(|(p, v)| (p - half + res, p + half + res, v))
        .filter(|&(s, _, _)| s > 0)
        .collect();

    // endZ: subtract the chromosome mean, as printed by awk
    if !out.is_empty() {
        let mean = awk(out.iter().map(|l| l.2).sum::<f64>() / out.len() as f64);
        for l in out.iter_mut() {
            l.2 = awk(l.2 - mean);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Merging chromosomes; the shifter; seed re-initialisation
// ---------------------------------------------------------------------------

/// chrom -> (start, end) -> value; duplicate keys averaged, as the merge awk does.
type Track = HashMap<String, BTreeMap<(i64, i64), f64>>;

fn merge_chroms(per: Vec<(String, Vec<Line>)>) -> Track {
    let mut t: Track = HashMap::new();
    for (chrom, lines) in per {
        let mut acc: BTreeMap<(i64, i64), (f64, u32)> = BTreeMap::new();
        for (s, e, v) in lines {
            if v.is_nan() {
                continue;
            }
            let x = acc.entry((s, e)).or_insert((0.0, 0));
            x.0 += v;
            x.1 += 1;
        }
        let m = acc
            .into_iter()
            .map(|(k, (s, n))| (k, if n > 1 { awk(s / n as f64) } else { awk(s) }))
            .collect();
        t.insert(chrom, m);
    }
    t
}

/// numpy's pairwise_sum for float64 (umath/loops_utils.h), which np.mean
/// uses: sequential below 8 elements, eight interleaved accumulators up to
/// the 128-element block size, recursive halving above it.
fn np_pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut r = -0.0f64;
        for &x in a {
            r += x;
        }
        r
    } else if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - (n % 8) {
            for k in 0..8 {
                r[k] += a[i + k];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        np_pairwise_sum(&a[..n2]) + np_pairwise_sum(&a[n2..])
    }
}

/// run_midshifter: subtract a rolling mean over 3*(coarse/res) bins from every
/// bin of a zero-filled array. The coarse track is never read (B4).
fn midshifter(track: &Track, sizes: &ChromSizes, res: i64, coarse: i64) -> Track {
    let mut out: Track = HashMap::new();
    for (ci, chrom) in sizes.names.iter().enumerate() {
        let size = sizes.lengths[ci] as i64;
        let len = ((size + 2 * coarse) + res - 1) / res;
        let mut score = vec![0.0f64; len as usize];
        if let Some(rows) = track.get(chrom) {
            // file order: sorted by start then end; the last write to a bin wins
            for (&(s, _), &v) in rows {
                let p = s / res;
                if p >= 0 && p < len {
                    score[p as usize] = v;
                }
            }
        }
        let ws = 3 * (coarse / res);
        let mut m = BTreeMap::new();
        for i in 0..len {
            let lo = (i - ws / 2).max(0);
            let hi = (i + ws / 2).min(len);
            // np.mean of the slice: numpy's pairwise sum, so a window of
            // zeros averages to exactly 0 rather than a prefix-sum residue
            let mean = np_pairwise_sum(&score[lo as usize..hi as usize]) / (hi - lo) as f64;
            m.insert((i * res, i * res + res), score[i as usize] - mean);
        }
        out.insert(chrom.clone(), m);
    }
    out
}

/// separating_ABbins with switch=1: genes overlapping positive bins become the
/// A seeds and B-annotation bins overlapping negative bins the B seeds. On a
/// chromosome with no hit the original annotation is kept; if no chromosome
/// has a hit the awk "NR==FNR" idiom empties the set altogether.
fn separate(track: &Track, genes: &Ivs, bbins: &Ivs) -> (Ivs, Ivs) {
    let pick = |annot: &Ivs, positive: bool| -> Ivs {
        let mut hits: Ivs = HashMap::new();
        for (chrom, rows) in track {
            let Some(src) = annot.get(chrom) else { continue };
            let mut sel: Vec<(i64, i64)> =
                rows.iter().filter(|(_, &v)| if positive { v > 0.0 } else { v < 0.0 }).map(|(&k, _)| k).collect();
            if sel.is_empty() {
                continue;
            }
            sel.sort_unstable();
            let maxlen = sel.iter().map(|&(s, e)| e - s).max().unwrap_or(0);
            let starts: Vec<i64> = sel.iter().map(|&(s, _)| s).collect();
            let mut keep = Vec::new();
            for &(s, e) in src {
                // any selected interval with start in (s - maxlen, e) that overlaps [s, e)
                let lo = starts.partition_point(|&x| x <= s - maxlen);
                let hi = starts.partition_point(|&x| x < e);
                if sel[lo..hi].iter().any(|&(rs, re)| rs < e && s < re) {
                    keep.push((s, e));
                }
            }
            if !keep.is_empty() {
                hits.insert(chrom.clone(), keep);
            }
        }
        if hits.is_empty() {
            return HashMap::new();
        }
        let mut combined = hits.clone();
        for (chrom, v) in annot {
            if !hits.contains_key(chrom) {
                combined.insert(chrom.clone(), v.clone());
            }
        }
        combined
    };
    (pick(genes, true), pick(bbins, false))
}

// ---------------------------------------------------------------------------
// Eigenvector initialisation: run_EigenVector, check_and_flip_sign, find_best_pc
// ---------------------------------------------------------------------------

fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma) * (x - ma);
        sbb += (y - mb) * (y - mb);
    }
    let r = sab / (saa * sbb).sqrt();
    r.is_finite().then_some(r)
}

/// Returns the A (positive) and B (negative) 100 kb intervals of the chosen PC.
fn eigen_seeds(cfg: &V1Config, chrom: &str, size: i64) -> Result<Option<(Vec<(i64, i64)>, Vec<(i64, i64)>, String)>, String> {
    let res = cfg.eigen_res as i64;
    let mut f = HicFile::open(&cfg.hic).map_err(|e| e.to_string())?;
    let Ok(recs) = f.straw_cis(ValueKind::Oe, "VC_SQRT", chrom, res as u32) else {
        return Ok(None); // as above: no PCs, so no eigenvector seed for this chromosome
    };
    let recs: Vec<_> = recs.into_iter().filter(|c| c.value.is_finite()).collect();
    if recs.is_empty() {
        return Ok(None);
    }
    let n = recs.iter().map(|c| c.x.max(c.y) as i64 / res).max().unwrap() as usize + 1;
    let mut mat = vec![0.0f64; n * n];
    for c in &recs {
        let (i, j) = (c.x as usize / res as usize, c.y as usize / res as usize);
        mat[i * n + j] = c.value as f64;
        mat[j * n + i] = c.value as f64;
    }
    let keep: Vec<usize> = (0..n).filter(|&j| (0..n).map(|i| mat[i * n + j]).sum::<f64>() > 0.0).collect();
    let k = keep.len();
    if k < 3 {
        return Ok(None);
    }
    // np.corrcoef over rows
    let mut z = vec![0.0f64; k * k];
    for (a, &i) in keep.iter().enumerate() {
        let row: Vec<f64> = keep.iter().map(|&j| mat[i * n + j]).collect();
        let mu = row.iter().sum::<f64>() / k as f64;
        let sd = (row.iter().map(|x| (x - mu) * (x - mu)).sum::<f64>()).sqrt();
        for b in 0..k {
            z[a * k + b] = (row[b] - mu) / sd;
        }
    }
    let mut corr = DMatrix::<f64>::zeros(k, k);
    for a in 0..k {
        for b in a..k {
            let r: f64 = (0..k).map(|t| z[a * k + t] * z[b * k + t]).sum();
            corr[(a, b)] = r;
            corr[(b, a)] = r;
        }
    }
    if corr.iter().any(|x| !x.is_finite()) {
        return Ok(None);
    }
    let eig = corr.symmetric_eigen();
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&x, &y| eig.eigenvalues[y].partial_cmp(&eig.eigenvalues[x]).unwrap());

    let genes = cfg.genes.get(chrom).cloned().unwrap_or_default();
    let num_bins = (size / res) as usize;
    let mut gene_d = vec![0.0f64; num_bins];
    for &(s, e) in &genes {
        if s < 0 || e < 0 {
            continue;
        }
        for i in (s / res)..=(e / res) {
            if (i as usize) < num_bins {
                gene_d[i as usize] += 1.0;
            }
        }
    }
    // gcbinned: B annotation snapped to eigenres, end capped at the chromosome
    let mut gc_d = vec![0.0f64; num_bins];
    for &(s, e) in cfg.bbins.get(chrom).map(|v| v.as_slice()).unwrap_or(&[]) {
        if s < 0 || e < 0 {
            continue;
        }
        let st = (s / res) * res;
        let en = ((e / res) * res).min(size);
        if st >= size {
            continue;
        }
        for i in (st / res)..=(en / res) {
            if (i as usize) < num_bins {
                gc_d[i as usize] += 1.0;
            }
        }
    }
    let overlaps_gene = |s: i64, e: i64| genes.iter().any(|&(gs, ge)| gs < e && s < ge);
    let gene_bin: Vec<bool> = keep.iter().map(|&i| overlaps_gene(i as i64 * res, i as i64 * res + res)).collect();

    // PC files in glob order (PC1, PC10, PC2, ..., PC9); strict ">" keeps the first maximum
    let npc = 10.min(k);
    let mut glob: Vec<usize> = (1..=npc).collect();
    glob.sort_by_key(|p| p.to_string());
    let mut best: Option<(f64, usize, Vec<f64>)> = None;
    for &pcn in &glob {
        let col = order[pcn - 1];
        let mut v: Vec<f64> = (0..k).map(|a| eig.eigenvectors[(a, col)]).collect();
        // check_and_flip_sign: count gene-overlapping bins of each sign
        let (mut npos, mut nneg) = (0u64, 0u64);
        for a in 0..k {
            if gene_bin[a] {
                if v[a] > 0.0 {
                    npos += 1;
                } else if v[a] < 0.0 {
                    nneg += 1;
                }
            }
        }
        // the awk compares "neg,pos"; with only one class present the empty
        // second field makes the comparison true and the PC is flipped
        let flip = match (nneg > 0, npos > 0) {
            (true, true) => nneg > npos,
            (false, false) => false,
            _ => true,
        };
        if flip {
            v.iter_mut().for_each(|x| *x = -*x);
        }
        let mut pc = vec![0.0f64; num_bins];
        for (a, &i) in keep.iter().enumerate() {
            if i < num_bins {
                pc[i] = awk(v[a]);
            }
        }
        let (Some(gr), Some(br)) = (pearson(&pc, &gene_d), pearson(&pc, &gc_d)) else { continue };
        let total = gr - br;
        if best.as_ref().map_or(total > -100.0, |b| total > b.0) {
            best = Some((total, pcn, v));
        }
    }
    let Some((total, pcn, v)) = best else { return Ok(None) };
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for (idx, &i) in keep.iter().enumerate() {
        let s = i as i64 * res;
        if v[idx] > 0.0 {
            a.push((s, s + res));
        } else if v[idx] < 0.0 {
            b.push((s, s + res));
        }
    }
    Ok(Some((a, b, format!("PC{pcn} (gene - gc = {total:.4})"))))
}

// ---------------------------------------------------------------------------
// The walk
// ---------------------------------------------------------------------------

/// bedtools groupBy's printed mean: 10 significant digits (measured on
/// bedtools v2.29.0: a mean of 1/3 prints as 0.3333333333).
fn prec10(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    format!("{:.9e}", x).parse().unwrap_or(x)
}

/// One bedGraph as bash holds it between steps: rows in file order.
type Rows = Vec<(String, i64, i64, f64)>;

fn rows_sorted(t: &Track, sizes: &ChromSizes) -> Rows {
    let mut out = Vec::new();
    for c in &sizes.names {
        if let Some(m) = t.get(c) {
            for (&(s, e), &v) in m {
                out.push((c.clone(), s, e, v));
            }
        }
    }
    out
}

/// run_endshifter: re-centre `high` on `coarse` through overlapping windows of
/// 2 coarse bins left and 3 right, the three-step intersect/groupBy chain.
fn endshifter(high: &Rows, coarse: &Rows, sizes: &ChromSizes, coarser: i64) -> Rows {
    let by_chrom = |r: &Rows| -> HashMap<String, Vec<(i64, i64, f64)>> {
        let mut m: HashMap<String, Vec<(i64, i64, f64)>> = HashMap::new();
        for (c, s, e, v) in r {
            m.entry(c.clone()).or_default().push((*s, *e, *v));
        }
        m
    };
    // all rows overlapping [s, e) of a list sorted by start, within maxlen
    fn hits(v: &[(i64, i64, f64)], starts: &[i64], maxlen: i64, s: i64, e: i64) -> Vec<usize> {
        let lo = starts.partition_point(|&x| x <= s - maxlen);
        let hi = starts.partition_point(|&x| x < e);
        (lo..hi).filter(|&k| v[k].0 < e && s < v[k].1).collect()
    }
    let hb = by_chrom(high);
    let cb = by_chrom(coarse);
    let mut out = Vec::new();
    for (ci, chrom) in sizes.names.iter().enumerate() {
        let size = sizes.lengths[ci] as i64;
        let (Some(h), Some(c)) = (hb.get(chrom), cb.get(chrom)) else { continue };
        let mut hs = h.clone();
        hs.sort_by_key(|x| (x.0, x.1));
        let hst: Vec<i64> = hs.iter().map(|x| x.0).collect();
        let hmax = hs.iter().map(|x| x.1 - x.0).max().unwrap_or(0);
        // step 1: windows and the mean high value inside each
        let mut wins: Vec<(i64, i64, f64)> = Vec::new();
        let mut i = 0;
        while i <= size {
            let (s, e) = ((i - 2 * coarser).max(0), (i + 3 * coarser).min(size));
            let k = hits(&hs, &hst, hmax, s, e);
            if !k.is_empty() {
                wins.push((s, e, prec10(k.iter().map(|&j| hs[j].2).sum::<f64>() / k.len() as f64)));
            }
            i += coarser;
        }
        wins.sort_by_key(|x| (x.0, x.1));
        let wst: Vec<i64> = wins.iter().map(|x| x.0).collect();
        let wmax = wins.iter().map(|x| x.1 - x.0).max().unwrap_or(0);
        // step 2: per coarse row, (window mean - coarse value) for each window it meets
        let mut diffs: Vec<(i64, i64, f64)> = Vec::new();
        for &(s, e, v) in c {
            for j in hits(&wins, &wst, wmax, s, e) {
                diffs.push((s, e, awk(wins[j].2 - v)));
            }
        }
        diffs.sort_by_key(|x| (x.0, x.1));
        let dst: Vec<i64> = diffs.iter().map(|x| x.0).collect();
        let dmax = diffs.iter().map(|x| x.1 - x.0).max().unwrap_or(0);
        // step 3: per high row, subtract the mean shift of the rows it meets
        for &(s, e, v) in h {
            let k = hits(&diffs, &dst, dmax, s, e);
            if !k.is_empty() {
                let m = prec10(k.iter().map(|&j| diffs[j].2).sum::<f64>() / k.len() as f64);
                out.push((chrom.clone(), s, e, awk(v - m)));
            }
        }
    }
    out
}

/// run_smoothing: 3-point mean over the file sorted by chromosome and start,
/// pasted against itself shifted up and down one line. It runs across
/// chromosome boundaries, and the last line averages against an empty field.
fn smoothing(rows: &Rows, sizes: &ChromSizes) -> Rows {
    let order: HashMap<&str, usize> = sizes.names.iter().enumerate().map(|(i, c)| (c.as_str(), i)).collect();
    let mut r = rows.clone();
    r.sort_by_key(|x| (order.get(x.0.as_str()).copied().unwrap_or(usize::MAX), x.1));
    let n = r.len();
    (0..n)
        .map(|i| {
            let right = if i + 1 < n { r[i + 1].3 } else { 0.0 };
            let left = if i == 0 { r[0].3 } else { r[i - 1].3 };
            (r[i].0.clone(), r[i].1, r[i].2, awk((r[i].3 + right + left) / 3.0))
        })
        .collect()
}

/// Debug: with CRUSH_V1_DUMP=<dir>, write each intermediate track under the
/// name bash gives it in its temp directory, for step-by-step comparison.
fn dump(name: &str, t: &Track, sizes: &ChromSizes) {
    let Ok(dir) = std::env::var("CRUSH_V1_DUMP") else { return };
    let _ = std::fs::create_dir_all(&dir);
    let rows = rows_sorted(t, sizes);
    let _ = std::fs::write(
        Path::new(&dir).join(format!("{name}.bedgraph")),
        rows.iter().map(|(c, s, e, v)| format!("{c}\t{s}\t{e}\t{}\n", fmt_g6(*v))).collect::<String>(),
    );
}

fn window_for(window_kb: f64, res: i64) -> i64 {
    if window_kb == 0.0 {
        1
    } else {
        (window_kb / (res as f64 / 1000.0)) as i64 + 1
    }
}

/// The final merge: each rung's line is spread over `myres` bins after
/// dividing its value by its width, then averaged per bin.
fn final_merge(rungs: &[Track], myres: i64) -> Track {
    let mut acc: HashMap<String, BTreeMap<(i64, i64), (f64, u32)>> = HashMap::new();
    for t in rungs {
        for (chrom, rows) in t {
            let m = acc.entry(chrom.clone()).or_default();
            for (&(s, e), &v) in rows {
                let a = awk(s as f64 / myres as f64);
                let b = awk(e as f64 / myres as f64);
                let val = awk(v / (e - s) as f64);
                let mut i = a;
                while i < b {
                    let ii = i.trunc() as i64;
                    let x = m.entry((ii * myres, (ii + 1) * myres)).or_insert((0.0, 0));
                    x.0 += val;
                    x.1 += 1;
                    i += 1.0;
                }
            }
        }
    }
    acc.into_iter()
        .map(|(c, m)| (c, m.into_iter().map(|(k, (s, n))| (k, awk(s / n as f64))).collect()))
        .collect()
}

/// GIave over every bin, the ×100/GIave rescale, then drop bins that run past
/// the chromosome end: the mergedCrush_<res>.bedgraph bash writes.
fn finalize(track: &Track, sizes: &ChromSizes) -> Rows {
    let n: usize = track.values().map(|m| m.len()).sum();
    let gi = awk(track.values().flat_map(|m| m.values()).map(|v| v.abs()).sum::<f64>() / n as f64);
    let mut out = Vec::new();
    for (ci, chrom) in sizes.names.iter().enumerate() {
        let size = sizes.lengths[ci] as i64;
        let Some(m) = track.get(chrom) else { continue };
        for (&(s, e), &v) in m {
            if e > 1 && s < size && e <= size {
                out.push((chrom.clone(), s, e, awk(v / (gi / 100.0))));
            }
        }
    }
    out
}

fn write_rows(path: &Path, rows: &Rows) -> Result<(), String> {
    let f = File::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut w = BufWriter::new(f);
    writeln!(w, "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-150:150 autoScale off")
        .map_err(|e| e.to_string())?;
    for (c, s, e, v) in rows {
        writeln!(w, "{c}\t{s}\t{e}\t{}", fmt_g6(*v)).map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// awk's printed form of a number.
fn fmt_g6(x: f64) -> String {
    if x == x.trunc() && x.abs() < 9.0e15 {
        return format!("{}", x as i64);
    }
    let s = format!("{:.5e}", x);
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    if (-5..6).contains(&exp) {
        let decimals = (5 - exp).max(0) as usize;
        let t = format!("{:.*}", decimals, x);
        let t = if t.contains('.') { t.trim_end_matches('0').trim_end_matches('.').to_string() } else { t };
        t
    } else {
        let m = mant.trim_end_matches('0').trim_end_matches('.');
        format!("{m}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    }
}

pub fn run_v1(cfg: &V1Config, cpu: usize) -> Result<Vec<PathBuf>, String> {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(cpu.max(1)).build().map_err(|e| e.to_string())?;
    let hic = HicFile::open(&cfg.hic).map_err(|e| e.to_string())?;
    // listres.py: the file's resolutions in file order, within [res, coarsest]
    let ladder: Vec<i64> = hic
        .bp_resolutions
        .iter()
        .map(|&r| r as i64)
        .filter(|&r| r >= cfg.res as i64 && r <= cfg.coarsest_res as i64)
        .collect();
    drop(hic);
    if ladder.is_empty() {
        return Err("no resolution in the .hic falls inside [-r, -m]".into());
    }
    let maxres = *ladder.iter().max().unwrap();
    println!("  v1 ladder  : {ladder:?}");

    let chroms: Vec<(String, i64)> =
        cfg.sizes.names.iter().cloned().zip(cfg.sizes.lengths.iter().map(|&l| l as i64)).collect();

    // round-0 seeds
    let (eva, evb): (Ivs, Ivs) = if cfg.eigen_init {
        let res: Vec<_> = pool.install(|| {
            chroms.par_iter().map(|(c, s)| (c.clone(), eigen_seeds(cfg, c, *s))).collect()
        });
        let (mut a, mut b) = (HashMap::new(), HashMap::new());
        for (c, r) in res {
            match r? {
                Some((x, y, why)) => {
                    if cfg.verbose {
                        println!("      {c}: {why}");
                    }
                    a.insert(c.clone(), x);
                    b.insert(c, y);
                }
                None => println!("      {c}: no eigenvector seed"),
            }
        }
        (a, b)
    } else {
        (cfg.genes.clone(), cfg.bbins.clone())
    };

    let score_all = |res: i64, window: i64, a: &Ivs, b: &Ivs| -> Result<Track, String> {
        let per: Vec<Result<(String, Vec<Line>), String>> = pool.install(|| {
            chroms
                .par_iter()
                .map(|(c, size)| {
                    let prep = prepare(&cfg.hic, c, *size, res, cfg.threshold)?;
                    let ak = bin_keys(a.get(c), res);
                    let bk = bin_keys(b.get(c), res);
                    Ok((c.clone(), oppocheck(&prep, &ak, &bk, *size, res, window)))
                })
                .collect()
        });
        Ok(merge_chroms(per.into_iter().collect::<Result<_, _>>()?))
    };

    let mut seeds: (Ivs, Ivs) = (HashMap::new(), HashMap::new());
    let (mut coarse_res, mut crtmp) = (maxres, maxres);
    let mut written = Vec::new();
    let mut finals: Vec<(i64, Rows)> = Vec::new();
    std::fs::create_dir_all(&cfg.outdir).map_err(|e| e.to_string())?;

    for i in 0..ladder.len() {
        let myres = ladder[i];
        let window = window_for(cfg.window_kb, myres);
        println!("  v1 rung {myres} bp (window {window})");

        // process_resolution
        let (a, b) = if i == 0 { (&eva, &evb) } else { (&seeds.0, &seeds.1) };
        let full = score_all(myres, window, a, b)?;
        dump(&format!("Crush_{myres}"), &full, &cfg.sizes);
        if i > 0 {
            if coarse_res % myres == 0 {
                crtmp = coarse_res;
            } else if crtmp % myres == 0 {
                coarse_res = crtmp;
            }
            let sh = midshifter(&full, &cfg.sizes, myres, coarse_res);
            seeds = separate(&sh, &cfg.genes, &cfg.bbins);
            coarse_res = myres;
        } else {
            seeds = separate(&full, &cfg.genes, &cfg.bbins);
        }

        // reprocess_resolutions_with_shifter
        if i > 0 {
            let (mut crr, mut crtmp_r) = (maxres, maxres);
            let mut rungs = Vec::with_capacity(i + 1);
            for j in 0..=i {
                let prev = ladder[j];
                let re = score_all(prev, window, &seeds.0, &seeds.1)?;
                dump(&format!("Crush_Re_{prev}"), &re, &cfg.sizes);
                if j > 0 {
                    if crr % prev == 0 {
                        crtmp_r = crr;
                    } else if crtmp_r % prev == 0 {
                        crr = crtmp_r;
                    }
                    let sh = midshifter(&re, &cfg.sizes, prev, crr);
                    seeds = separate(&sh, &cfg.genes, &cfg.bbins);
                    crr = prev;
                } else {
                    seeds = separate(&re, &cfg.genes, &cfg.bbins);
                }
                rungs.push(re);
            }
            finals.push((myres, finalize(&final_merge(&rungs, myres), &cfg.sizes)));
        }
    }

    // Final Shifter (--reshift 1): each rung from the third on is re-centred
    // on the previous, already re-centred, rung; --smoothing 1 smooths each.
    if cfg.reshift {
        for k in 1..finals.len() {
            let coarser = finals[k - 1].0;
            let mut shifted = endshifter(&finals[k].1, &finals[k - 1].1, &cfg.sizes, coarser);
            if cfg.smoothing {
                shifted = smoothing(&shifted, &cfg.sizes);
            }
            finals[k].1 = shifted;
        }
    }
    for (res, rows) in &finals {
        let path = cfg.outdir.join(format!("{}mergedCrush_{res}.bedgraph", cfg.outpre));
        write_rows(&path, rows)?;
        written.push(path);
    }
    Ok(written)
}
