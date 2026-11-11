//! Post-hoc check: does the contact map actually support the track CRUSH just
//! wrote?
//!
//! CRUSH rescales its output to a fixed mean magnitude, so it emits a signed,
//! compartment-shaped track for any input whatsoever. "A track appeared" is
//! therefore not evidence that compartments exist. This module measures whether
//! the partition the track proposes is respected by the contacts:
//!
//!     strength = (AA + BB) / 2AB
//!
//! over the top and bottom quartiles of the track, on observed/expected values
//! with a per-separation expected. 1.0 means the partition carries no contact
//! preference at all.
//!
//! FOUR THINGS MAKE THE NUMBER MEAN SOMETHING.
//!
//! 1. A FIXED SCORING SCALE. Binning a track before quantiling smooths it, and
//!    smoother partitions score higher: the same call scores 1.55 scored at
//!    5 kb and 1.83 at 500 kb. Scoring always at FLAG_SCALE, whatever -r the run
//!    used, is what makes one threshold portable across resolutions.
//!
//! 2. THE ZEROS ARE COUNTED. The expected at each separation divides by the
//!    number of POSSIBLE mappable pairs, not the number of pairs present in the
//!    file, and each category mean likewise. Averaging over stored records only
//!    silently drops empty pairs -- and cross-compartment pairs are the ones
//!    most often empty, so dropping them inflates AB and destroys the contrast.
//!
//! 3. DETRENDING. A track that drifts monotonically along the chromosome is
//!    "respected by the map" for a trivial reason: one end of a chromosome
//!    contacts the other end less. Mitotic tracks do exactly this. Subtracting
//!    a wide moving average leaves only the alternating component, which is what
//!    "compartment" means.
//!
//! 4. A NULL FROM THE TRACK ITSELF. The track is circularly shifted along the
//!    chromosome: same length, same values, same autocorrelation -- just as
//!    compartment-shaped as the real thing -- but no longer registered to the
//!    map. That is the honest baseline, and it lands at 1.0.
//!
//! Orientation is irrelevant throughout: flipping a track swaps AA with BB.
//!
//! CALIBRATION. Thresholds come from maps with a known answer: the inactive X
//! (no compartments) against the active X and an Xist-KO inactive X, across
//! CRUSH calls spanning -r 5 kb to 500 kb; and a mitotic time course
//! (prometaphase flat, G2 and non-synchronised compartmentalised) across three
//! chromosomes. Everything without compartments scored <= 1.12 and everything
//! with them >= 1.54. A pachytene map -- whose compartments are unusually
//! short-range -- sits at 1.27, between the two, and is reported as
//! inconclusive rather than forced to a verdict.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::bed::ChromSizes;
use crate::hic::{HicFile, ValueKind};

/// Normalisations tried, in order, for the check itself.
///
/// Deliberately NOT the run's -N. On an unbalanced matrix this statistic
/// collapses onto coverage rather than compartments -- an unnormalised mitotic
/// series scored interphase BELOW prometaphase, i.e. backwards -- so the check
/// insists on a balanced matrix whatever the calling used, and falls back to the
/// run's own choice only if the file stores nothing better.
const FLAG_NORMS: [&str; 4] = ["KR", "SCALE", "VC_SQRT", "VC"];

/// Scoring scale. Fixed on purpose -- see note 1 above.
pub const FLAG_SCALE: u64 = 100_000;
/// Contacts closer than this are excluded: at short separation everything
/// contacts everything, which dilutes the contrast without being informative.
const MIN_SEP: u64 = 1_000_000;
/// Moving-average width for detrending (note 3).
const DETREND_WIN: u64 = 10_000_000;
/// Circular shifts used to build the null.
const N_NULL: usize = 16;
/// Below this, warn: the map does not support the track.
pub const WARN_BELOW: f64 = 1.20;
/// At or above this, the track is supported.
pub const PASS_AT: f64 = 1.40;
/// Chromosomes with fewer usable bins than this are skipped as uninformative.
const MIN_BINS: usize = 50;

pub struct ChromFlag {
    pub chrom: String,
    pub strength: f64,
    pub null: f64,
    pub z: f64,
    pub bins: usize,
}

pub struct FlagReport {
    pub scale: u64,
    pub per_chrom: Vec<ChromFlag>,
    pub median: f64,
}

impl FlagReport {
    pub fn verdict(&self) -> &'static str {
        if self.median < WARN_BELOW {
            "unsupported"
        } else if self.median < PASS_AT {
            "inconclusive"
        } else {
            "supported"
        }
    }
}

/// Bin a track onto an `n`-bin grid at `scale`, averaging within each bin.
fn bin_track(intervals: &[(u64, u64, f64)], scale: u64, n: usize) -> Vec<f64> {
    let mut sum = vec![0.0f64; n];
    let mut cnt = vec![0u32; n];
    for &(s, _, v) in intervals {
        if !v.is_finite() {
            continue;
        }
        let i = (s / scale) as usize;
        if i < n {
            sum[i] += v;
            cnt[i] += 1;
        }
    }
    (0..n)
        .map(|i| {
            if cnt[i] > 0 {
                sum[i] / cnt[i] as f64
            } else {
                f64::NAN
            }
        })
        .collect()
}

/// Subtract a wide moving average, NaN-aware (note 3).
fn detrend(v: &[f64], scale: u64) -> Vec<f64> {
    let n = v.len();
    let w = ((DETREND_WIN / scale) as usize).max(3);
    let half = w / 2;
    // Prefix sums over finite entries so the window is O(1) per bin.
    let mut psum = vec![0.0f64; n + 1];
    let mut pcnt = vec![0u32; n + 1];
    for i in 0..n {
        let ok = v[i].is_finite();
        psum[i + 1] = psum[i] + if ok { v[i] } else { 0.0 };
        pcnt[i + 1] = pcnt[i] + u32::from(ok);
    }
    (0..n)
        .map(|i| {
            if !v[i].is_finite() {
                return f64::NAN;
            }
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(n);
            let c = pcnt[hi] - pcnt[lo];
            if c == 0 {
                f64::NAN
            } else {
                v[i] - (psum[hi] - psum[lo]) / c as f64
            }
        })
        .collect()
}

/// Top and bottom quartile masks over bins that are both called and mappable.
fn quartile_masks(v: &[f64], valid: &[bool]) -> Option<(Vec<bool>, Vec<bool>)> {
    let mut vals: Vec<f64> = v
        .iter()
        .zip(valid)
        .filter(|(x, ok)| x.is_finite() && **ok)
        .map(|(x, _)| *x)
        .collect();
    if vals.len() < MIN_BINS {
        return None;
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let lo = vals[vals.len() / 4];
    let hi = vals[vals.len() * 3 / 4];
    let top = v
        .iter()
        .zip(valid)
        .map(|(x, ok)| *ok && x.is_finite() && *x >= hi)
        .collect();
    let bot = v
        .iter()
        .zip(valid)
        .map(|(x, ok)| *ok && x.is_finite() && *x <= lo)
        .collect();
    Some((top, bot))
}

/// Per-separation pair counts: mappable pairs, and pairs within each category.
/// Exact, and the reason the zeros are counted correctly (note 2).
fn lag_counts(
    top: &[bool],
    bot: &[bool],
    valid: &[bool],
    n: usize,
    min_bin: usize,
) -> (Vec<u64>, Vec<u64>, Vec<u64>, Vec<u64>) {
    let (mut poss, mut aa, mut bb, mut ab) =
        (vec![0u64; n], vec![0u64; n], vec![0u64; n], vec![0u64; n]);
    for i in 0..n {
        if !valid[i] {
            continue;
        }
        let (ti, bi) = (top[i], bot[i]);
        for j in (i + min_bin)..n {
            if !valid[j] {
                continue;
            }
            let d = j - i;
            poss[d] += 1;
            if ti && top[j] {
                aa[d] += 1;
            } else if bi && bot[j] {
                bb[d] += 1;
            } else if (ti && bot[j]) || (bi && top[j]) {
                ab[d] += 1;
            }
        }
    }
    (poss, aa, bb, ab)
}

/// (AA+BB)/2AB for one partition against pre-binned contacts.
#[allow(clippy::too_many_arguments)]
fn strength(
    xs: &[u32],
    ys: &[u32],
    vs: &[f32],
    n: usize,
    top: &[bool],
    bot: &[bool],
    valid: &[bool],
    min_bin: usize,
    poss: &[u64],
    obs_tot: &[f64],
) -> Option<f64> {
    let (_, caa, cbb, cab) = lag_counts(top, bot, valid, n, min_bin);
    let mut oaa = vec![0.0f64; n];
    let mut obb = vec![0.0f64; n];
    let mut oab = vec![0.0f64; n];
    for k in 0..xs.len() {
        let (i, j) = (xs[k] as usize, ys[k] as usize);
        if i >= n || j >= n || !valid[i] || !valid[j] {
            continue;
        }
        let d = j.abs_diff(i);
        if d < min_bin {
            continue;
        }
        let val = vs[k] as f64;
        if top[i] && top[j] {
            oaa[d] += val;
        } else if bot[i] && bot[j] {
            obb[d] += val;
        } else if (top[i] && bot[j]) || (bot[i] && top[j]) {
            oab[d] += val;
        }
    }
    let mean = |obs: &[f64], cnt: &[u64]| -> Option<f64> {
        let (mut num, mut den) = (0.0f64, 0.0f64);
        for d in min_bin..n {
            if poss[d] == 0 || cnt[d] == 0 {
                continue;
            }
            let exp = obs_tot[d] / poss[d] as f64;
            if exp <= 0.0 {
                continue;
            }
            num += obs[d] / exp;
            den += cnt[d] as f64;
        }
        (den > 0.0).then_some(num / den)
    };
    let (a, b, x) = (mean(&oaa, &caa)?, mean(&obb, &cbb)?, mean(&oab, &cab)?);
    (x > 0.0).then(|| (a + b) / (2.0 * x))
}

/// Score every chromosome and reduce to a median.
pub fn evaluate(
    hic: &std::path::Path,
    norm: &str,
    sizes: &ChromSizes,
    track: &HashMap<String, Vec<(u64, u64, f64)>>,
    available: &[u32],
) -> Option<FlagReport> {
    // Fixed scale where the file has it; otherwise the nearest stored
    // resolution at or above it, so the number stays comparable.
    let scale = available
        .iter()
        .map(|&r| r as u64)
        .filter(|&r| r >= FLAG_SCALE)
        .min()
        .or_else(|| available.iter().map(|&r| r as u64).max())?;
    let min_bin = ((MIN_SEP / scale) as usize).max(1);

    let mut per_chrom: Vec<ChromFlag> = sizes
        .names
        .par_iter()
        .filter_map(|chrom| {
            let intervals = track.get(chrom)?;
            let mut f = HicFile::open(hic).ok()?;
            let contacts = FLAG_NORMS
                .iter()
                .copied()
                .chain(std::iter::once(norm))
                .find_map(|nm| {
                    f.straw_cis(ValueKind::Observed, nm, chrom, scale as u32)
                        .ok()
                        .filter(|c| !c.is_empty())
                })?;
            let len = f.chromosome(chrom)?.length as u64;
            let n = (len / scale) as usize + 1;
            if n < MIN_BINS {
                return None;
            }

            let mut xs = Vec::with_capacity(contacts.len());
            let mut ys = Vec::with_capacity(contacts.len());
            let mut vs = Vec::with_capacity(contacts.len());
            let mut valid = vec![false; n];
            for c in &contacts {
                let (i, j) = ((c.x as u64 / scale) as usize, (c.y as u64 / scale) as usize);
                if i >= n || j >= n || !c.value.is_finite() || c.value <= 0.0 {
                    continue;
                }
                xs.push(i as u32);
                ys.push(j as u32);
                vs.push(c.value);
                valid[i] = true;
                valid[j] = true;
            }
            if xs.is_empty() {
                return None;
            }

            let binned = bin_track(intervals, scale, n);
            let det = detrend(&binned, scale);
            let (top, bot) = quartile_masks(&det, &valid)?;

            // Shared across the real partition and every null.
            let (poss, _, _, _) = lag_counts(&top, &bot, &valid, n, min_bin);
            let mut obs_tot = vec![0.0f64; n];
            for k in 0..xs.len() {
                let d = (ys[k] as usize).abs_diff(xs[k] as usize);
                if d >= min_bin {
                    obs_tot[d] += vs[k] as f64;
                }
            }

            let s = strength(
                &xs, &ys, &vs, n, &top, &bot, &valid, min_bin, &poss, &obs_tot,
            )?;

            // Null: circular shifts of the detrended track (note 4). Offsets are
            // spread deterministically over the middle of the chromosome so the
            // report is reproducible run to run.
            let mut nulls = Vec::with_capacity(N_NULL);
            for k in 0..N_NULL {
                let off = n / 10 + (k + 1) * (n * 8 / 10) / (N_NULL + 1);
                let rolled: Vec<f64> = (0..n).map(|i| det[(i + n - off % n) % n]).collect();
                if let Some((t2, b2)) = quartile_masks(&rolled, &valid) {
                    if let Some(q) = strength(
                        &xs, &ys, &vs, n, &t2, &b2, &valid, min_bin, &poss, &obs_tot,
                    ) {
                        nulls.push(q);
                    }
                }
            }
            let (mu, sd) = if nulls.is_empty() {
                (1.0, f64::NAN)
            } else {
                let m = nulls.iter().sum::<f64>() / nulls.len() as f64;
                let var = nulls.iter().map(|q| (q - m) * (q - m)).sum::<f64>()
                    / nulls.len() as f64;
                (m, var.sqrt())
            };

            Some(ChromFlag {
                chrom: chrom.clone(),
                strength: s,
                null: mu,
                z: if sd > 0.0 { (s - mu) / sd } else { f64::NAN },
                bins: valid.iter().filter(|b| **b).count(),
            })
        })
        .collect();

    if per_chrom.is_empty() {
        return None;
    }
    per_chrom.sort_by(|a, b| a.chrom.cmp(&b.chrom));
    let mut ss: Vec<f64> = per_chrom.iter().map(|c| c.strength).collect();
    ss.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if ss.len() % 2 == 1 {
        ss[ss.len() / 2]
    } else {
        0.5 * (ss[ss.len() / 2 - 1] + ss[ss.len() / 2])
    };
    Some(FlagReport {
        scale,
        per_chrom,
        median,
    })
}
