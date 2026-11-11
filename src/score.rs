//! The CRUSH score itself: distance normalisation, row Z-scoring, and the
//! A-versus-B interaction preference for every bin of one chromosome.
//!
//! This mirrors the bash pipeline's `task()` + `process_oppocheck_statement()`
//! stage by stage, but keeps everything in memory. The bash version writes the
//! fully symmetrised contact list to disk (`symmonfile`) and reads it back
//! twice; here the row statistics are accumulated in a single pass and the
//! symmetric partner is generated on the fly, which is where most of the
//! speed-up comes from.

use std::collections::HashMap;

use crate::hic::Contact;

/// Per-bin CRUSH result before the sliding-window average.
#[derive(Debug, Clone)]
pub struct BinScore {
    pub bin_start: u64,
    pub score: f64,
}

/// Welch-test inputs for one bin, retained so p-values can be computed later.
#[derive(Debug, Clone, Copy, Default)]
pub struct Moments {
    pub sum: f64,
    pub count: f64,
    pub sumsq: f64,
}

impl Moments {
    fn push(&mut self, v: f64) {
        self.sum += v;
        self.count += 1.0;
        self.sumsq += v * v;
    }
    pub fn mean(&self) -> f64 {
        if self.count == 0.0 {
            0.0
        } else {
            self.sum / self.count
        }
    }
    /// Population SD, matching the bash pipeline's `sqrt(E[x^2] - E[x]^2)`.
    pub fn sd(&self) -> f64 {
        if self.count == 0.0 {
            return 0.0;
        }
        let m = self.mean();
        let v = self.sumsq / self.count - m * m;
        if v <= 0.0 {
            0.0
        } else {
            v.sqrt()
        }
    }
}

pub struct ChromScoreInput<'a> {
    pub contacts: &'a [Contact],
    pub chrom_size: u64,
    pub res: u64,
    /// Distance filters, in base pairs (0 = disabled).
    pub min_distance: u64,
    pub max_distance: u64,
    /// Upper clamp applied to distance-normalised values.
    pub threshold: f64,
    /// Bin starts seeded as A.
    pub a_bins: &'a HashMap<u64, u32>,
    /// Bin starts seeded as B.
    pub b_bins: &'a HashMap<u64, u32>,
}

pub struct ChromScoreOutput {
    pub scores: Vec<BinScore>,
    pub a_moments: HashMap<u64, Moments>,
    pub b_moments: HashMap<u64, Moments>,
    /// Which bins carried any contact, indexed by `bin_start / res`.
    ///
    /// This is CRUSH's own row marginal — the same `row_sum` used for the row
    /// Z-scores, which is zero exactly when a bin has no contacts after the
    /// distance filters. Centromeres, telomeres and unmappable regions are
    /// false. The scoring stage already ignores these bins (`row_sd <= 0`), but
    /// the shifter works on a dense grid and would otherwise invent a score for
    /// them out of the neighbouring rolling mean.
    pub covered: Vec<bool>,
}

/// Distance-normalise, row-Z-score, and compute each bin's A-vs-B preference.
pub fn score_chromosome(input: &ChromScoreInput<'_>) -> ChromScoreOutput {
    let res = input.res;
    let n_bins = (input.chrom_size / res) as usize + 1;

    // ---- 1. Expected value per genomic distance -----------------------------
    // expected[d] = (sum of counts at distance d) / (number of bin pairs at d)
    // The bash pipeline uses (chrom_size/res - d) as the pair count, which is a
    // float division; reproduce that exactly rather than using n_bins - d.
    let mut dist_sum = vec![0f64; n_bins + 1];
    let bins_f = input.chrom_size as f64 / res as f64;

    let keep = |c: &Contact| -> bool {
        let d = (c.y as u64).abs_diff(c.x as u64);
        if d < input.min_distance {
            return false;
        }
        if input.max_distance > 0 && d > input.max_distance {
            return false;
        }
        true
    };

    for c in input.contacts.iter().filter(|c| keep(c)) {
        let d = ((c.y as u64).abs_diff(c.x as u64) / res) as usize;
        if d < dist_sum.len() {
            dist_sum[d] += c.value as f64;
        }
    }

    let expected: Vec<f64> = dist_sum
        .iter()
        .enumerate()
        .map(|(d, &s)| {
            let pairs = bins_f - d as f64;
            if pairs > 0.0 {
                s / pairs
            } else {
                0.0
            }
        })
        .collect();

    // ---- 2. O/E with pseudocount, clamped at `threshold` --------------------
    // The bash pipeline computes ($4+1)/(expected+1) and drops the diagonal.
    let mut oe: Vec<(u32, u32, f64)> = Vec::with_capacity(input.contacts.len());
    for c in input.contacts.iter().filter(|c| keep(c)) {
        if c.x == c.y {
            continue; // diagonal excluded ("if ($2 != $3)")
        }
        let d = ((c.y as u64).abs_diff(c.x as u64) / res) as usize;
        let e = expected.get(d).copied().unwrap_or(0.0);
        let mut v = (c.value as f64 + 1.0) / (e + 1.0);
        if v >= input.threshold {
            v = input.threshold;
        }
        oe.push((c.x, c.y, v));
    }

    // ---- 3. Row mean / SD over the symmetrised matrix -----------------------
    // Denominator is the full bin count, not the number of observed entries:
    // absent contacts count as structural zeros. Matches `-v rowbins`.
    let rowbins = (input.chrom_size / res) as f64 + 1.0;
    let mut row_sum = vec![0f64; n_bins];
    let mut row_sqsum = vec![0f64; n_bins];

    let mut accum = |bin: u32, v: f64| {
        let i = (bin as u64 / res) as usize;
        if i < n_bins {
            row_sum[i] += v;
            row_sqsum[i] += v * v;
        }
    };
    for &(x, y, v) in &oe {
        accum(x, v);
        accum(y, v); // symmetric partner, never materialised
    }

    let row_mean: Vec<f64> = row_sum.iter().map(|s| s / rowbins).collect();
    let row_sd: Vec<f64> = row_sqsum
        .iter()
        .zip(&row_mean)
        .map(|(sq, m)| {
            let var = sq / rowbins - m * m;
            if var <= 0.0 {
                0.0
            } else {
                var.sqrt()
            }
        })
        .collect();

    // ---- 4. Row Z-scores, partitioned by the partner bin's A/B state --------
    // For a contact (i, j) the bash pipeline emits both (i,j,z) and (j,i,z),
    // z-scored by the SECOND bin's row statistics, then keeps rows whose FIRST
    // bin is an A (or B) seed and aggregates the z by the second bin. Net
    // effect: bin j accumulates the z-scores of its contacts with A (or B)
    // seed bins. Generate both orientations without building the list.
    let mut a_moments: HashMap<u64, Moments> = HashMap::new();
    let mut b_moments: HashMap<u64, Moments> = HashMap::new();

    let mut push = |partner: u32, target: u32, v: f64,
                    a_m: &mut HashMap<u64, Moments>,
                    b_m: &mut HashMap<u64, Moments>| {
        let tj = (target as u64 / res) as usize;
        if tj >= n_bins || row_sd[tj] <= 0.0 {
            return;
        }
        let z = (v - row_mean[tj]) / row_sd[tj];
        let pbin = (partner as u64 / res) * res;
        let tbin = target as u64;
        if input.a_bins.contains_key(&pbin) {
            a_m.entry(tbin).or_default().push(z);
        }
        if input.b_bins.contains_key(&pbin) {
            b_m.entry(tbin).or_default().push(z);
        }
    };

    for &(x, y, v) in &oe {
        push(x, y, v, &mut a_moments, &mut b_moments);
        push(y, x, v, &mut a_moments, &mut b_moments);
    }

    // ---- 5. Diagonal pseudo-count -------------------------------------------
    // `zerofile` contributes (i, i, 0.1) for every bin, so a seed bin gets a
    // 0.1 z-score against itself. Reproduced verbatim.
    let mut i = 0u64;
    while i <= input.chrom_size {
        if input.a_bins.contains_key(&i) {
            a_moments.entry(i).or_default().push(0.1);
        }
        if input.b_bins.contains_key(&i) {
            b_moments.entry(i).or_default().push(0.1);
        }
        i += res;
    }

    // ---- 6. CRUSH score ------------------------------------------------------
    // score(j) = (Asum*Acount - Bsum*Bcount) / total_count
    let ab_sum: f64 = a_moments.values().map(|m| m.count).sum::<f64>()
        + b_moments.values().map(|m| m.count).sum::<f64>();

    let mut scores: Vec<BinScore> = Vec::new();
    if ab_sum > 0.0 {
        let mut all_bins: Vec<u64> = a_moments.keys().chain(b_moments.keys()).copied().collect();
        all_bins.sort_unstable();
        all_bins.dedup();

        for bin in all_bins {
            let a = a_moments.get(&bin).copied().unwrap_or_default();
            let b = b_moments.get(&bin).copied().unwrap_or_default();
            let value = (a.sum * a.count - b.sum * b.count) / ab_sum;
            if value != 0.0 {
                scores.push(BinScore {
                    bin_start: bin,
                    score: value,
                });
            }
        }
    }

    ChromScoreOutput {
        scores,
        a_moments,
        b_moments,
        covered: row_sum.iter().map(|&s| s > 0.0).collect(),
    }
}

/// Trailing-window mean matching the bash sliding-window awk, including its
/// coordinate recentring. `window` is the width in bins.
pub fn sliding_window(scores: &[BinScore], res: u64, window: usize) -> Vec<(u64, u64, f64)> {
    let w = window.max(1);
    let mut out = Vec::with_capacity(scores.len());
    let mut ring = vec![0f64; w];
    let mut sum = 0f64;
    let mut count = 0f64;

    for (n, s) in scores.iter().enumerate() {
        let nr = n + 1; // awk's NR is 1-based
        let m = nr % w;
        if nr <= w {
            count += 1.0;
        } else {
            sum -= ring[m];
        }
        sum += s.score;
        ring[m] = s.score;

        // Centre the window, then apply the +res bin offset the bash pipeline
        // adds when it converts a point into an interval.
        let centred = s.bin_start as i64 - ((res * w as u64) / 2) as i64;
        let start = centred - (res / 2) as i64 + res as i64;
        let end = centred + (res / 2) as i64 + res as i64;
        if start > 0 {
            out.push((start as u64, end as u64, sum / count));
        }
    }
    out
}
