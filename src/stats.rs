//! Welch statistic, normal survival function, and Benjamini-Hochberg FDR.
//!
//! The bash pipeline shells out to scipy and statsmodels for these. They are a
//! few dozen lines each, so implementing them directly removes two Python
//! process launches per chromosome and the statsmodels dependency entirely.

use crate::score::Moments;

/// Complementary error function, Numerical Recipes rational approximation.
/// Relative error < 1.2e-7, which is well inside the precision that survives
/// the text round-trip to bedGraph.
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 2.0 / (2.0 + z);
    let ty = 4.0 * t - 2.0;

    const C: [f64; 28] = [
        -1.3026537197817094,
        6.419_697_923_564_902e-1,
        1.9476473204185836e-2,
        -9.561_514_786_808_631e-3,
        -9.46595344482036e-4,
        3.66839497852761e-4,
        4.2523324806907e-5,
        -2.0278578112534e-5,
        -1.624290004647e-6,
        1.303655835580e-6,
        1.5626441722e-8,
        -8.5238095915e-8,
        6.529054439e-9,
        5.059343495e-9,
        -9.91364156e-10,
        -2.27365122e-10,
        9.6467911e-11,
        2.394038e-12,
        -6.886027e-12,
        8.94487e-13,
        3.13092e-13,
        -1.12708e-13,
        3.81e-16,
        7.106e-15,
        -1.523e-15,
        -9.4e-17,
        1.21e-16,
        -2.8e-17,
    ];

    let (mut d, mut dd) = (0.0f64, 0.0f64);
    for &c in C.iter().rev().take(C.len() - 1) {
        let tmp = d;
        d = ty * d - dd + c;
        dd = tmp;
    }
    let ans = t * (-z * z + 0.5 * (C[0] + ty * d) - dd).exp();
    if x >= 0.0 {
        ans
    } else {
        2.0 - ans
    }
}

/// Upper-tail probability of the standard normal, i.e. `scipy.stats.norm.sf`.
pub fn norm_sf(x: f64) -> f64 {
    0.5 * erfc(x / std::f64::consts::SQRT_2)
}

/// Two-sample Welch statistic between the A and B moment sets for one bin.
///
/// Returns `None` when either group has a single observation, matching the
/// bash guard `($3 != 1) && (c[$1] != 1)`.
pub fn welch_p(a: &Moments, b: &Moments) -> Option<f64> {
    if a.count <= 1.0 || b.count <= 1.0 {
        return None;
    }
    let (sa, sb) = (a.sd(), b.sd());
    let denom = ((sa * sa) / a.count + (sb * sb) / b.count).sqrt();
    if denom == 0.0 || !denom.is_finite() {
        return None;
    }
    let t = (a.mean() - b.mean()) / denom;
    if !t.is_finite() {
        return None;
    }
    Some(norm_sf(t.abs()))
}

/// Benjamini-Hochberg step-up FDR, equivalent to
/// `statsmodels.stats.multitest.multipletests(..., method='fdr_bh')`.
/// Returns q-values in the input order.
pub fn benjamini_hochberg(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    if n == 0 {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| p[i].partial_cmp(&p[j]).unwrap_or(std::cmp::Ordering::Equal));

    let mut q = vec![0f64; n];
    let mut running = f64::INFINITY;
    // Walk from largest p to smallest, enforcing monotonicity.
    for rank in (1..=n).rev() {
        let idx = order[rank - 1];
        let val = p[idx] * n as f64 / rank as f64;
        running = running.min(val);
        q[idx] = running.min(1.0);
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn norm_sf_matches_known_values() {
        // scipy.stats.norm.sf reference values
        assert!((norm_sf(0.0) - 0.5).abs() < 1e-9);
        assert!((norm_sf(1.0) - 0.158_655_253_9).abs() < 1e-7);
        assert!((norm_sf(1.959_963_985) - 0.025).abs() < 1e-7);
        assert!((norm_sf(3.0) - 0.001_349_898_0).abs() < 1e-7);
        assert!((norm_sf(5.0) - 2.866_515_719e-7).abs() < 1e-12);
    }

    #[test]
    fn bh_matches_statsmodels() {
        // multipletests([0.01,0.02,0.03,0.04,0.05], method='fdr_bh')[1]
        let p = [0.01, 0.02, 0.03, 0.04, 0.05];
        let q = benjamini_hochberg(&p);
        let expect = [0.05, 0.05, 0.05, 0.05, 0.05];
        for (a, b) in q.iter().zip(expect.iter()) {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
    }

    #[test]
    fn bh_is_monotone_and_unsorted_safe() {
        let p = [0.9, 0.001, 0.5, 0.02];
        let q = benjamini_hochberg(&p);
        assert!(q.iter().all(|v| (0.0..=1.0).contains(v)));
        // smallest p must get the smallest q
        assert!(q[1] < q[3] && q[3] < q[2] && q[2] <= q[0]);
    }
}
