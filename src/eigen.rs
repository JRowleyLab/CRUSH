//! Eigenvector initialisation of the A/B axis.
//!
//! Same strategy as the bash pipeline: build the O/E matrix at `eigenres`,
//! Pearson-correlate it, take the top principal components, orient each by gene
//! content, then pick the component most correlated with gene density and least
//! correlated with the B-seed (low-GC) density.
//!
//! Two differences from the bash implementation, both deliberate:
//!   * only the leading `n_components` eigenpairs are needed, so this uses
//!     subspace iteration rather than a full symmetric eigendecomposition —
//!     O(n^2 k) instead of O(n^3);
//!   * empty rows are dropped before correlating (as in bash) but the mapping
//!     back to genomic coordinates is kept exact.

use std::collections::HashMap;

use crate::hic::Contact;

pub struct EigenResult {
    /// (bin_start, value) for the selected principal component.
    pub track: Vec<(u64, f64)>,
    pub chosen_pc: usize,
    pub gene_corr: f64,
    pub gc_corr: f64,
}

/// Pearson correlation; `None` if either input is constant.
pub fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    let n = a.len().min(b.len());
    if n < 2 {
        return None;
    }
    let ma = a[..n].iter().sum::<f64>() / n as f64;
    let mb = b[..n].iter().sum::<f64>() / n as f64;
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let (x, y) = (a[i] - ma, b[i] - mb);
        num += x * y;
        da += x * x;
        db += y * y;
    }
    if da == 0.0 || db == 0.0 {
        return None;
    }
    Some(num / (da.sqrt() * db.sqrt()))
}

/// Row-wise Pearson correlation matrix of a dense square matrix.
fn correlation_matrix(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = m.len();
    // Pre-centre and pre-scale each row so the correlation is a plain dot product.
    let mut norm: Vec<Vec<f64>> = Vec::with_capacity(n);
    for row in m {
        let mean = row.iter().sum::<f64>() / row.len() as f64;
        let mut r: Vec<f64> = row.iter().map(|v| v - mean).collect();
        let ss = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if ss > 0.0 {
            for v in r.iter_mut() {
                *v /= ss;
            }
        }
        norm.push(r);
    }

    let mut c = vec![vec![0f64; n]; n];
    for i in 0..n {
        c[i][i] = 1.0;
        for j in (i + 1)..n {
            let dot: f64 = norm[i].iter().zip(&norm[j]).map(|(a, b)| a * b).sum();
            // NaN (constant row) is treated as zero correlation, as np.nan_to_num does.
            let v = if dot.is_finite() { dot } else { 0.0 };
            c[i][j] = v;
            c[j][i] = v;
        }
    }
    c
}

/// Leading `k` eigenvectors of a symmetric matrix via orthogonal (subspace)
/// iteration with Gram-Schmidt re-orthogonalisation.
fn top_eigenvectors(a: &[Vec<f64>], k: usize, iters: usize) -> Vec<Vec<f64>> {
    let n = a.len();
    let k = k.min(n);
    if k == 0 {
        return Vec::new();
    }

    // Deterministic starting basis: no RNG, so runs are reproducible.
    let mut q: Vec<Vec<f64>> = (0..k)
        .map(|c| {
            (0..n)
                .map(|r| {
                    let x = ((r * 2654435761usize).wrapping_add(c * 40503) % 1000) as f64;
                    x / 1000.0 - 0.5 + if r % (c + 1) == 0 { 0.25 } else { 0.0 }
                })
                .collect()
        })
        .collect();
    orthonormalise(&mut q);

    for _ in 0..iters {
        let mut z: Vec<Vec<f64>> = Vec::with_capacity(k);
        for v in &q {
            let mut av = vec![0f64; n];
            for (i, row) in a.iter().enumerate() {
                let mut s = 0.0;
                for (j, &aij) in row.iter().enumerate() {
                    s += aij * v[j];
                }
                av[i] = s;
            }
            z.push(av);
        }
        orthonormalise(&mut z);
        q = z;
    }

    // Order by Rayleigh quotient, descending, so index 0 is PC1.
    let mut with_val: Vec<(f64, Vec<f64>)> = q
        .into_iter()
        .map(|v| {
            let mut av = vec![0f64; n];
            for (i, row) in a.iter().enumerate() {
                av[i] = row.iter().zip(&v).map(|(x, y)| x * y).sum();
            }
            let lambda: f64 = v.iter().zip(&av).map(|(x, y)| x * y).sum();
            (lambda, v)
        })
        .collect();
    with_val.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap_or(std::cmp::Ordering::Equal));
    with_val.into_iter().map(|(_, v)| v).collect()
}

fn orthonormalise(vecs: &mut Vec<Vec<f64>>) {
    for i in 0..vecs.len() {
        for j in 0..i {
            let dot: f64 = vecs[i].iter().zip(&vecs[j]).map(|(a, b)| a * b).sum();
            for t in 0..vecs[i].len() {
                vecs[i][t] -= dot * vecs[j][t];
            }
        }
        let nrm: f64 = vecs[i].iter().map(|v| v * v).sum::<f64>().sqrt();
        if nrm > 1e-12 {
            for v in vecs[i].iter_mut() {
                *v /= nrm;
            }
        }
    }
}

pub struct EigenInput<'a> {
    /// O/E contacts at `res` for this chromosome.
    pub contacts: &'a [Contact],
    pub res: u64,
    pub chrom_size: u64,
    /// Bin occupancy of the gene annotation (initial A).
    pub gene_bins: &'a HashMap<u64, u32>,
    /// Bin occupancy of the B seeds (low-GC bins).
    pub gc_bins: &'a HashMap<u64, u32>,
    pub n_components: usize,
    /// Refuse to allocate a dense matrix bigger than this, in bytes.
    pub max_matrix_bytes: usize,
}

pub fn compute_eigen(input: &EigenInput<'_>) -> Result<EigenResult, String> {
    let res = input.res;
    let max_bin = input
        .contacts
        .iter()
        .map(|c| c.x.max(c.y) as u64)
        .max()
        .unwrap_or(0);
    let n_bins = (max_bin / res) as usize + 1;
    if n_bins < 2 {
        return Err("fewer than two bins with data".into());
    }

    let need = n_bins
        .checked_mul(n_bins)
        .and_then(|v| v.checked_mul(std::mem::size_of::<f64>()))
        .ok_or_else(|| "matrix size overflow".to_string())?;
    if need > input.max_matrix_bytes {
        return Err(format!(
            "eigenvector matrix would need {:.1} GB ({} bins at {} bp); \
             use a coarser --eigenres",
            need as f64 / (1024.0 * 1024.0 * 1024.0),
            n_bins,
            res
        ));
    }

    let mut m = vec![vec![0f64; n_bins]; n_bins];
    for c in input.contacts {
        let (i, j) = ((c.x as u64 / res) as usize, (c.y as u64 / res) as usize);
        if i < n_bins && j < n_bins {
            m[i][j] = c.value as f64;
            if i != j {
                m[j][i] = c.value as f64;
            }
        }
    }

    // Drop all-zero rows/columns before correlating, as the bash version does.
    let keep: Vec<usize> = (0..n_bins)
        .filter(|&i| m[i].iter().any(|&v| v != 0.0))
        .collect();
    if keep.len() < 2 {
        return Err("fewer than two non-empty bins".into());
    }

    let sub: Vec<Vec<f64>> = keep
        .iter()
        .map(|&i| keep.iter().map(|&j| m[i][j]).collect())
        .collect();
    drop(m);

    let corr = correlation_matrix(&sub);
    let pcs = top_eigenvectors(&corr, input.n_components, 64);

    // Densities on the retained bins only, so all vectors line up.
    let gene_d: Vec<f64> = keep
        .iter()
        .map(|&i| *input.gene_bins.get(&((i as u64) * res)).unwrap_or(&0) as f64)
        .collect();
    let gc_d: Vec<f64> = keep
        .iter()
        .map(|&i| *input.gc_bins.get(&((i as u64) * res)).unwrap_or(&0) as f64)
        .collect();

    let mut best: Option<EigenResult> = None;
    for (idx, pc) in pcs.iter().enumerate() {
        // Orient so that gene-rich bins are positive.
        let mut v = pc.clone();
        let pos: f64 = v
            .iter()
            .zip(&gene_d)
            .filter(|(x, _)| **x > 0.0)
            .map(|(_, g)| *g)
            .sum();
        let neg: f64 = v
            .iter()
            .zip(&gene_d)
            .filter(|(x, _)| **x < 0.0)
            .map(|(_, g)| *g)
            .sum();
        if neg > pos {
            for x in v.iter_mut() {
                *x = -*x;
            }
        }

        let (Some(gc_r), Some(gene_r)) = (pearson(&v, &gc_d), pearson(&v, &gene_d)) else {
            continue;
        };
        let total = gene_r - gc_r;
        if best
            .as_ref()
            .map_or(true, |b| total > b.gene_corr - b.gc_corr)
        {
            best = Some(EigenResult {
                track: keep
                    .iter()
                    .zip(&v)
                    .map(|(&i, &val)| ((i as u64) * res, val))
                    .collect(),
                chosen_pc: idx + 1,
                gene_corr: gene_r,
                gc_corr: gc_r,
            });
        }
    }

    best.ok_or_else(|| "no principal component had a defined correlation".to_string())
}
