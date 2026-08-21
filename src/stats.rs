use crate::bins::BinMatrix;
use noodles_sam::Header;

pub struct RefSummary {
    pub name: String,
    pub length: u32,
    pub n_bins: usize,
    pub total_aligned_bases: u64,
    pub mean_depth: f64,
    pub min_depth: f64,
    pub max_depth: f64,
    pub median_depth: f64,
    pub stdev_depth: f64,
    pub breadth_ge1: f64,
    pub breadth_ge5: f64,
    pub breadth_ge10: f64,
    pub breadth_ge20: f64,
    pub breadth_ge30: f64,
}

pub struct Summary {
    pub per_ref: Vec<RefSummary>,
    pub total: RefSummary,
}

pub fn summarize(header: &Header, matrix: &BinMatrix, bin_size: u32) -> Summary {
    let names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    let per_ref: Vec<RefSummary> = matrix
        .refs
        .iter()
        .zip(names.iter())
        .map(|(rb, name)| ref_summary(name, rb.length, &rb.counts, bin_size))
        .collect();

    let mut total_bins = Vec::new();
    let mut total_len: u64 = 0;
    for rb in &matrix.refs {
        total_bins.extend_from_slice(&rb.counts);
        total_len += rb.length as u64;
    }
    let total = ref_summary(
        "total",
        total_len.min(u32::MAX as u64) as u32,
        &total_bins,
        bin_size,
    );

    Summary { per_ref, total }
}

fn ref_summary(name: &str, length: u32, counts: &[u64], bin_size: u32) -> RefSummary {
    let n_bins = counts.len();
    if n_bins == 0 {
        return RefSummary {
            name: name.to_string(),
            length,
            n_bins: 0,
            total_aligned_bases: 0,
            mean_depth: 0.0,
            min_depth: 0.0,
            max_depth: 0.0,
            median_depth: 0.0,
            stdev_depth: 0.0,
            breadth_ge1: 0.0,
            breadth_ge5: 0.0,
            breadth_ge10: 0.0,
            breadth_ge20: 0.0,
            breadth_ge30: 0.0,
        };
    }

    let mut depths: Vec<f64> = counts
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            let span = bin_span_len(i, n_bins, length, bin_size) as f64;
            if span > 0.0 {
                c as f64 / span
            } else {
                0.0
            }
        })
        .collect();

    let total_aligned_bases: u64 = counts.iter().sum();
    let mean_depth = total_aligned_bases as f64 / length.max(1) as f64;

    let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut sumsq = 0.0;
    for &d in &depths {
        if d < min {
            min = d;
        }
        if d > max {
            max = d;
        }
        sumsq += (d - mean_depth).powi(2);
    }
    let stdev_depth = (sumsq / n_bins as f64).sqrt();

    let mid = n_bins / 2;
    let (_, m, _) = depths.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap());
    let median_depth = *m;

    let mut c1 = 0usize;
    let mut c5 = 0usize;
    let mut c10 = 0usize;
    let mut c20 = 0usize;
    let mut c30 = 0usize;
    for &d in &depths {
        if d >= 1.0 {
            c1 += 1;
        }
        if d >= 5.0 {
            c5 += 1;
        }
        if d >= 10.0 {
            c10 += 1;
        }
        if d >= 20.0 {
            c20 += 1;
        }
        if d >= 30.0 {
            c30 += 1;
        }
    }
    let n = n_bins as f64;

    RefSummary {
        name: name.to_string(),
        length,
        n_bins,
        total_aligned_bases,
        mean_depth,
        min_depth: if min.is_finite() { min } else { 0.0 },
        max_depth: if max.is_finite() { max } else { 0.0 },
        median_depth,
        stdev_depth,
        breadth_ge1: c1 as f64 / n,
        breadth_ge5: c5 as f64 / n,
        breadth_ge10: c10 as f64 / n,
        breadth_ge20: c20 as f64 / n,
        breadth_ge30: c30 as f64 / n,
    }
}

fn bin_span_len(i: usize, n_bins: usize, length: u32, bin_size: u32) -> u32 {
    if i + 1 < n_bins {
        bin_size
    } else {
        length - (i as u32 * bin_size)
    }
}
