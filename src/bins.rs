/// Per-reference accumulator: bins[i] holds the count of aligned reference bases
/// whose coordinate falls inside bin i (i.e. summed depth contributions).
/// Mean depth = bins[i] as f64 / bin_size, except the last (possibly short) bin
/// which uses its actual length.
pub struct BinMatrix {
    pub bin_size: u32,
    pub refs: Vec<RefBins>,
}

pub struct RefBins {
    pub length: u32,
    pub counts: Vec<u64>,
}

impl BinMatrix {
    pub fn new(reference_lengths: &[u32], bin_size: u32) -> Self {
        assert!(bin_size > 0, "bin_size must be > 0");
        let refs = reference_lengths
            .iter()
            .map(|&len| {
                let n = ((len as u64 + bin_size as u64 - 1) / bin_size as u64) as usize;
                RefBins {
                    length: len,
                    counts: vec![0u64; n],
                }
            })
            .collect();
        Self { bin_size, refs }
    }

    /// Add `len` reference-consumed bases starting at 0-based `start` on reference `ref_id`,
    /// distributing across overlapping bins. Out-of-range portions are clipped silently.
    pub fn add_interval(&mut self, ref_id: usize, start: u32, len: u32) {
        let Some(r) = self.refs.get_mut(ref_id) else {
            return;
        };
        if len == 0 || start >= r.length {
            return;
        }
        let end = start.saturating_add(len).min(r.length); // exclusive
        let bin_size = self.bin_size;
        let mut pos = start;
        let mut bin_idx = (pos / bin_size) as usize;
        while pos < end {
            let bin_end = ((bin_idx as u32).saturating_add(1)) * bin_size;
            let stop = bin_end.min(end);
            r.counts[bin_idx] += (stop - pos) as u64;
            pos = stop;
            bin_idx += 1;
        }
    }

    pub fn bin_span(&self, ref_idx: usize, bin_idx: usize) -> (u32, u32) {
        let start = bin_idx as u32 * self.bin_size;
        let end = ((bin_idx as u32 + 1) * self.bin_size).min(self.refs[ref_idx].length);
        (start, end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_bin_full_coverage() {
        let mut m = BinMatrix::new(&[1000], 100);
        m.add_interval(0, 0, 100);
        assert_eq!(m.refs[0].counts[0], 100);
        assert_eq!(m.refs[0].counts[1], 0);
    }

    #[test]
    fn spans_three_bins() {
        let mut m = BinMatrix::new(&[1000], 100);
        m.add_interval(0, 50, 200); // 50..250
        assert_eq!(m.refs[0].counts[0], 50); // 50..100
        assert_eq!(m.refs[0].counts[1], 100); // 100..200
        assert_eq!(m.refs[0].counts[2], 50); // 200..250
    }

    #[test]
    fn clipped_at_chrom_end() {
        let mut m = BinMatrix::new(&[150], 100);
        m.add_interval(0, 120, 100); // 120..220 -> clipped to 120..150
        assert_eq!(m.refs[0].counts[1], 30);
    }
}
