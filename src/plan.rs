/// A half-open genomic stripe owned by a single worker.
/// Reads whose `alignment_start` (0-based) falls in `[start, end)` are this
/// stripe's responsibility; their CIGAR may extend past `end`, and that's fine
/// — it just lands in bins owned by the next stripe's accumulator and gets
/// summed during the final reduction.
#[derive(Clone, Debug)]
pub struct Region {
    #[allow(dead_code)]
    pub ref_id: usize,
    pub name: String,
    pub start: u32,
    pub end: u32,
}

pub fn plan(ref_seqs: &[(String, u32)], threads: usize, stripes_per_thread: usize) -> Vec<Region> {
    let threads = threads.max(1);
    let spt = stripes_per_thread.max(1);
    let total: u64 = ref_seqs.iter().map(|(_, l)| *l as u64).sum();
    let target_stripes = (threads as u64) * (spt as u64);
    let stripe_size = if total == 0 {
        1
    } else {
        total.div_ceil(target_stripes).max(1) as u32
    };

    let mut out = Vec::new();
    for (ref_id, (name, length)) in ref_seqs.iter().enumerate() {
        if *length == 0 {
            continue;
        }
        let mut pos = 0u32;
        while pos < *length {
            let end = pos.saturating_add(stripe_size).min(*length);
            out.push(Region {
                ref_id,
                name: name.clone(),
                start: pos,
                end,
            });
            pos = end;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stripes_cover_each_reference_exactly_once() {
        let refs = vec![("chr1".to_string(), 10_000u32), ("chr2".to_string(), 3_500)];
        let regions = plan(&refs, 4, 4);
        // Per ref, stripes must be contiguous, non-overlapping, covering [0, length).
        for (rid, (_, len)) in refs.iter().enumerate() {
            let mut stripes: Vec<&Region> = regions.iter().filter(|r| r.ref_id == rid).collect();
            stripes.sort_by_key(|r| r.start);
            assert_eq!(stripes.first().unwrap().start, 0);
            assert_eq!(stripes.last().unwrap().end, *len);
            for w in stripes.windows(2) {
                assert_eq!(w[0].end, w[1].start);
            }
        }
    }

    #[test]
    fn zero_length_references_skipped() {
        let refs = vec![("empty".to_string(), 0u32), ("chr1".to_string(), 100)];
        let regions = plan(&refs, 2, 4);
        assert!(regions.iter().all(|r| r.ref_id == 1));
    }
}
