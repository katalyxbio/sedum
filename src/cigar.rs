use noodles_sam::alignment::record::cigar::op::Kind;

/// Visit each reference-consuming aligned interval as (ref_start_0based, length).
///
/// `count_splice` controls whether CIGAR N (skipped reference, typical of spliced
/// RNA-seq alignments) contributes to depth. Soft/hard clips and insertions never do.
pub fn visit_ref_intervals<I, F>(start_0based: u32, ops: I, count_splice: bool, mut emit: F)
where
    I: IntoIterator<Item = (Kind, u32)>,
    F: FnMut(u32, u32),
{
    let mut pos = start_0based;
    for (kind, len) in ops {
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                if len > 0 {
                    emit(pos, len);
                }
                pos = pos.saturating_add(len);
            }
            Kind::Skip => {
                if count_splice && len > 0 {
                    emit(pos, len);
                }
                pos = pos.saturating_add(len);
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
}
