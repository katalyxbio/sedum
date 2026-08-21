use crate::bins::BinMatrix;
use crate::cigar::visit_ref_intervals;
use crate::debug::ScanCounters;
use crate::filters::Filter;
use crate::plan::{self, Region};
use crate::progress::Progress;
use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_core::Region as CoreRegion;
use noodles_sam::Header;
use std::path::Path;
use std::sync::{Arc, Mutex};

const RECORD_BATCH: u64 = 4096;

pub fn scan_single_threaded(
    path: &Path,
    bin_size: u32,
    filter: &Filter,
    progress: &Progress,
) -> Result<(Header, BinMatrix, ScanCounters)> {
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("opening BAM {}", path.display()))?;
    let header = reader.read_header().context("reading BAM header")?;
    let ref_lengths = ref_lengths(&header);
    let mut matrix = BinMatrix::new(&ref_lengths, bin_size);
    let mut counters = ScanCounters::default();

    let mut record = bam::Record::default();
    let mut local: u64 = 0;
    while reader.read_record(&mut record)? != 0 {
        process_record(&record, filter, None, &mut matrix, &mut counters);
        local += 1;
        if local == RECORD_BATCH {
            progress.inc_records(local);
            local = 0;
        }
    }
    progress.inc_records(local);
    Ok((header, matrix, counters))
}

pub fn scan_parallel(
    path: &Path,
    bin_size: u32,
    filter: &Filter,
    threads: usize,
    stripes_per_thread: usize,
    progress: &Arc<Progress>,
) -> Result<(Header, BinMatrix, ScanCounters)> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("opening indexed BAM {}", path.display()))?;
    let header = reader.read_header().context("reading BAM header")?;

    let ref_seqs: Vec<(String, u32)> = header
        .reference_sequences()
        .iter()
        .map(|(name, rs)| (name.to_string(), usize::from(rs.length()) as u32))
        .collect();
    let ref_lengths: Vec<u32> = ref_seqs.iter().map(|(_, l)| *l).collect();

    let regions = plan::plan(&ref_seqs, threads, stripes_per_thread);
    progress.set_stripes_total(regions.len() as u64);
    let queue: Arc<Mutex<Vec<Region>>> = Arc::new(Mutex::new(regions));
    let header_ref = &header;
    let path_ref = path;
    let filter_owned = *filter;
    let ref_lengths_ref = &ref_lengths;

    let shards: Result<Vec<(BinMatrix, ScanCounters)>> = std::thread::scope(|s| {
        let mut handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            let q = Arc::clone(&queue);
            let prog = Arc::clone(progress);
            let h = s.spawn(move || -> Result<(BinMatrix, ScanCounters)> {
                let mut shard = BinMatrix::new(ref_lengths_ref, bin_size);
                let mut counters = ScanCounters::default();
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(path_ref)
                    .with_context(|| {
                        format!("worker opening indexed BAM {}", path_ref.display())
                    })?;
                let _ = reader.read_header()?;
                let mut local: u64 = 0;

                loop {
                    let region = {
                        let mut g = q.lock().unwrap();
                        g.pop()
                    };
                    let Some(region) = region else { break };
                    if region.end <= region.start {
                        prog.inc_stripe();
                        continue;
                    }
                    let region_str = format!(
                        "{}:{}-{}",
                        region.name,
                        region.start as usize + 1,
                        region.end as usize
                    );
                    let q_region: CoreRegion = region_str
                        .parse()
                        .with_context(|| format!("parsing region {region_str}"))?;
                    let query = reader
                        .query(header_ref, &q_region)
                        .with_context(|| format!("querying {region_str}"))?;
                    for rec_res in query {
                        let record = rec_res?;
                        process_record(
                            &record,
                            &filter_owned,
                            Some((region.start, region.end)),
                            &mut shard,
                            &mut counters,
                        );
                        local += 1;
                        if local == RECORD_BATCH {
                            prog.inc_records(local);
                            local = 0;
                        }
                    }
                    prog.inc_stripe();
                }
                prog.inc_records(local);
                Ok((shard, counters))
            });
            handles.push(h);
        }
        handles
            .into_iter()
            .map(|h| h.join().expect("worker thread panicked"))
            .collect()
    });
    let shards = shards?;

    let mut merged = BinMatrix::new(&ref_lengths, bin_size);
    let mut counters = ScanCounters::default();
    for (shard, shard_counters) in shards {
        counters.merge(&shard_counters);
        for (mref, sref) in merged.refs.iter_mut().zip(shard.refs.iter()) {
            for (a, b) in mref.counts.iter_mut().zip(sref.counts.iter()) {
                *a += *b;
            }
        }
    }
    Ok((header, merged, counters))
}

fn ref_lengths(header: &Header) -> Vec<u32> {
    header
        .reference_sequences()
        .iter()
        .map(|(_, rs)| usize::from(rs.length()) as u32)
        .collect()
}

#[inline]
fn process_record(
    record: &bam::Record,
    filter: &Filter,
    stripe: Option<(u32, u32)>,
    matrix: &mut BinMatrix,
    counters: &mut ScanCounters,
) {
    counters.records_visited += 1;
    let flags = u16::from(record.flags());
    if !filter.flags_ok(flags) {
        counters.filtered_flags += 1;
        return;
    }
    let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
    if !filter.mapq_ok(mapq) {
        counters.filtered_mapq += 1;
        return;
    }
    let Some(Ok(ref_id)) = record.reference_sequence_id() else {
        counters.unmapped += 1;
        return;
    };
    let Some(Ok(start)) = record.alignment_start() else {
        counters.unmapped += 1;
        return;
    };
    let start_0 = (usize::from(start) - 1) as u32;
    if let Some((rstart, rend)) = stripe {
        if start_0 < rstart || start_0 >= rend {
            counters.out_of_stripe += 1;
            return;
        }
    }
    counters.kept += 1;
    let cigar = record.cigar();
    let ops_iter = cigar
        .iter()
        .filter_map(|res| res.ok().map(|op| (op.kind(), op.len() as u32)));
    visit_ref_intervals(start_0, ops_iter, filter.count_splice, |s, l| {
        counters.aligned_bases += l as u64;
        matrix.add_interval(ref_id, s, l);
    });
}
