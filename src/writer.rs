use crate::bins::BinMatrix;
use crate::stats::{RefSummary, Summary};
use anyhow::{Context, Result};
use noodles_sam::Header;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn write_bins(prefix: &str, header: &Header, matrix: &BinMatrix, bin_size: u32) -> Result<()> {
    let path = format!("{prefix}.bins.bed");
    let mut w = BufWriter::new(File::create(&path).with_context(|| format!("creating {path}"))?);

    let names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(n, _)| n.to_string())
        .collect();

    for (ref_idx, rb) in matrix.refs.iter().enumerate() {
        let name = &names[ref_idx];
        for (bin_idx, &count) in rb.counts.iter().enumerate() {
            let (start, end) = matrix.bin_span(ref_idx, bin_idx);
            let span = (end - start) as f64;
            let mean = if span > 0.0 { count as f64 / span } else { 0.0 };
            writeln!(w, "{name}\t{start}\t{end}\t{mean:.4}")?;
        }
    }
    let _ = bin_size;
    Ok(())
}

pub fn write_summary(prefix: &str, summary: &Summary) -> Result<()> {
    let path = format!("{prefix}.summary.tsv");
    let mut w = BufWriter::new(File::create(&path).with_context(|| format!("creating {path}"))?);

    writeln!(
        w,
        "chrom\tlength\tn_bins\ttotal_aligned_bases\tmean\tmin\tmax\tmedian\tstdev\tbreadth_ge1\tbreadth_ge5\tbreadth_ge10\tbreadth_ge20\tbreadth_ge30"
    )?;
    for s in &summary.per_ref {
        write_row(&mut w, s)?;
    }
    write_row(&mut w, &summary.total)?;
    Ok(())
}

fn write_row<W: Write>(w: &mut W, s: &RefSummary) -> Result<()> {
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        s.name,
        s.length,
        s.n_bins,
        s.total_aligned_bases,
        s.mean_depth,
        s.min_depth,
        s.max_depth,
        s.median_depth,
        s.stdev_depth,
        s.breadth_ge1,
        s.breadth_ge5,
        s.breadth_ge10,
        s.breadth_ge20,
        s.breadth_ge30,
    )?;
    Ok(())
}
