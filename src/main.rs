mod bins;
mod cigar;
mod cli;
mod filters;
mod plan;
mod progress;
mod reader;
mod stats;
mod writer;

use anyhow::Result;
use clap::Parser;

const NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

fn main() -> Result<()> {
    let args = cli::Cli::parse();
    run(args)
}

fn run(args: cli::Cli) -> Result<()> {
    let started = std::time::Instant::now();
    let filter = filters::Filter::from_cli(&args);
    let threads = args.effective_threads();
    let use_parallel = !args.no_parallel && threads > 1 && bai_present(&args.bam);

    eprintln!("{NAME} v{VERSION}");
    eprintln!("  input        {}", args.bam.display());
    eprintln!("  bin size     {} bp", args.bin_size);
    eprintln!(
        "  filters      min_mapq={}, exclude=0x{:x}, include=0x{:x}, count_splice={}",
        filter.min_mapq, filter.exclude_flags, filter.include_flags, filter.count_splice
    );
    if use_parallel {
        eprintln!(
            "  scan         parallel, threads={threads}, stripes/thread={}",
            args.stripes_per_thread
        );
    } else {
        eprintln!("  scan         single-threaded");
    }
    eprintln!("  output       {}.bins.bed, {}.summary.tsv", args.output_prefix, args.output_prefix);
    eprintln!();

    let prog = progress::Progress::new(!args.no_progress);
    let reporter = progress::Reporter::spawn(std::sync::Arc::clone(&prog));

    let (header, matrix) = if use_parallel {
        reader::scan_parallel(
            &args.bam,
            args.bin_size,
            &filter,
            threads,
            args.stripes_per_thread,
            &prog,
        )?
    } else {
        reader::scan_single_threaded(&args.bam, args.bin_size, &filter, &prog)?
    };

    drop(reporter); // flushes the final progress line before we start writing.

    let summary = stats::summarize(&header, &matrix, args.bin_size);
    writer::write_bins(&args.output_prefix, &header, &matrix, args.bin_size)?;
    writer::write_summary(&args.output_prefix, &summary)?;

    print_final_summary(&summary, started.elapsed());
    Ok(())
}

fn print_final_summary(summary: &stats::Summary, elapsed: std::time::Duration) {
    let t = &summary.total;
    let n_refs_with_coverage = summary
        .per_ref
        .iter()
        .filter(|r| r.total_aligned_bases > 0)
        .count();
    eprintln!();
    eprintln!("genome summary");
    eprintln!(
        "  references           {} ({} with coverage)",
        summary.per_ref.len(),
        n_refs_with_coverage
    );
    eprintln!("  total length         {} bp", t.length);
    eprintln!("  total bins           {}", t.n_bins);
    eprintln!("  aligned bases        {}", t.total_aligned_bases);
    eprintln!("  mean depth           {:.4}x", t.mean_depth);
    eprintln!(
        "  bin depth            min={:.4}  median={:.4}  max={:.4}  stdev={:.4}",
        t.min_depth, t.median_depth, t.max_depth, t.stdev_depth
    );
    eprintln!(
        "  breadth (>= 1x/5x/10x/20x/30x)  {:.2}% / {:.2}% / {:.2}% / {:.2}% / {:.2}%",
        100.0 * t.breadth_ge1,
        100.0 * t.breadth_ge5,
        100.0 * t.breadth_ge10,
        100.0 * t.breadth_ge20,
        100.0 * t.breadth_ge30
    );
    eprintln!("  wall time            {:.2}s", elapsed.as_secs_f64());
}

fn bai_present(bam: &std::path::Path) -> bool {
    let with_ext = {
        let mut p = bam.to_path_buf();
        let new_ext = match p.extension().and_then(|e| e.to_str()) {
            Some(ext) => format!("{ext}.bai"),
            None => "bai".to_string(),
        };
        p.set_extension(new_ext);
        p
    };
    let sibling = bam.with_extension("bai");
    with_ext.exists() || sibling.exists()
}
