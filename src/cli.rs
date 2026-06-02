use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "sedum", version, about = "Fast per-bin depth calculator")]
pub struct Cli {
    /// Input BAM (indexed BAM recommended; SAM accepted)
    pub bam: PathBuf,

    /// Bin size in base pairs
    #[arg(short = 'b', long, default_value_t = 500)]
    pub bin_size: u32,

    /// Minimum mapping quality
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Exclude reads with any of these SAM flags set (hex or decimal).
    /// Default 0x704 = UNMAP|SECONDARY|QCFAIL|DUP
    #[arg(long, default_value = "0x704", value_parser = parse_flag)]
    pub exclude_flags: u16,

    /// Require these SAM flags set (hex or decimal). Default 0
    #[arg(long, default_value = "0", value_parser = parse_flag)]
    pub include_flags: u16,

    /// Count reference skips (CIGAR N) toward depth (off by default; RNA-seq use case)
    #[arg(long, default_value_t = false)]
    pub count_splice: bool,

    /// Output prefix; writes {prefix}.bins.bed and {prefix}.summary.tsv
    #[arg(short = 'o', long, default_value = "sedum")]
    pub output_prefix: String,

    /// Worker threads. 0 = auto (num_cpus). 1 = single-threaded sequential scan.
    #[arg(short = 't', long, default_value_t = 0)]
    pub threads: usize,

    /// Genomic stripes per worker for load balancing (parallel mode only).
    #[arg(long, default_value_t = 16)]
    pub stripes_per_thread: usize,

    /// Force single-threaded scan even if a BAI index is present.
    #[arg(long, default_value_t = false)]
    pub no_parallel: bool,

    /// Disable the live progress line on stderr.
    #[arg(long, default_value_t = false)]
    pub no_progress: bool,
}

impl Cli {
    pub fn effective_threads(&self) -> usize {
        if self.threads == 0 {
            num_cpus::get().max(1)
        } else {
            self.threads
        }
    }
}

fn parse_flag(s: &str) -> Result<u16, String> {
    let s = s.trim();
    let parsed = if let Some(rest) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u16::from_str_radix(rest, 16)
    } else {
        s.parse::<u16>()
    };
    parsed.map_err(|e| format!("invalid flag '{s}': {e}"))
}
