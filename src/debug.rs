use std::path::Path;
use std::time::Duration;

/// Per-scan counters, always collected (the increments are negligible next to
/// the CIGAR walk) and only rendered when `--debug-stats` is set. In parallel
/// mode each worker keeps its own set and they are summed during reduction.
#[derive(Default, Clone, Copy)]
pub struct ScanCounters {
    /// Records handed to `process_record`. In parallel mode a read that
    /// overlaps N stripes is visited N times, so this exceeds the file's
    /// record count; `out_of_stripe` accounts for the extra visits.
    pub records_visited: u64,
    /// Dropped by the exclude/include SAM flag masks.
    pub filtered_flags: u64,
    /// Passed the flag masks but below `--min-mapq`.
    pub filtered_mapq: u64,
    /// No reference id or no alignment start (unmapped).
    pub unmapped: u64,
    /// Alignment start fell outside the owning stripe (parallel mode only).
    pub out_of_stripe: u64,
    /// Records that contributed depth.
    pub kept: u64,
    /// Reference bases added to bins (sum of emitted interval lengths).
    pub aligned_bases: u64,
}

impl ScanCounters {
    pub fn merge(&mut self, o: &ScanCounters) {
        self.records_visited += o.records_visited;
        self.filtered_flags += o.filtered_flags;
        self.filtered_mapq += o.filtered_mapq;
        self.unmapped += o.unmapped;
        self.out_of_stripe += o.out_of_stripe;
        self.kept += o.kept;
        self.aligned_bases += o.aligned_bases;
    }
}

pub struct MemoryUsage {
    pub peak_rss: u64,
    pub current_rss: u64,
}

/// Peak (VmHWM) and current (VmRSS) resident set size in bytes, read from
/// `/proc/self/status`. Returns `None` on non-Linux or if parsing fails.
pub fn memory_usage() -> Option<MemoryUsage> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let mut peak = None;
    let mut current = None;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            peak = parse_kb(rest);
        } else if let Some(rest) = line.strip_prefix("VmRSS:") {
            current = parse_kb(rest);
        }
    }
    Some(MemoryUsage {
        peak_rss: peak?,
        current_rss: current?,
    })
}

fn parse_kb(s: &str) -> Option<u64> {
    // e.g. "  157000 kB"
    let kb: u64 = s.split_whitespace().next()?.parse().ok()?;
    Some(kb * 1024)
}

/// Cumulative CPU time (user + system) consumed by the whole process,
/// summed across all threads. Diff two samples to attribute CPU to a phase.
/// Returns `None` on unsupported platforms.
#[cfg(unix)]
pub fn cpu_time() -> Option<Duration> {
    use std::mem::MaybeUninit;
    let mut usage = MaybeUninit::<libc::rusage>::uninit();
    // SAFETY: getrusage only writes into the provided rusage on success.
    let ret = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    if ret != 0 {
        return None;
    }
    let usage = unsafe { usage.assume_init() };
    Some(timeval(usage.ru_utime) + timeval(usage.ru_stime))
}

#[cfg(not(unix))]
pub fn cpu_time() -> Option<Duration> {
    None
}

#[cfg(unix)]
fn timeval(t: libc::timeval) -> Duration {
    Duration::new(t.tv_sec.max(0) as u64, (t.tv_usec.max(0) as u32) * 1000)
}

/// Everything the debug report needs, assembled by `main` after the run.
pub struct DebugReport<'a> {
    pub bam_path: &'a Path,
    pub bam_bytes: Option<u64>,
    pub bai_bytes: Option<u64>,
    pub parallel: bool,
    pub threads: usize,
    pub stripes: u64,
    pub bin_size: u32,
    pub counters: ScanCounters,
    pub scan: Duration,
    /// CPU time (user + system, all threads) consumed during the scan phase,
    /// if measurable. Diagnoses CPU- vs IO-bound scans.
    pub scan_cpu: Option<Duration>,
    pub summarize: Duration,
    pub write: Duration,
    pub total: Duration,
}

pub fn print_report(r: &DebugReport) {
    let scan_s = r.scan.as_secs_f64().max(1e-9);
    let c = &r.counters;

    eprintln!();
    eprintln!("==================== debug stats ====================");

    eprintln!("input");
    eprintln!("  path                 {}", r.bam_path.display());
    eprintln!("  bam size             {}", opt_bytes(r.bam_bytes));
    eprintln!("  bai size             {}", opt_bytes(r.bai_bytes));
    if r.parallel {
        eprintln!(
            "  scan mode            parallel ({} threads, {} stripes)",
            r.threads, r.stripes
        );
    } else {
        eprintln!("  scan mode            single-threaded");
    }
    eprintln!("  bin size             {} bp", r.bin_size);

    eprintln!("records");
    if r.parallel {
        eprintln!(
            "  visited              {}   (queries revisit reads overlapping >1 stripe)",
            c.records_visited
        );
    } else {
        eprintln!("  visited              {}", c.records_visited);
    }
    eprintln!("  kept                 {}", c.kept);
    eprintln!("  filtered: flags      {}", c.filtered_flags);
    eprintln!("  filtered: mapq       {}", c.filtered_mapq);
    eprintln!("  unmapped/no-start    {}", c.unmapped);
    if r.parallel {
        eprintln!("  out-of-stripe        {}", c.out_of_stripe);
    }
    eprintln!("  aligned bases        {}", c.aligned_bases);

    eprintln!("timing");
    eprintln!("  scan                 {:.2}s", r.scan.as_secs_f64());
    eprintln!("  summarize            {:.2}s", r.summarize.as_secs_f64());
    eprintln!("  write                {:.2}s", r.write.as_secs_f64());
    eprintln!("  total                {:.2}s", r.total.as_secs_f64());

    eprintln!("throughput (over scan)");
    eprintln!(
        "  records              {:.0} rec/s",
        c.records_visited as f64 / scan_s
    );
    eprintln!(
        "  aligned bases        {}/s",
        human_bytes((c.aligned_bases as f64 / scan_s) as u64)
    );
    if let Some(bytes) = r.bam_bytes {
        eprintln!(
            "  compressed input     {}/s",
            human_bytes((bytes as f64 / scan_s) as u64)
        );
    }

    if let Some(cpu) = r.scan_cpu {
        let cpu_s = cpu.as_secs_f64();
        let cores_busy = cpu_s / scan_s;
        eprintln!("cpu (over scan)");
        eprintln!("  user + system        {cpu_s:.1} core-s");
        if r.parallel && r.threads > 0 {
            eprintln!(
                "  utilization          {cores_busy:.1} cores busy ({:.0}% of {} threads)",
                100.0 * cores_busy / r.threads as f64,
                r.threads
            );
        } else {
            eprintln!("  utilization          {cores_busy:.1} cores busy");
        }
    }

    match memory_usage() {
        Some(m) => {
            eprintln!("memory");
            eprintln!("  peak RSS             {}", human_bytes(m.peak_rss));
            eprintln!("  current RSS          {}", human_bytes(m.current_rss));
        }
        None => eprintln!("memory                 unavailable on this platform"),
    }

    eprintln!("=====================================================");
}

fn opt_bytes(b: Option<u64>) -> String {
    match b {
        Some(n) => human_bytes(n),
        None => "unknown".to_string(),
    }
}

fn human_bytes(n: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    if n < 1024 {
        return format!("{n} B");
    }
    let mut v = n as f64;
    let mut i = 0;
    while v >= 1024.0 && i < UNITS.len() - 1 {
        v /= 1024.0;
        i += 1;
    }
    format!("{v:.2} {}", UNITS[i])
}
