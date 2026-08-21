# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Sedum is a fast, parallel per-bin depth calculator for aligned BAM files, written in Rust. It reads an indexed BAM, partitions the genome into stripes, scans each stripe on a worker thread with a private bin accumulator, and writes per-bin mean depth (`{prefix}.bins.bed`) plus a per-chromosome summary (`{prefix}.summary.tsv`). It is benchmarked against and matches `mosdepth` output bit-for-bit at two-decimal precision, while using a fraction of the memory.

## Commands

```bash
cargo build --release          # optimized binary at target/release/sedum
cargo build                    # debug build
cargo test                     # run unit tests (in bins.rs, plan.rs)
cargo test <name>              # run a single test by name substring, e.g. cargo test spans_three_bins
cargo clippy --all-targets     # lint
cargo fmt                      # format

# Run it
target/release/sedum sample.bam -b 500 --min-mapq 20 -o sample.cov
```

`.cargo/config.toml` sets `-C target-cpu=native` by default so LLVM can auto-vectorize the bin reduction and `add_interval` hot loops and libdeflate picks the widest SIMD path. To produce a **portable** binary, override this:

```bash
RUSTFLAGS="" cargo build --release --no-default-features
```

`--no-default-features` drops the `libdeflate` feature (falls back to miniz_oxide for BGZF inflate). Portable release artifacts are built via cargo-dist (`[workspace.metadata.dist]` in Cargo.toml, `.github/workflows/`).

## Architecture

The pipeline is a fan-out/reduce over genomic stripes. Data flows `main.rs` → `reader.rs` → `stats.rs` → `writer.rs`, with the parallelism strategy decided in `main::run`.

- **`main.rs`** — CLI entry, prints the run banner and final genome summary. Chooses the parallel path only when `!no_parallel && threads > 1 && bai_present()`; otherwise falls back to a sequential scan (also the only path for SAM or index-less BAM).
- **`cli.rs`** — clap arg parsing. SAM flags parse as hex (`0x704`) or decimal via `parse_flag`. Default exclude mask `0x704` = UNMAP|SECONDARY|QCFAIL|DUP, matching mosdepth.
- **`plan.rs`** — `plan()` splits references into `Region` stripes. Target stripe size = `genome_len / (threads * stripes_per_thread)`. Stripes are contiguous, non-overlapping, and cover each reference exactly once.
- **`reader.rs`** — the core. `scan_parallel` spawns N threads (via `std::thread::scope`) that pop `Region`s off a shared `Mutex<Vec<Region>>` work queue (LIFO). **Each worker opens its own `IndexedReader`** (independent BGZF stream) and accumulates into a **private `BinMatrix` shard**. Shards are summed element-wise in a final reduction. `process_record` applies the filter, then walks the CIGAR into bin increments.
- **`bins.rs`** — `BinMatrix` / `RefBins`: `counts[i]` = total aligned reference bases landing in bin `i` (summed depth, not per-base). `add_interval` distributes an interval across overlapping bins, clipping at chromosome end.
- **`cigar.rs`** — `visit_ref_intervals` emits reference-consuming intervals. `M`/`=`/`X`/`D` contribute; `N` (skip) only with `--count-splice`; `I`/`S`/`H`/`P` never do.
- **`filters.rs`** — `Filter::flags_ok(flags)` (`(flags & exclude)==0 && (flags & include)==include`) and `Filter::mapq_ok(mapq)` (`mapq >= min_mapq`), kept as two predicates so `process_record` can attribute each rejection to its own `ScanCounters` bucket. `Copy` so it's cheaply shared across threads.
- **`stats.rs`** — `summarize` computes per-ref and total min/max/mean/median/stdev + breadth at 1/5/10/20/30x. Mean depth uses total aligned bases / reference length; per-bin depth uses each bin's actual span (last bin may be short).
- **`progress.rs`** — atomic counters (Relaxed) batched every `RECORD_BATCH` (4096) records; a `Reporter` thread renders the live stderr line and is dropped before writing to flush.
- **`debug.rs`** — `--debug-stats` support: `ScanCounters` (filter-reason breakdown + bases, always collected, summed across worker shards), `memory_usage()` (peak/current RSS from `/proc/self/status`, Linux-only), `cpu_time()` (cumulative user+system across all threads via `getrusage`, unix-only — `main::run` samples it around the scan to report cores-busy and diagnose CPU- vs IO-bound scans), and `print_report`. `main::run` times each phase (scan/summarize/write) and prints the report only when the flag is set. Note: in parallel mode `records_visited` counts per-stripe query hits, so a read overlapping N stripes is visited N times — the extra visits show up as `out_of_stripe`; `kept` and `aligned_bases` are the true totals.

### Key correctness invariant

Every aligned read is owned by **exactly one stripe** — the one containing its 0-based `alignment_start` (enforced in `process_record` via the `stripe` bounds check). A read whose CIGAR spills past `stripe.end` still contributes fully: the spillover lands in bins that the next stripe's shard also touches, and the final reduction sums all shards. When changing the striping or reduction logic, preserve this: double-counting or dropping reads at stripe boundaries is the main failure mode.

Memory scales with thread count (one `BinMatrix` shard per worker, ~50 MB/worker for human at 500 bp). This is intentional — it's the tradeoff that avoids a per-base depth array.

## Notes

- Licensed **CC BY-NC-ND 4.0** (`license = "CC-BY-NC-ND-4.0"` in `Cargo.toml`, matching README/LICENSE). The LICENSE file is authoritative. `publish = false` — this is a private repo and cargo-dist is configured with no publish jobs.
- Not yet implemented (mentioned in README as planned): CRAM support, `--regions BED`. SAM and index-less BAM force the single-threaded path.
- `PLAN.md` holds the original design doc.
