# Sedum — Implementation Plan

A fast, parallel per-bin depth calculator for aligned BAM/SAM (and CRAM-ready) files, written in Rust.

## 1. Crate layout

```
sedum/
├── Cargo.toml
├── src/
│   ├── main.rs           # CLI entry, wires modules
│   ├── cli.rs            # clap-derive argument parsing
│   ├── plan.rs           # Region planning (genome partitioning from BAI)
│   ├── reader.rs         # BAM/SAM/CRAM reader abstraction + parallel BGZF
│   ├── filters.rs        # MAPQ, flag, dup/secondary/supplementary filters
│   ├── cigar.rs          # CIGAR → reference span iteration
│   ├── bins.rs           # Per-chromosome bin accumulators
│   ├── aggregate.rs      # Merge per-worker shards
│   ├── stats.rs          # Summary statistics (min/max/mean/median/stdev/Nx)
│   └── writer.rs         # BED/TSV per-bin + summary writers (bgzipped optional)
└── tests/                # Integration tests on toy BAMs
```

## 2. Dependencies (pure-Rust first, htslib opt-in)

- `noodles` ecosystem — `noodles-bam`, `noodles-sam`, `noodles-bgzf`, `noodles-csi`/`noodles-bai`, `noodles-cram` (optional feature)
- `crossbeam-channel` — bounded MPMC channels for the producer/consumer pipeline
- `rayon` — data-parallel reductions for aggregation/summary
- `clap` (derive) — CLI
- `anyhow` + `thiserror` — error handling
- `indicatif` — progress bars (optional `--progress` flag)
- `flate2` / `bgzip` writer for optional gzip output
- Cargo feature `htslib` to swap in `rust-htslib` for users wanting the C decoder

## 3. CLI surface (clap)

```
sedum <BAM> [--regions BED] [--bin-size N] [--min-mapq Q]
            [--exclude-flags 0x704] [--include-flags 0]
            [--threads T] [--decompress-threads D]
            [--output-prefix PREFIX]
            [--no-per-base] [--fast-mode] [--fasta REF (for CRAM)]
            [--format {bed,tsv,bedgraph}]
```

Sensible defaults:
- `--bin-size`: 500 bp (mosdepth uses none by default; we default to 500 because the user wants per-bin output)
- `--threads`: `num_cpus::get()`
- `--decompress-threads`: `max(1, threads/2)` — see §5
- `--exclude-flags`: `0x704` (UNMAP|SECONDARY|QCFAIL|DUP) — samtools-flagstat-friendly default
- `--min-mapq`: 0

## 4. Region planning (`plan.rs`)

1. Load BAI (or CSI). Fail-fast if missing and SAM input cannot be partitioned.
2. Build `Vec<Region>` where each region is `(ref_id, start, end, virtual_offset_chunks)`.
3. Strategy:
   - For each reference: split into stripes of `≈ chrom_len / (threads × oversubscription)` bp, target ~16 stripes per thread for load balancing.
   - Use BAI's 16 kb linear index to translate each stripe into BGZF virtual-offset ranges → workers seek directly, no scanning.
4. Optional `--regions BED` intersects with the partition before scheduling.

## 5. Threading model — the decompression/worker balance

The bottleneck on modern NVMe is BGZF inflate, not IO. Two-stage pipeline:

```
[BAM file]
    │
    ├── D × BGZF decompression threads (noodles-bgzf::MultithreadedReader OR
    │      one indexed reader per region with internal mt inflate)
    │
    ├── bounded channel of decoded record batches (Vec<Record>, ~4k records)
    │
    └── W × worker threads → per-thread bin shard (chrom → Vec<u32>)
            │
            └── on completion: rayon reduce into the global accumulator
```

Key choices:

- **Two pools, not one.** Decompression is CPU-bound on inflate; workers are CPU-bound on CIGAR walking. Putting both in rayon's global pool causes head-of-line blocking when a long CIGAR stalls a decompression slot.
- **Per-region readers** (preferred over a single shared multithreaded reader) because each worker can seek to its BAI chunk and stream linearly. This sidesteps the lock on a shared reader and gives near-linear scaling up to disk bandwidth.
- **Decomp threads default to threads/2**, capped by measured benefit: noodles-bgzf scales well to ~4 inflate threads per stream. Expose `--decompress-threads` for tuning.
- **Bounded channel (capacity ~ 2 × W)** applies natural backpressure so a slow worker doesn't balloon memory.
- **Batch records** (4096 per message) to amortize channel + allocation cost.

## 6. Per-bin accumulation (`bins.rs`)

Per worker, per chromosome: `Vec<u32>` length `ceil(chrom_len / bin_size)`.

For each record passing filters:
1. Walk CIGAR, emitting `(ref_start, ref_end)` aligned intervals (M, =, X consume ref+query; D, N consume ref only — N skipped from depth by default, configurable with `--count-splice`).
2. For each interval, compute first/last bin and:
   - Full-overlap bins: `bin[i] += bin_size`
   - Partial: `bin[i] += overlap_len`
3. Store **summed aligned bases per bin**; convert to mean depth at write time (`sum / bin_size`).

Optimizations:
- Skip records with `ref_id < 0` early.
- Use `u32` per bin (4 GB max sum per bin — safe for bin_size ≤ 1 Mb and reasonable depths). Promote to `u64` via feature flag for ultra-deep targeted sequencing.
- SIMD-free; rely on the compiler. Hot loop should be branch-predictable.

## 7. Aggregation (`aggregate.rs`)

After channel close:
- Reduce N worker shards into one global `BinMatrix` using rayon's `par_iter().reduce()`, summing element-wise per chromosome.
- Memory ceiling: `(n_bins × 4) × W`. For human + 500 bp + 16 workers ≈ 200 MB — acceptable.

## 8. Summary statistics (`stats.rs`)

Per chromosome and genome-wide:
- `n_bins`, `total_aligned_bases`, `mean_depth`, `min`, `max`, `median`, `stdev`
- Coverage breadth: fraction of bins with depth ≥ {1, 5, 10, 20, 30}
- Optional Nx-style metrics (`--nx 50,90`)

Compute via streaming over the global bin matrix; median via in-place `select_nth_unstable` on a clone per chromosome (parallel with rayon).

## 9. Output (`writer.rs`)

- `{prefix}.bins.bed[.gz]` — `chrom\tstart\tend\tmean_depth` (BED3+1; alt `bedgraph`)
- `{prefix}.summary.tsv` — per-chrom stats + final `total` row
- Optional `--per-base` writes a separate `{prefix}.per-base.bed.gz` (mosdepth-parity); off by default since the user asked for bins.
- Use `BufWriter` over a bgzip-capable writer when output ends in `.gz`.
- Write summary last so it can include final aggregates without buffering bins in memory.

## 10. Correctness & tests

- Golden-file tests against `samtools depth -a` and `mosdepth` on a small NA12878 chr22 BAM (vendored or downloaded in CI).
- Property test: random synthetic BAM via `noodles-sam` writer → recompute depth in pure Rust naïvely → must match.
- Edge cases: empty BAM, unmapped-only, CIGAR with N/S/H, bin boundaries, chromosomes shorter than one bin, multiple references with same name guard.

## 11. Benchmarks

- `criterion` benches on a fixed BAM measuring throughput vs `--threads ∈ {1,2,4,8,16}` and `--decompress-threads ∈ {1,2,4,8}`.
- Track wall-clock against `mosdepth` and `samtools depth -a -@ N` on the same hardware; record in `BENCHMARKS.md`.

## 12. Milestones

1. Skeleton crate, CLI, single-threaded reader → correct bin output vs mosdepth on chr22.
2. Add filters, summary writer, golden tests.
3. Introduce per-region BAI planner; switch to worker pool with shared decoders.
4. Split decompression pool from worker pool; tune channel/batch sizes via criterion.
5. CRAM support behind `cram` feature; htslib backend behind `htslib` feature.
6. Documentation, README, packaging (cargo install + GitHub Actions release binaries).
