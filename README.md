# Sedum

![Sedum Logo](resources/sedum_logo.png)

A fast, parallel per-bin depth calculator for aligned BAM files, written in Rust.

Sedum reads an indexed BAM, partitions the genome into stripes, scans each stripe on its own worker thread with a private bin accumulator, and writes per-bin mean depth plus a per-chromosome summary (min/max/mean/median/stdev and breadth at 1/5/10/20/30x).

## Why another depth tool

Sedum is built around a simple observation: for the common per-bin output, you do not need to keep a per-base depth array in memory. By accumulating directly into bins and partitioning by BAI stripe, each worker gets a small private shard, and the final reduction is a parallel element-wise sum.

On a 0.75x ONT WGS BAM at 500 bp bins, 2 worker threads:

| | Sedum | mosdepth |
|---|---|---|
| Wall time | 52.0 s | 128.5 s |
| Peak RSS | 157 MB | 1.99 GB |
| CPU% | 153% | 107% |

Sedum's per-bin mean depth matches mosdepth bit-for-bit at mosdepth's two-decimal output precision; total aligned bases per chromosome are exactly equal.

## Install

Requires a Rust toolchain (1.74+).

```bash
git clone <repo> sedum
cd sedum
cargo build --release
# binary at target/release/sedum
```

The build uses libdeflate for BGZF decompression and `-C target-cpu=native` for the host CPU. To produce a portable binary:

```bash
RUSTFLAGS="" cargo build --release --no-default-features
```

## Usage

```
sedum <BAM> [OPTIONS]

  -b, --bin-size <N>             Bin size in bp [default: 500]
  -t, --threads <T>              Worker threads, 0 = num_cpus [default: 0]
      --stripes-per-thread <N>   Stripes per worker for load balancing [default: 16]
      --min-mapq <Q>             Minimum mapping quality [default: 0]
      --exclude-flags <F>        SAM flags to exclude (hex or dec) [default: 0x704]
      --include-flags <F>        SAM flags required [default: 0]
      --count-splice             Count CIGAR N (skipped reference) toward depth
  -o, --output-prefix <P>        Prefix for outputs [default: sedum]
      --no-parallel              Force single-threaded scan
      --no-progress              Disable the live progress line
```

Default exclude mask `0x704` matches mosdepth's defaults: UNMAP | SECONDARY | QCFAIL | DUP.

### Examples

```bash
# Default: 500 bp bins, all cores, BAI auto-detected.
sedum sample.bam

# Larger bins, MAPQ 20 filter, custom output prefix.
sedum sample.bam -b 10000 --min-mapq 20 -o sample.cov

# Force sequential scan (e.g. SAM input, or BAM without BAI).
sedum sample.bam --no-parallel
```

## Output

Two files are written next to your chosen `--output-prefix`:

### `{prefix}.bins.bed`

Plain BED3+1, one row per bin, mean depth as `f64`:

```
chr1    0       500     0.7140
chr1    500     1000    1.0820
chr1    1000    1500    1.2400
...
```

The last bin of each reference may be shorter than `--bin-size`; the reported depth uses that bin's actual span as the denominator, so values stay correct at chromosome ends.

### `{prefix}.summary.tsv`

One row per reference plus a final `total` row:

```
chrom   length  n_bins  total_aligned_bases  mean    min   max     median  stdev   breadth_ge1  breadth_ge5  breadth_ge10  breadth_ge20  breadth_ge30
chr1    248956422  497913  267832462         1.0758  0.0   593.3   1.0     5.95    0.5699       0.0101       0.0012        0.0005        0.0004
...
total   3099922541 6199940 3087022690        0.9958  0.0   1565.0  1.0     3.63    0.5690       0.0079       0.0009        0.0004        0.0003
```

Breadth columns are the fraction of bins with mean depth >= the threshold.

## Design

```
[ BAM + BAI ]
      |
      |  BAI-driven stripe planner
      |  target stripe size = genome_len / (threads x stripes_per_thread)
      v
[ work queue of Region(ref_id, start, end) ]
      |
      |  shared Mutex<Vec<Region>>, workers pop from the back
      v
[ N worker threads ]
   - each opens its own IndexedReader (independent BGZF stream)
   - reader.query(region) -> records that overlap the stripe
   - count records whose alignment_start falls inside [start, end)
   - CIGAR walk -> add_interval into a private BinMatrix shard
      |
      v
[ parallel reduction: sum N shards element-wise -> final BinMatrix ]
      v
[ writer: per-bin BED + summary TSV ]
```

Key correctness invariant: every aligned read is owned by exactly one stripe (the one containing its `alignment_start`). Reads whose CIGAR spills past `stripe.end` still contribute fully because the spillover lands in bins also touched by the next stripe's shard, and the final reduction is a sum.

Performance choices:

- BGZF inflate uses libdeflate (SSE2/AVX2 internally).
- Per-worker shard is `Vec<u64>` per reference; for human at 500 bp this is ~50 MB per worker, dominating memory.
- The hot path increments a u64 with a branch-light loop; the compiler auto-vectorizes the final reduction under `-C target-cpu=native`.
- Progress counters use Relaxed atomics with 4096-record batching so the reporter thread costs nothing measurable.

## Filters

A record contributes to depth iff:

- `(flags & exclude_flags) == 0`
- `(flags & include_flags) == include_flags`
- `mapq >= min_mapq`
- It is mapped (`reference_sequence_id` and `alignment_start` both present)

CIGAR ops that consume reference and contribute to depth: `M`, `=`, `X`, `D`. CIGAR `N` (skipped reference) is excluded by default; pass `--count-splice` for RNA-seq.

## Limitations

- BAM only at present. CRAM and SAM-without-index force the sequential path.
- No region restriction flag yet (`--regions BED` is planned).
- One BinMatrix shard per worker means memory scales with thread count. At 32 threads and 500 bp bins on a human genome, expect ~1.5 GB peak.

## License

MIT.
