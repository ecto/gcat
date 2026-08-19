# gcat

`cat` for genomes.

```bash
gcat stats reads.fq
gcat grep -p ACGTACGT reads.fq | gcat stats
gcat kmer -k 21 --top 20 genome.fa
```

gcat is a genomics compute platform built on [tang](https://github.com/ecto/tang).
Small sharp tools over streams, on a substrate that runs the same everywhere —
a laptop, a datacenter GPU, a browser, or an embedded box next to the sequencer.

## Why

Genomics compute is locked to specific hardware. The state-of-the-art nanopore
basecaller depends on a closed-source, NVIDIA-only kernel library, which means
practical high-speed basecalling is available on a narrow class of machines and
nowhere else. Meanwhile the numerical core of the field — REML over relatedness
matrices, SVD over genotype matrices, sequence-model inference — is spread
across tools that each rebuild the same primitives.

tang already solves the portability half: one expression graph, compiled to
Metal, Vulkan, DX12, or CPU, `no_std` throughout. gcat is that substrate pointed
at genomics.

## Crates

| Crate | What it does |
|-------|--------------|
| **`gcat`** | Core types — `Base`, 2-bit packed `Seq`, rolling `Kmer`. `no_std`. |
| **`gcat-io`** | Streaming FASTA/FASTQ readers and writers. |
| **`gcat-signal`** | Raw nanopore signal — med/MAD scaling, overlapped windowing, tensor handoff. `no_std`. |
| **`gcat-cli`** | The `gcat` binary. |

`gcat` and `gcat-signal` build for `wasm32-unknown-unknown`, which is the
standing check that the portability claim is real and not aspirational.

## Status

Stage 1: substrate. The core types, I/O, and the signal front-end are here and
tested. The compute path is wired far enough to prove it holds — `gcat-signal`
emits tang tensors under `--features tang` — but no model runs yet.

Stage 2 is the proof: a portable basecaller. That is the sharpest unmet need in
the field and the thing this stack is unusually suited to, since it is a
convolutional/recurrent network over a 1-D signal, which tang already has the
layers for.

The blocker there was CTC — the loss that trains a sequence model when nobody
has told you which input sample produced which base. tang now has it
(`ctc_loss`, `ctc_loss_grad`, `greedy_decode`, `beam_decode`), and
`crates/gcat-signal/tests/basecall.rs` runs the whole seam: synthesize a trace,
scale it, build a tensor, take a CTC gradient, decode back to DNA.

That test needs `--features tang` **and** a tang checkout carrying CTC. Until
that lands on tang's default branch, the default build skips it.

## Test data

```bash
./scripts/fetch-data.sh
```

Pulls four reference genomes from NCBI (~4.8 MB): phiX174, lambda phage,
SARS-CoV-2, and E. coli K-12. Lambda is the nanopore calibration strand and
phiX174 is the Illumina spike-in control, so both are what a sequencer emits
when you ask it to prove it works.

Integration tests assert gcat's output against the published length and GC
content for each accession. They skip when `data/` is absent, so a fresh clone
still runs green.

| | length | GC |
|---|---|---|
| phiX174 | 5,386 | 44.76% |
| lambda | 48,502 | 49.86% |
| SARS-CoV-2 | 29,903 | 37.97% |
| E. coli K-12 | 4,641,652 | 50.79% |

## Install

```bash
cargo install --path crates/gcat-cli
```

tang is unpublished, so `gcat-signal --features tang` expects a checkout at
`../tang`.

## Commands

| | |
|---|---|
| `gcat cat` | Concatenate records (the default — `gcat reads.fq` works) |
| `gcat stats` | Count, bases, N50, GC, mean Q |
| `gcat head -n N` | First N records |
| `gcat grep -p PAT` | Records matching a pattern in id or sequence |
| `gcat rc` | Reverse complement, quality reversed with it |
| `gcat kmer -k K` | Canonical k-mer counts |

Reads stdin when given no files or `-`, writes stdout, and exits quietly on a
closed pipe. FASTA and FASTQ are detected from the first record.

## Design

**Ambiguity is data, not a parse error.** `N` and friends get their own error
variant, distinct from genuine garbage, and `from_ascii_lossy` reports how many
positions it substituted so callers can decide whether a read is worth keeping.

**Quality averages in probability space.** Phred scores are logarithmic;
averaging them directly overstates read quality. `mean_quality` converts,
averages, and converts back.

**Med/MAD, not mean/stddev.** Nanopore traces are full of real outliers —
stalls, blockages, adapter. The median does not chase them. Both are
implemented; the test suite asserts the difference.

**A file that will not open is an error.** Reading it as an empty stream would
make `gcat stats typo.fa` report zero records and exit 0, which is worse than
useless.

**Degenerate inputs return identity, not NaN.** A constant trace has zero
deviation. Dividing by it would poison every tensor downstream, so scaling falls
back to 1.0.

## License

MIT
