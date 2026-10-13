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
| **`gcat-io`** | Streaming FASTA/FASTQ readers and writers; GFF3, BED and bedGraph annotation. |
| **`gcat-signal`** | Raw nanopore signal — med/MAD scaling, overlapped windowing, tensor handoff. `no_std`. |
| **`gcat-viz`** | Visualization kernels — chaos game, GC skew, 2-D histograms. `no_std`. |
| **`gcat-align`** | Seed-and-chain pairwise alignment, for comparing a region across species. `no_std`. |
| **`gcat-wasm`** | Browser bindings. Streaming analysis, zero-copy pixel buffers. |
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

## The browser demo

```bash
./scripts/fetch-data.sh
wasm-pack build crates/gcat-wasm --target web --out-dir ../../site/pkg --release
python3 -m http.server 8731 --directory site
```

Pure black, with the palette from campedersen.com — gold as the interface
accent, and bases on the Sanger convention drawn from the same ramp: A lime,
C sky, G orange, T red. Gold is deliberately not a base colour, because a
letter sharing its colour with the selection highlight makes both harder to
read.

An app, not a page: title bar, view rail, status bar, keyboard shortcuts, and
panes that fill the viewport. `1`–`5` switch views, `/` finds a motif, `g` jumps
to a base, `esc` clears. Palette and last view persist.

Drop a FASTA or FASTQ anywhere in the window. It is read with the File API and
processed in WebAssembly — there is no server to send it to. E. coli's 4.6 Mb
streams through in about 120 ms.

**A hex editor for DNA.** Offsets down the left, bases in the middle coloured
by nucleotide, GC and Phred quality in the gutter. Rows are virtualized, so a
4.6 Mb genome is 77,000 rows of which only the ~28 on screen exist in the DOM.
Drag to select and get length, GC, and reverse complement; type a motif to
search it, `N` acting as a wildcard. The minimap along the bottom shows GC
deviation across the whole sequence — click it to jump.

Searching lambda for `GAATTC` returns 5 hits, which is exactly how many EcoRI
sites lambda phage has.

Three more views, all computed in the same pass, and a dot plot:

**Chaos game k-mer map.** Each base pulls a point halfway toward its corner of
a square. The usual description makes this sound like a random walk, but after
`k` jumps the position depends only on the last `k` bases — so a 2^k square is
exactly the histogram of every k-mer, laid out so that k-mers sharing a suffix
sit near each other. Zooming in is looking at longer k-mers. E. coli fills its
9-mer grid with visible structure; random sequence renders as flat noise.

**Cumulative GC skew.** Where the curve turns marks where DNA replication
starts. On E. coli it puts the origin at 3,920,000 bp against a published oriC
of 3,925,744 — within one window, from counting G and C. No alignment, no
reference, no database.

**Read length against quality.** A 2-D histogram rather than a scatter, so it
holds up at ten million reads. Answers the question neither axis answers alone:
are the long reads also the good ones?

**Dot plot.** Every 12-mer the open sequence shares with itself, or with a
species opened after it, binned to a square: a diagonal where they agree, a gap
or step where one lost sequence, parallel lines for repeats, red lines against
the grain for inversions. Hover a pixel for both positions. Cat against itself
shows no repeat linking the two ends of the orange deletion, so it was a clean
break rather than a slip between repeats.

**Annotation and comparison.** Open a GFF3, BED or bedGraph with a sequence,
or after it, and a track appears above the editor: genes to scale with exons
tall where they code for protein, BED regions bracketed with an allele lane
that splices each one out, and bedGraph read depth beneath a zoom inset. Open
a second FASTA and it is aligned onto the first in the browser, drawn as a
percent-identity lane. A header naming a range (`NC_018741.3:110430079-…`, as
NCBI and faidx write it) puts everything in chromosome coordinates, including
the editor's offsets. Hover anything on the track for what it is.

**Why orange cats are orange.** The orange-cat sample is 9 kb of cat chrX
around *ARHGAP36*, carrying the 5,074 bp deletion that two groups found in 2025
switches the gene on in pigment cells. The deletion is struck through and
bracketed. Below it, coverage from tiling the orange cat behind the Fca126
assembly into reads drops to zero over exactly the published coordinates.
Tiger, dog, human and mouse lanes show most of the deleted stretch is shared
across mammals. A calico panel beside the track simulates the X-inactivation
that turns one allele into patches.

The wasm bundle is 121 KB, 58 KB gzipped.

## Test data

```bash
./scripts/fetch-data.sh
```

Pulls four reference genomes from NCBI (~4.8 MB): phiX174, lambda phage,
SARS-CoV-2, and E. coli K-12. Plus the orange-cat sample: the cat chrX window,
its annotation, the orange Fca126 sequence it simulates coverage from, and the
matching stretch in four other species (~175 KB). Lambda is the nanopore calibration strand and
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

**Pixels cross the boundary as a pointer.** The kernels write RGBA straight
into WebAssembly memory and JavaScript wraps that buffer in an `ImageData`, so a
megapixel image costs one pass and no serialisation.

**Kernels are histograms on purpose.** Bounded integer accumulators updated by a
local rule — the shape that maps to a GPU storage buffer and an atomic
increment without a rewrite, for when this moves to `tang-gpu`.

**A file that will not open is an error.** Reading it as an empty stream would
make `gcat stats typo.fa` report zero records and exit 0, which is worse than
useless.

**Alignment fills only short gaps.** Between seeds the aligner runs a global
alignment, which lines up about half of any two unrelated stretches by chance.
Filling gaps past 500 bases invents conservation, so longer ones stay blank.

**Degenerate inputs return identity, not NaN.** A constant trace has zero
deviation. Dividing by it would poison every tensor downstream, so scaling falls
back to 1.0.

## License

MIT
