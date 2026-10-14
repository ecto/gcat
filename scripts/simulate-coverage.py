#!/usr/bin/env python3
"""Simulated read depth of one genome over another, as bedGraph.

    simulate-coverage.py REFERENCE.fa SAMPLE.fa NAME > out.bedgraph

Tiles SAMPLE into error-free 150 bp reads every few bases and places each on
REFERENCE by exact 31-mer matches, counting a base as covered when a k-mer
that spans it matches uniquely. A read that crosses a deletion in SAMPLE
lands as two partial alignments, the way a real aligner would soft-clip it,
so depth ramps down into each breakpoint and is zero across the gap.

The reference header's range suffix (NC_018741.3:110430079-110439152) puts
the output in chromosome coordinates, matching the annotation files.
"""
import sys

READ, STEP, K = 150, 4, 31


def read_fasta(path):
    with open(path) as f:
        header = f.readline()[1:].split()[0]
        return header, "".join(l.strip() for l in f).upper()


def main(ref_path, sample_path, name):
    ref_id, ref = read_fasta(ref_path)
    _, sample = read_fasta(sample_path)
    seqid, _, rng = ref_id.partition(":")
    offset = int(rng.split("-")[0]) - 1 if rng else 0

    seen = {}
    for i in range(len(ref) - K + 1):
        kmer = ref[i : i + K]
        seen[kmer] = -1 if kmer in seen else i  # -1 marks a repeat
    depth = [0] * len(ref)
    for r in range(0, len(sample) - READ + 1, STEP):
        covered = set()
        for j in range(r, r + READ - K + 1):
            at = seen.get(sample[j : j + K], -1)
            if at >= 0:
                covered.update(range(at, at + K))
        for p in covered:
            depth[p] += 1

    print(f'track type=bedGraph name="{name}"')
    start = 0
    for p in range(1, len(ref) + 1):
        if p == len(ref) or depth[p] != depth[start]:
            print(f"{seqid}\t{offset + start}\t{offset + p}\t{depth[start]}")
            start = p


if __name__ == "__main__":
    main(*sys.argv[1:4])
