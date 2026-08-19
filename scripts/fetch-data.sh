#!/usr/bin/env bash
# Fetch reference genomes used for manual testing and integration tests.
# Public sequences from NCBI Nucleotide. Total ~4.8 MB.
set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p data

fetch() {
    local acc="$1" out="data/$2" label="$3"
    if [ -s "$out" ]; then
        printf '  %-14s %s (cached)\n' "$2" "$label"
        return
    fi
    curl -sS --fail --retry 3 --retry-delay 2 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" \
        -o "$out"
    # NCBI answers rate limits with HTTP 200 and a JSON body, so check content.
    if head -c 1 "$out" | grep -q '{'; then
        rm -f "$out"
        echo "  $2: rate limited by NCBI — wait a few seconds and re-run" >&2
        exit 1
    fi
    printf '  %-14s %s\n' "$2" "$label"
    # Stay under NCBI's 3 requests/second limit for unauthenticated clients.
    sleep 1
}

echo "fetching reference genomes into data/"
fetch NC_001422.1 phix174.fa "phiX174, 5.4 kb — Illumina spike-in control"
fetch NC_001416.1 lambda.fa  "Lambda phage, 48.5 kb — nanopore calibration strand"
fetch NC_045512.2 sars2.fa   "SARS-CoV-2, 29.9 kb"
fetch NC_000913.3 ecoli.fa   "E. coli K-12 MG1655, 4.6 Mb"
echo "done"
