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
# The web demo serves its samples from site/data.
mkdir -p site/data
cp data/*.fa site/data/ 2>/dev/null || true

# A synthetic FASTQ so the read-quality panel has something to show. Long reads
# are deliberately worse, which is the correlation the length-vs-quality plot
# exists to reveal.
python3 - <<'PYGEN'
import random, os
random.seed(7)
out = 'site/data/sample.fq'
if os.path.exists(out) and os.path.getsize(out) > 0:
    raise SystemExit
bases = 'ACGT'
with open(out, 'w') as f:
    for i in range(4000):
        length = int(200 * (1.6 ** random.uniform(0, 8)))
        centre = max(4, min(28, int(26 - length / 3000)))
        seq = ''.join(random.choice(bases) for _ in range(length))
        qual = ''.join(chr(33 + max(2, int(random.gauss(centre, 2.5)))) for _ in range(length))
        f.write(f"@read{i} len={length}\n{seq}\n+\n{qual}\n")
PYGEN

echo "done"
