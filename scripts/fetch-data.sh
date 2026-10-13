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
# Orange coat in cats: 9 kb of chrX (Felis_catus_9.0) around ARHGAP36, with
# its gene annotation and the 5,074 bp deletion that makes cats orange
# (NC_018741.3:g.110432079_110437152del; Toh et al. and Kaelin et al., 2025).
if [ ! -s data/orange-cat.fa ]; then
    curl -sS --fail --retry 3 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_018741.3&rettype=fasta&retmode=text&seq_start=110430079&seq_stop=110439152" \
        -o data/orange-cat.fa
    sleep 1
fi
if [ ! -s data/orange-cat.gff3 ]; then
    # NCBI returns the whole chromosome's annotation regardless of range, so
    # keep only what overlaps the 120 kb around the gene.
    curl -sS --fail --retry 3 \
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NC_018741.3" |
        awk -F'\t' '/^##gff-version/ || (NF >= 9 && $3 != "region" && $3 != "cDNA_match" && $5 >= 110380000 && $4 <= 110500000)' \
        > data/orange-cat.gff3
fi
printf 'NC_018741.3\t110432078\t110437152\torange\n' > data/orange-cat.bed
# The cat behind the F.catus_Fca126_mat1.0 assembly is orange, so its sequence
# over the same stretch is an orange allele. Tiled into reads and mapped back,
# it gives the coverage an orange male's sequencing would: none over the deletion.
if [ ! -s data/orange-cat-fca126.fa ]; then
    sleep 1
    curl -sS --fail --retry 3 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_058386.1&rettype=fasta&retmode=text&seq_start=107716000&seq_stop=107728000" \
        -o data/orange-cat-fca126.fa
fi
# The same stretch in other species, for the comparison lanes. Each is 14 kb
# ending 39 kb before that species' ARHGAP36 3' end, which is where the cat
# window lines up (annotated start sites vary too much between species to
# anchor on). The page aligns them itself.
while read -r sp acc end; do
    out="data/orange-cat-$sp.fa"
    [ -s "$out" ] && continue
    sleep 1
    curl -sS --fail --retry 3 \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text&seq_start=$((end - 53000))&seq_stop=$((end - 39000))" \
        -o "$out"
done <<'SPECIES'
tiger NC_056677.1 107761871
dog NC_051843.1 104084083
human NC_000023.11 131089884
mouse NC_000086.8 48589120
SPECIES
python3 -I scripts/simulate-coverage.py data/orange-cat.fa data/orange-cat-fca126.fa "orange (Fca126)" \
    > data/orange-cat.bedgraph
printf '  %-14s %s\n' orange-cat.fa "cat chrX, 9 kb — ARHGAP36, with .gff3 genes and .bed deletion"
# The web demo serves its samples from site/data.
mkdir -p site/data
cp data/*.fa data/*.gff3 data/*.bed data/*.bedgraph site/data/ 2>/dev/null || true

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
