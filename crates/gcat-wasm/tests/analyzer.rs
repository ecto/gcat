//! The analyzer runs natively too, so its stream handling is testable without
//! a browser. Chunk boundaries are the interesting part: a byte stream splits
//! wherever the network or the File API feels like it.

use gcat_wasm::{Analyzer, Annotations};

fn feed(data: &[u8], chunk: usize) -> Analyzer {
    let mut a = Analyzer::new(4, 10);
    for part in data.chunks(chunk.max(1)) {
        a.push(part);
    }
    a.finish();
    a
}

const FASTA: &[u8] = b">chr1 test\nACGTACGTAC\nGGCCTTAANN\n>chr2\nTTTTGGGGCC\n";
const FASTQ: &[u8] = b"@r1\nACGTACGTGG\n+\nIIIIIIIIII\n@r2\nGGGGCCCCAA\n+\nIIIIIIII##\n";

#[test]
fn parses_fasta_regardless_of_chunk_boundaries() {
    // Every chunk size must give the same answer, including one byte at a time.
    let reference = feed(FASTA, FASTA.len());
    for chunk in [1usize, 2, 3, 7, 13, 31] {
        let a = feed(FASTA, chunk);
        assert_eq!(a.records(), reference.records(), "records at chunk {chunk}");
        assert_eq!(a.bases(), reference.bases(), "bases at chunk {chunk}");
        assert_eq!(a.seq_len(), reference.seq_len(), "seq_len at chunk {chunk}");
        assert_eq!(
            a.gc_fraction(),
            reference.gc_fraction(),
            "gc at chunk {chunk}"
        );
    }
}

#[test]
fn parses_fastq_regardless_of_chunk_boundaries() {
    let reference = feed(FASTQ, FASTQ.len());
    for chunk in [1usize, 5, 11, 23] {
        let a = feed(FASTQ, chunk);
        assert_eq!(a.records(), reference.records(), "records at chunk {chunk}");
        assert_eq!(a.bases(), reference.bases(), "bases at chunk {chunk}");
        assert!(a.is_fastq());
        assert!(a.has_quality());
    }
}

#[test]
fn counts_records_and_bases() {
    let a = feed(FASTA, 8);
    assert_eq!(a.records(), 2.0);
    assert_eq!(a.bases(), 30.0);
    // Two Ns are ambiguous and excluded from the GC denominator.
    assert_eq!(a.ambiguous(), 2.0);
    assert_eq!(a.seq_len(), 30.0);
}

#[test]
fn retains_sequence_verbatim_including_ambiguity() {
    let mut a = feed(FASTA, 8);
    let n = a.read_bases(0, 30);
    assert_eq!(n, 30);
    // Reaching into the scratch buffer the way the editor does.
    let ptr = a.scratch_ptr();
    let got: Vec<u8> = (0..n).map(|i| unsafe { *ptr.add(i) }).collect();
    assert_eq!(&got[..20], b"ACGTACGTACGGCCTTAANN");
    assert_eq!(
        &got[18..20],
        b"NN",
        "ambiguity must survive to the editor, not be swapped for a real base"
    );
}

#[test]
fn reading_past_the_end_is_short_not_fatal() {
    let mut a = feed(FASTA, 8);
    assert_eq!(a.read_bases(25, 100), 5);
    assert_eq!(a.read_bases(30, 10), 0);
    assert_eq!(a.read_bases(9999, 10), 0);
}

#[test]
fn gc_of_a_span_ignores_ambiguity() {
    let a = feed(FASTA, 8);
    // "GGCC" is entirely GC.
    assert_eq!(a.gc_of(10, 4), 1.0);
    // "TTAA" has none.
    assert_eq!(a.gc_of(14, 4), 0.0);
    // A span of only Ns has no defined GC.
    assert_eq!(a.gc_of(18, 2), -1.0);
}

#[test]
fn gc_profile_covers_the_sequence() {
    let a = feed(FASTA, 8);
    let p = a.gc_profile(5);
    assert_eq!(p.len(), 5);
    assert!(p.iter().all(|&v| v == -1.0 || (0.0..=1.0).contains(&v)));
}

#[test]
fn finds_motifs_and_wildcards() {
    let mut a = feed(FASTA, 8);
    assert_eq!(a.find_motif("ACGT", 0), 0.0);
    assert_eq!(a.find_motif("ACGT", 1), 4.0);
    assert_eq!(a.find_motif("GGCC", 0), 10.0);
    assert_eq!(a.find_motif("ZZZZ", 0), -1.0);
    // N is a wildcard, so this matches ACGT and ACGG and so on.
    assert_eq!(a.find_motif("ACGN", 0), 0.0);
    assert_eq!(a.count_motif("ACGT"), 2.0);
    // Overlaps count: one TT in chr1's TTAA, and the TTTT in chr2 yields three.
    assert_eq!(a.count_motif("TT"), 4.0);
}

#[test]
fn motif_search_is_case_insensitive() {
    let mut a = feed(FASTA, 8);
    assert_eq!(a.find_motif("acgt", 0), 0.0);
}

#[test]
fn motif_search_handles_degenerate_input() {
    let mut a = feed(FASTA, 8);
    assert_eq!(a.find_motif("", 0), -1.0);
    assert_eq!(a.find_motif(&"A".repeat(500), 0), -1.0);
    // Searching from past the end must not panic.
    assert_eq!(a.find_motif("ACGT", 9999), -1.0);
}

#[test]
fn revcomp_of_a_span() {
    let a = feed(FASTA, 8);
    assert_eq!(a.revcomp_of(0, 4), "ACGT");
    assert_eq!(a.revcomp_of(10, 4), "GGCC");
    assert_eq!(a.revcomp_of(14, 4), "TTAA");
    // Ambiguity comes back as N rather than a fabricated base.
    assert_eq!(a.revcomp_of(18, 2), "NN");
    assert_eq!(a.revcomp_of(9999, 4), "");
}

#[test]
fn quality_is_retained_for_fastq() {
    let mut a = feed(FASTQ, 7);
    assert!(a.has_quality());
    let n = a.read_quality(0, 10);
    assert_eq!(n, 10);
    let ptr = a.scratch_ptr();
    let got: Vec<u8> = (0..n).map(|i| unsafe { *ptr.add(i) }).collect();
    assert_eq!(&got, b"IIIIIIIIII");
}

#[test]
fn fasta_has_no_quality() {
    let a = feed(FASTA, 8);
    assert!(!a.has_quality());
    assert!(!a.is_fastq());
    assert_eq!(a.mean_quality(), -1.0);
}

#[test]
fn length_stats_are_per_record() {
    let a = feed(FASTA, 8);
    let s = a.length_stats();
    assert_eq!(s[0], 2.0, "two records");
    assert_eq!(s[1], 10.0, "shortest is chr2");
    assert_eq!(s[4], 20.0, "longest is chr1");
}

#[test]
fn empty_input_does_not_panic() {
    let mut a = Analyzer::new(4, 10);
    a.finish();
    assert_eq!(a.records(), 0.0);
    assert_eq!(a.seq_len(), 0.0);
    assert_eq!(a.gc_fraction(), 0.0);
    assert_eq!(a.find_motif("ACGT", 0), -1.0);
    assert_eq!(a.length_stats(), vec![0.0; 6]);
    assert_eq!(a.gc_profile(4).len(), 0);
}

#[test]
fn bare_sequence_without_a_header_still_parses() {
    let a = feed(b"ACGTACGTAC\nACGTACGTAC\n", 6);
    assert_eq!(a.bases(), 20.0);
    assert_eq!(a.seq_len(), 20.0);
}

#[test]
fn captures_the_first_header_only() {
    let a = feed(FASTA, 5);
    assert_eq!(a.header(), "chr1 test");
    let q = feed(FASTQ, 5);
    assert_eq!(q.header(), "r1");
    // No header at all is not an error.
    let bare = feed(b"ACGTACGT\n", 3);
    assert_eq!(bare.header(), "");
}

#[test]
fn header_range_sets_the_coordinate_offset() {
    let mut a = Analyzer::new(4, 100);
    a.push(b">NC_018741.3:110430079-110439152 Felis catus\nACGT\n");
    a.finish();
    assert_eq!(a.seq_name(), "NC_018741.3");
    assert_eq!(a.seq_offset(), 110430078.0);
}

#[test]
fn annotations_query_by_window() {
    let mut ann = Annotations::new();
    let gff = "chrX\t.\tgene\t101\t500\t.\t-\t.\tID=g;Name=G\n\
               chrX\t.\tmRNA\t101\t500\t.\t-\t.\tID=t;Parent=g\n\
               chrX\t.\texon\t101\t200\t.\t-\t.\tParent=t\n\
               chrX\t.\texon\t401\t500\t.\t-\t.\tParent=t\n\
               chrX\t.\tCDS\t150\t450\t.\t-\t0\tParent=t\n";
    assert_eq!(ann.add_gff3(gff).unwrap(), 1);
    assert_eq!(ann.add_bed("chrX 250 300 del\nchr1 0 10\n").unwrap(), 2);

    assert_eq!(
        ann.genes_in("chrX", 0.0, 1000.0),
        vec![100.0, 500.0, -1.0, 149.0, 450.0, 2.0, 100.0, 200.0, 400.0, 500.0]
    );
    assert_eq!(ann.gene_names_in("chrX", 0.0, 1000.0), vec!["G".to_string()]);
    assert!(ann.genes_in("chrX", 500.0, 600.0).is_empty(), "end is exclusive");
    assert!(ann.genes_in("chr1", 0.0, 1000.0).is_empty(), "seqid must match");
    assert_eq!(ann.regions_in("chrX", 0.0, 1000.0), vec![250.0, 300.0, 0.0]);
    assert_eq!(ann.region_names_in("chrX", 0.0, 1000.0), vec!["del".to_string()]);
}
