//! Checks against real reference genomes.
//!
//! These run only when `scripts/fetch-data.sh` has been run; without the data
//! they skip rather than fail, so a fresh clone still gets a green suite.

use gcat::Seq;
use gcat_io::Reader;
use std::io::BufReader;
use std::path::PathBuf;

fn data(name: &str) -> Option<PathBuf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(name);
    p.exists().then_some(p)
}

fn load(name: &str) -> Option<(String, Seq)> {
    let path = data(name)?;
    let file = std::fs::File::open(&path).expect("data file readable");
    let mut r = Reader::new(BufReader::new(file));
    let rec = r.next_record().expect("parses").expect("has a record");
    Some((rec.id.to_string(), rec.to_seq().expect("unambiguous")))
}

/// Published length and GC content for each reference.
///
/// These are the canonical figures for these accessions; if gcat disagrees,
/// gcat is wrong.
const REFERENCES: &[(&str, usize, f64)] = &[
    ("phix174.fa", 5386, 44.76),
    ("lambda.fa", 48502, 49.86),
    ("sars2.fa", 29903, 37.97),
    ("ecoli.fa", 4_641_652, 50.79),
];

#[test]
fn lengths_and_gc_match_published_values() {
    let mut checked = 0;
    for &(file, len, gc) in REFERENCES {
        let Some((id, seq)) = load(file) else {
            continue;
        };
        assert_eq!(seq.len(), len, "{file} ({id}) length");
        let got = seq.gc_content().expect("non-empty") * 100.0;
        assert!(
            (got - gc).abs() < 0.01,
            "{file} GC: expected {gc:.2}%, got {got:.2}%"
        );
        checked += 1;
    }
    if checked == 0 {
        eprintln!("skipped: run scripts/fetch-data.sh to enable genome checks");
    }
}

#[test]
fn revcomp_roundtrips_on_a_real_genome() {
    let Some((_, seq)) = load("ecoli.fa") else {
        eprintln!("skipped: run scripts/fetch-data.sh");
        return;
    };
    assert_eq!(seq.revcomp().revcomp(), seq);
}

#[test]
fn packing_uses_a_quarter_of_the_bytes() {
    let Some((_, seq)) = load("ecoli.fa") else {
        eprintln!("skipped: run scripts/fetch-data.sh");
        return;
    };
    assert_eq!(seq.packed_bytes(), seq.len().div_ceil(4));
    assert!(seq.packed_bytes() < seq.len() / 3);
}

#[test]
fn every_kmer_of_lambda_is_distinct() {
    let Some((_, seq)) = load("lambda.fa") else {
        eprintln!("skipped: run scripts/fetch-data.sh");
        return;
    };
    // Lambda has no repeats long enough to duplicate a 21-mer, so the count of
    // distinct canonical 21-mers equals the number of windows.
    let kmers: std::collections::HashSet<u64> =
        seq.kmers(21).map(|k| k.canonical().bits()).collect();
    assert_eq!(kmers.len(), seq.len() - 21 + 1);
}
