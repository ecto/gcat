//! Kernel checks against real genomes. Skips without `scripts/fetch-data.sh`.

use gcat_io::Reader;
use gcat_viz::{Cgr, GcSkew};
use std::io::BufReader;
use std::path::PathBuf;

fn load(name: &str) -> Option<Vec<u8>> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(name);
    if !path.exists() {
        eprintln!("skipped: run scripts/fetch-data.sh");
        return None;
    }
    let file = std::fs::File::open(&path).expect("readable");
    let mut r = Reader::new(BufReader::new(file));
    let rec = r.next_record().expect("parses").expect("a record");
    Some(rec.seq.to_vec())
}

/// The headline claim: counting letters locates the replication origin.
///
/// E. coli K-12 MG1655 has a well-mapped oriC at 3,925,744 bp. The cumulative
/// GC skew minimum should land on it without alignment, reference, or database.
#[test]
fn gc_skew_finds_the_ecoli_replication_origin() {
    let Some(seq) = load("ecoli.fa") else { return };

    let window = 10_000;
    let mut skew = GcSkew::new(window);
    skew.push_ascii(&seq);
    skew.finish();

    const PUBLISHED_ORIC: usize = 3_925_744;
    const PUBLISHED_TER: usize = 1_590_000;

    let origin = skew.origin().expect("a curve");
    let terminus = skew.terminus().expect("a curve");

    let err = origin.abs_diff(PUBLISHED_ORIC);
    assert!(
        err <= window,
        "predicted oriC {origin} is {err} bp from the published {PUBLISHED_ORIC}"
    );
    // The terminus is a broader feature and less sharply defined; allow more room.
    assert!(
        terminus.abs_diff(PUBLISHED_TER) <= 20 * window,
        "predicted terminus {terminus} is far from the published {PUBLISHED_TER}"
    );
    assert!(origin > terminus, "origin and terminus came out swapped");
}

#[test]
fn skew_curve_covers_the_whole_genome() {
    let Some(seq) = load("ecoli.fa") else { return };
    let mut skew = GcSkew::new(10_000);
    skew.push_ascii(&seq);
    skew.finish();
    let last = skew.points().last().expect("points");
    assert!(
        last.offset + 10_000 >= seq.len() - 10_000,
        "curve stops short of the end"
    );
    // E. coli spans about 7.8 over 464 windows. A sequence with no
    // replication-strand bias would random-walk near zero, so anything past ~5
    // is the V rather than noise.
    let (lo, hi) = skew.range().expect("range");
    assert!(
        hi - lo > 5.0,
        "expected a pronounced V, got a span of {}",
        hi - lo
    );
}

/// Different genomes have visibly different k-mer occupancy — the property that
/// makes chaos game images distinguishable at a glance.
#[test]
fn chaos_game_separates_genomes_by_sparsity() {
    let (Some(ecoli), Some(phix)) = (load("ecoli.fa"), load("phix174.fa")) else {
        return;
    };

    let mut big = Cgr::new(8);
    big.push_ascii(&ecoli);
    let mut small = Cgr::new(8);
    small.push_ascii(&phix);

    // 4.6 Mb covers essentially every 8-mer; 5.4 kb cannot.
    assert!(
        big.sparsity() < 0.01,
        "E. coli should fill the grid, {} empty",
        big.sparsity()
    );
    assert!(
        small.sparsity() > 0.5,
        "phiX cannot cover 65536 cells with 5386 bases, {} empty",
        small.sparsity()
    );
    assert_eq!(big.total() as usize, ecoli.len() - 8 + 1);
}
