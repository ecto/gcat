use gcat_align::{align, Fate};

/// A deterministic pseudo-random sequence, so tests do not depend on luck.
fn random_seq(n: usize, seed: u64) -> Vec<u8> {
    let mut s = seed;
    (0..n)
        .map(|_| {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(s >> 33) as usize % 4]
        })
        .collect()
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn count(f: &[Fate], which: Fate) -> usize {
    f.iter().filter(|&&x| x == which).count()
}

#[test]
fn identical_sequences_match_everywhere() {
    let r = random_seq(2000, 1);
    let a = align(&r, &r);
    assert_eq!(count(&a.fates, Fate::Match), 2000);
    assert!(!a.reverse);
}

#[test]
fn a_deletion_in_the_query_reads_as_gaps() {
    let r = random_seq(3000, 2);
    let q: Vec<u8> = r[..1000].iter().chain(&r[1400..]).copied().collect();
    let a = align(&r, &q);
    // Bases either side of the cut can match by chance, which slides the gap
    // a few places with an equally good score; its size is what is fixed.
    assert_eq!(count(&a.fates, Fate::Gap), 400);
    assert!((990..1410).filter(|&i| a.fates[i] == Fate::Gap).count() == 400);
    assert_eq!(count(&a.fates, Fate::Match), 2600);
}

#[test]
fn point_mutations_read_as_mismatches() {
    let r = random_seq(3000, 3);
    let mut q = r.clone();
    for i in (100..2900).step_by(100) {
        q[i] = if q[i] == b'A' { b'C' } else { b'A' };
    }
    let a = align(&r, &q);
    assert_eq!(count(&a.fates, Fate::Mismatch), 28);
    assert_eq!(a.fates[100], Fate::Mismatch);
}

#[test]
fn reverse_strand_is_found() {
    let r = random_seq(2000, 4);
    let a = align(&r, &revcomp(&r));
    assert!(a.reverse);
    assert_eq!(count(&a.fates, Fate::Match), 2000);
}

#[test]
fn flanking_query_sequence_is_ignored_and_unrelated_is_unaligned() {
    let r = random_seq(2000, 5);
    let q: Vec<u8> = random_seq(5000, 6).into_iter().chain(r.iter().copied()).chain(random_seq(5000, 7)).collect();
    assert_eq!(count(&align(&r, &q).fates, Fate::Match), 2000);
    assert_eq!(count(&align(&r, &random_seq(2000, 8)).fates, Fate::Unaligned), 2000);
}

#[test]
fn identity_is_windowed_and_blank_where_unaligned() {
    let r = random_seq(3000, 9);
    let q: Vec<u8> = r[..1000].iter().chain(&r[1400..]).copied().collect();
    let id = align(&r, &q).identity(50);
    assert_eq!(id[500], Some(100.0));
    assert_eq!(id[1200], Some(0.0), "inside the deletion every base is lost");
    let blank = align(&r, &random_seq(3000, 10)).identity(50);
    assert!(blank.iter().all(Option::is_none));
}

#[test]
fn dot_plot_draws_the_diagonal_and_an_inversion() {
    use gcat_align::dot_plot;
    let a = random_seq(4000, 12);
    let n = 8;
    let same = dot_plot(&a, &a, n);
    for i in 0..n {
        assert!(same[i * n + i] > 0, "diagonal cell {i} empty");
    }
    let off: u32 = (0..n).flat_map(|y| (0..n).filter(move |&x| x != y).map(move |x| (y, x))).map(|(y, x)| same[y * n + x]).sum();
    let on: u32 = (0..n).map(|i| same[i * n + i]).sum();
    assert!(off * 50 < on, "random sequence should be nearly empty off the diagonal");

    // Invert the second half of b: its matches move to the reverse grid,
    // running against the diagonal.
    let b: Vec<u8> = a[..2000].iter().copied().chain(revcomp(&a[2000..])).collect();
    let plot = dot_plot(&a, &b, n);
    let rev = &plot[n * n..];
    assert!(rev[(n - 1) * n + n / 2] > 0, "start of the inverted half pairs with the end of b");
    assert!(rev[(n / 2) * n + n - 1] > 0);
    assert_eq!(rev[0], 0);
}
