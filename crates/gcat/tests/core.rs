use gcat::{Base, Kmer, Seq, SeqError};

#[test]
fn packs_four_bases_per_byte() {
    let s: Seq = "ACGTACGTA".parse().unwrap();
    assert_eq!(s.len(), 9);
    assert_eq!(s.packed_bytes(), 3);
}

#[test]
fn roundtrips_through_ascii() {
    let src = "ACGTTTGACCAGTGCA";
    let s: Seq = src.parse().unwrap();
    assert_eq!(s.to_string(), src);
    assert_eq!(String::from_utf8(s.to_ascii()).unwrap(), src);
}

#[test]
fn indexes_every_position() {
    let src = "ACGTACGTACGTAC";
    let s: Seq = src.parse().unwrap();
    for (i, c) in src.bytes().enumerate() {
        assert_eq!(s.get(i), Base::from_ascii(c), "mismatch at {i}");
    }
    assert_eq!(s.get(src.len()), None);
}

#[test]
fn revcomp_is_an_involution() {
    let s: Seq = "AACCGGTTACGT".parse().unwrap();
    assert_eq!(s.revcomp().revcomp(), s);
    assert_eq!(s.revcomp().to_string(), "ACGTAACCGGTT");
}

#[test]
fn lowercase_and_uracil_parse() {
    let s: Seq = "acguACGU".parse().unwrap();
    assert_eq!(s.to_string(), "ACGTACGT");
}

#[test]
fn ambiguity_is_distinguished_from_garbage() {
    assert_eq!(
        "ACGNT".parse::<Seq>(),
        Err(SeqError::Ambiguous {
            offset: 3,
            byte: b'N'
        })
    );
    assert_eq!(
        "ACG?T".parse::<Seq>(),
        Err(SeqError::Invalid {
            offset: 3,
            byte: b'?'
        })
    );
}

#[test]
fn lossy_parse_reports_substitutions() {
    let (s, n) = Seq::from_ascii_lossy(b"ACGNTNN", Base::A).unwrap();
    assert_eq!(n, 3);
    assert_eq!(s.to_string(), "ACGATAA");
    assert!(Seq::from_ascii_lossy(b"AC!G", Base::A).is_err());
}

#[test]
fn gc_content_counts_g_and_c() {
    let s: Seq = "GGCCAATT".parse().unwrap();
    assert_eq!(s.gc_content(), Some(0.5));
    assert_eq!(Seq::new().gc_content(), None);
}

#[test]
fn kmers_slide_one_base_at_a_time() {
    let s: Seq = "ACGTA".parse().unwrap();
    let got: Vec<String> = s.kmers(3).map(|k| k.to_string()).collect();
    assert_eq!(got, ["ACG", "CGT", "GTA"]);
}

#[test]
fn kmer_count_and_size_hint_agree() {
    let s: Seq = "ACGTACGTACGT".parse().unwrap();
    for k in 1..=8 {
        let expect = s.len() - k + 1;
        assert_eq!(s.kmers(k).size_hint(), (expect, Some(expect)));
        assert_eq!(s.kmers(k).count(), expect, "k={k}");
    }
}

#[test]
fn degenerate_k_yields_nothing() {
    let s: Seq = "ACGT".parse().unwrap();
    assert_eq!(s.kmers(0).count(), 0);
    assert_eq!(s.kmers(5).count(), 0);
    assert_eq!(s.kmers(33).count(), 0);
}

#[test]
fn rolling_matches_rebuilding_from_scratch() {
    let s: Seq = "ACGTTGCAACGTTGCAACGT".parse().unwrap();
    let k = 7;
    let rolled: Vec<String> = s.kmers(k).map(|x| x.to_string()).collect();
    let rebuilt: Vec<String> = (0..=s.len() - k)
        .map(|i| {
            let bases: Vec<Base> = (i..i + k).map(|j| s.get(j).unwrap()).collect();
            Kmer::from_bases(&bases).unwrap().to_string()
        })
        .collect();
    assert_eq!(rolled, rebuilt);
}

#[test]
fn kmer_revcomp_matches_seq_revcomp() {
    let s: Seq = "ACGTTGCA".parse().unwrap();
    let km = Kmer::from_bases(&s.iter().collect::<Vec<_>>()).unwrap();
    assert_eq!(km.revcomp().to_string(), s.revcomp().to_string());
}

#[test]
fn canonical_is_strand_independent() {
    let fwd: Seq = "ACGTTGCAT".parse().unwrap();
    let rev = fwd.revcomp();
    let f: Vec<_> = fwd.kmers(5).map(|k| k.canonical()).collect();
    let mut r: Vec<_> = rev.kmers(5).map(|k| k.canonical()).collect();
    r.reverse();
    assert_eq!(f, r);
}

#[test]
fn max_k_does_not_overflow() {
    let src = "ACGT".repeat(9);
    let s: Seq = src.parse().unwrap();
    let ks: Vec<_> = s.kmers(32).collect();
    assert_eq!(ks.len(), 5);
    assert_eq!(ks[0].k(), 32);
    assert_eq!(ks[0].to_string().len(), 32);
    assert_eq!(ks[4].to_string(), &src[4..36]);
}
