use gcat::{Base, Seq};
use gcat_viz::{colorize, Cgr, GcSkew, Hist2d, Palette};

#[test]
fn cgr_cell_is_the_kmer() {
    // Order 3 means every cell is one 3-mer, and there are 64 of them.
    let mut grid = Cgr::new(3);
    assert_eq!(grid.size(), 8);
    assert_eq!(grid.counts().len(), 64);

    // Feed every 3-mer exactly once and each cell must hold exactly one count.
    for a in Base::ALL {
        for b in Base::ALL {
            for c in Base::ALL {
                grid.discontinuity();
                grid.extend([a, b, c]);
            }
        }
    }
    assert_eq!(grid.total(), 64);
    assert!(
        grid.counts().iter().all(|&c| c == 1),
        "each 3-mer should occupy its own distinct cell"
    );
    assert_eq!(grid.peak(), 1);
    assert_eq!(grid.sparsity(), 0.0);
}

#[test]
fn cgr_corners_are_the_homopolymers() {
    // Poly-A drives to one corner and stays there.
    let mut grid = Cgr::new(4);
    for _ in 0..100 {
        grid.push(Base::A);
    }
    let hits: Vec<usize> = grid
        .counts()
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, _)| i)
        .collect();
    assert_eq!(hits.len(), 1, "a homopolymer occupies exactly one cell");
    // A is the corner at column 0, and y=0 which renders at the bottom row.
    let size = grid.size();
    assert_eq!(
        hits[0],
        (size - 1) * size,
        "poly-A belongs in the bottom-left"
    );
}

#[test]
fn cgr_counts_every_window_once() {
    let seq: Seq = "ACGTACGTTGCAACGT".parse().unwrap();
    let order = 4;
    let mut grid = Cgr::new(order);
    grid.extend(seq.iter());
    assert_eq!(grid.total() as usize, seq.len() - order + 1);
}

#[test]
fn cgr_discontinuity_prevents_invented_kmers() {
    let mut joined = Cgr::new(4);
    joined.push_ascii(b"ACGTACGT");
    joined.push_ascii(b"TTTTGGGG");

    let mut split = Cgr::new(4);
    split.push_ascii(b"ACGTACGT");
    split.discontinuity();
    split.push_ascii(b"TTTTGGGG");

    // Joining spans the seam and manufactures 3 k-mers that exist in neither read.
    assert_eq!(joined.total(), split.total() + 3);
}

#[test]
fn cgr_treats_ambiguity_as_a_break() {
    let mut grid = Cgr::new(4);
    grid.push_ascii(b"ACGTNACGT");
    // Each 4-base run yields exactly one 4-mer; the N severs them.
    assert_eq!(grid.total(), 2);
}

#[test]
fn cgr_order_is_clamped() {
    assert_eq!(Cgr::new(0).order(), 1);
    assert_eq!(Cgr::new(999).order(), gcat_viz::MAX_ORDER);
}

#[test]
fn skew_is_flat_for_balanced_sequence() {
    let mut skew = GcSkew::new(4);
    skew.push_ascii(b"GCGCGCGCGCGC");
    skew.finish();
    for p in skew.points() {
        assert!(
            p.window.abs() < 1e-6,
            "balanced window skewed by {}",
            p.window
        );
    }
}

#[test]
fn skew_tracks_a_deliberate_switch() {
    // G-rich then C-rich: the cumulative curve should rise, then fall.
    let mut skew = GcSkew::new(10);
    skew.push_ascii(&b"G".repeat(50));
    skew.push_ascii(&b"C".repeat(50));
    skew.finish();

    let pts = skew.points();
    assert_eq!(pts.len(), 10);
    assert!(
        (pts[4].cumulative - 5.0).abs() < 1e-5,
        "G-rich half should climb"
    );
    assert!(
        pts[9].cumulative.abs() < 1e-5,
        "C-rich half should cancel it"
    );
    // The maximum sits at the switch.
    assert_eq!(skew.terminus(), Some(40));
}

#[test]
fn skew_survives_a_window_with_no_gc() {
    let mut skew = GcSkew::new(4);
    skew.push_ascii(b"AAAATTTTGGGG");
    skew.finish();
    assert!(skew.points().iter().all(|p| p.cumulative.is_finite()));
    assert_eq!(skew.points()[0].window, 0.0);
    assert_eq!(skew.points()[0].gc, 0.0);
    assert_eq!(skew.points()[2].gc, 1.0);
}

#[test]
fn skew_finish_closes_a_partial_window() {
    let mut skew = GcSkew::new(10);
    skew.push_ascii(b"GGGGG");
    assert_eq!(skew.points().len(), 0);
    skew.finish();
    assert_eq!(skew.points().len(), 1);
}

#[test]
fn hist_bins_and_clamps() {
    let mut h = Hist2d::new(10, 10, 1.0, 1000.0, 0.0, 10.0);
    h.push(1.0, 0.0);
    h.push(1000.0, 10.0);
    // Out of range in both directions: clamped to edges, never dropped.
    h.push(0.001, -50.0);
    h.push(1e9, 1e9);
    assert_eq!(h.total(), 4);

    let cols = h.cols();
    let rows = h.rows();
    // Low length, low quality -> bottom-left.
    assert_eq!(h.counts()[(rows - 1) * cols], 2);
    // High length, high quality -> top-right.
    assert_eq!(h.counts()[cols - 1], 2);
}

#[test]
fn hist_log_axis_puts_the_midpoint_in_the_middle() {
    let mut h = Hist2d::new(11, 1, 1.0, 10000.0, 0.0, 1.0);
    // 100 is the geometric midpoint of 1 and 10000.
    h.push(100.0, 0.5);
    let mid = h.counts()[5];
    assert_eq!(
        mid, 1,
        "geometric midpoint should land in the centre column"
    );
}

#[test]
fn colorize_fills_the_buffer() {
    let counts = [0u32, 5, 10];
    let mut rgba = [0u8; 12];
    colorize(&counts, 10, Palette::Ember, 0.5, &mut rgba);
    // Every pixel opaque.
    assert_eq!(rgba[3], 255);
    assert_eq!(rgba[7], 255);
    assert_eq!(rgba[11], 255);
    // Intensity increases with count.
    let lum = |i: usize| rgba[i * 4] as u32 + rgba[i * 4 + 1] as u32 + rgba[i * 4 + 2] as u32;
    assert!(lum(0) < lum(1), "zero should be darkest");
    assert!(lum(1) < lum(2), "peak should be brightest");
}

#[test]
fn colorize_handles_an_empty_grid() {
    let counts = [0u32; 4];
    let mut rgba = [7u8; 16];
    colorize(&counts, 0, Palette::Ember, 0.5, &mut rgba);
    assert!(
        rgba.chunks(4).all(|p| p[3] == 255),
        "must not divide by a zero peak"
    );
}

#[test]
fn colorize_respects_a_short_buffer() {
    let counts = [1u32; 8];
    let mut rgba = [0u8; 8];
    colorize(&counts, 1, Palette::Ember, 1.0, &mut rgba);
    // Only two pixels fit; the call must not run past the end.
    assert_eq!(rgba[3], 255);
    assert_eq!(rgba[7], 255);
}

/// The no_std ln/exp/powf are hand-rolled, so check them against std.
#[test]
fn gamma_curve_matches_std_powf() {
    for &gamma in &[0.25f32, 0.4, 0.5, 1.0, 2.0] {
        for i in 0..=100 {
            let x = i as f32 / 100.0;
            let counts = [(x * 1000.0) as u32];
            let mut rgba = [0u8; 4];
            colorize(&counts, 1000, Palette::Mono, gamma, &mut rgba);

            let expect_t = if counts[0] == 0 {
                0.0
            } else {
                (counts[0] as f32 / 1000.0).powf(gamma)
            };
            // Mono ramp is linear in luminance, so invert it to recover t.
            let got = rgba[0] as f32;
            let want = if expect_t < 0.5 {
                16.0 + (128.0 - 16.0) * (expect_t * 2.0)
            } else {
                128.0 + (245.0 - 128.0) * ((expect_t - 0.5) * 2.0)
            };
            assert!(
                (got - want).abs() <= 2.0,
                "gamma {gamma} at x={x}: got {got}, want {want}"
            );
        }
    }
}
