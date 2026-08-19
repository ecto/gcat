use gcat_signal::{Scaling, Signal};

fn ramp(n: usize) -> Signal {
    Signal::new((0..n).map(|i| i as f32).collect(), 4000)
}

#[test]
fn adc_conversion_applies_scale_and_offset() {
    let s = Signal::from_adc(&[0, 10, -10], 0.5, 2.0, 4000);
    assert_eq!(s.samples(), &[1.0, 6.0, -4.0]);
    assert_eq!(s.sample_rate(), 4000);
}

#[test]
fn duration_uses_sample_rate() {
    assert_eq!(ramp(4000).duration_secs(), Some(1.0));
    assert_eq!(Signal::new(vec![1.0], 0).duration_secs(), None);
}

#[test]
fn medmad_centres_on_the_median() {
    let mut s = Signal::new(vec![10.0, 11.0, 12.0, 13.0, 14.0], 4000);
    let p = s.rescale(Scaling::MedMad);
    assert_eq!(p.shift, 12.0);
    // Middle sample sits exactly at the centre after scaling.
    assert!(s.samples()[2].abs() < 1e-6);
}

#[test]
fn medmad_ignores_outliers_that_meanstd_chases() {
    // A stall: most of the trace is tight, one span is wildly off.
    let mut base = vec![10.0; 40];
    base.extend([900.0; 6]);

    let med = Scaling::MedMad.estimate(&base);
    let mean = Scaling::MeanStd.estimate(&base);

    assert_eq!(
        med.shift, 10.0,
        "median should sit on the bulk of the trace"
    );
    assert!(
        mean.shift > 100.0,
        "mean should be dragged by the outliers, got {}",
        mean.shift
    );
}

#[test]
fn constant_trace_does_not_produce_nan() {
    let mut s = Signal::new(vec![7.0; 32], 4000);
    let p = s.rescale(Scaling::MedMad);
    assert_eq!(p.scale, 1.0);
    assert!(s.samples().iter().all(|x| x.is_finite()));

    let mut s2 = Signal::new(vec![7.0; 32], 4000);
    s2.rescale(Scaling::MeanStd);
    assert!(s2.samples().iter().all(|x| x.is_finite()));
}

#[test]
fn empty_trace_is_handled() {
    let p = Scaling::MedMad.estimate(&[]);
    assert_eq!(p.scale, 1.0);
    assert_eq!(p.shift, 0.0);
}

#[test]
fn scaling_roundtrips() {
    let original: Vec<f32> = (0..64)
        .map(|i| (i as f32 * 0.7).sin() * 20.0 + 90.0)
        .collect();
    let mut work = original.clone();
    let p = Scaling::MedMad.estimate(&work);
    p.apply(&mut work);
    p.invert(&mut work);
    for (a, b) in original.iter().zip(&work) {
        assert!((a - b).abs() < 1e-3, "{a} vs {b}");
    }
}

#[test]
fn windows_cover_the_whole_trace() {
    let s = ramp(100);
    let ws: Vec<_> = s.windows(30, 10).collect();
    assert_eq!(ws[0].start(), 0);
    // stride = 30 - 10 = 20
    assert_eq!(ws[1].start(), 20);
    // The last window is short rather than dropped, so no signal is lost.
    let last = ws.last().unwrap();
    assert_eq!(last.start() + last.len(), 100);
}

#[test]
fn windows_of_exact_multiple_have_no_short_tail() {
    let s = ramp(90);
    let ws: Vec<_> = s.windows(30, 0).collect();
    assert_eq!(ws.len(), 3);
    assert!(ws.iter().all(|w| w.len() == 30));
}

#[test]
fn overlap_at_or_above_size_still_terminates() {
    let s = ramp(50);
    // Would be a zero stride if not clamped — must not hang or loop forever.
    let ws: Vec<_> = s.windows(10, 10).collect();
    assert!(ws.len() <= 50);
    let ws2: Vec<_> = s.windows(10, 99).collect();
    assert!(ws2.len() <= 50);
}

#[test]
fn degenerate_windows_yield_nothing() {
    assert_eq!(ramp(10).windows(0, 0).count(), 0);
    assert_eq!(Signal::new(vec![], 4000).windows(10, 2).count(), 0);
}

#[test]
fn window_contents_match_the_source() {
    let s = ramp(50);
    for w in s.windows(16, 4) {
        for (i, &v) in w.samples().iter().enumerate() {
            assert_eq!(v, (w.start() + i) as f32);
        }
    }
}

#[cfg(feature = "tang")]
#[test]
fn windows_stack_into_a_tang_tensor() {
    let s = Signal::new((0..100).map(|i| i as f32).collect(), 4000);
    let t = s.to_tensor(30, 10);
    // stride 20 over 100 samples -> starts at 0,20,40,60,80 = 5 windows
    assert_eq!(t.shape().dims(), &[5, 1, 30]);
    // First window's samples land contiguously at the front.
    assert_eq!(&t.data()[..4], &[0.0, 1.0, 2.0, 3.0]);
    // Second window starts at sample 20.
    assert_eq!(t.data()[30], 20.0);
    // The final window covers samples 80..100, so the last 10 slots are padding.
    let tail = &t.data()[4 * 30..];
    assert_eq!(tail[19], 99.0);
    assert!(
        tail[20..].iter().all(|&x| x == 0.0),
        "tail must be zero-padded"
    );
}

#[cfg(feature = "tang")]
#[test]
fn exact_multiple_needs_no_padding() {
    let s = Signal::new((0..90).map(|i| i as f32).collect(), 4000);
    let t = s.to_tensor(30, 0);
    assert_eq!(t.shape().dims(), &[3, 1, 30]);
    assert_eq!(t.data()[89], 89.0);
}
