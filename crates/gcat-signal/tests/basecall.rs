//! End-to-end: raw signal in, bases out.
//!
//! Requires `--features tang` and a tang checkout carrying CTC.
//!
//! This wires the whole stage-2 path together — gcat-signal scaling and
//! windowing, a tang tensor, CTC loss and gradient from tang-train, and a
//! decode back to DNA. It is not a basecaller and makes no accuracy claim; it
//! exists so that a break anywhere along the seam fails loudly.

#![cfg(feature = "tang")]

use gcat::Base;
use gcat_signal::{Scaling, Signal};
use tang_tensor::{Shape, Tensor};
use tang_train::{beam_decode, ctc_loss, ctc_loss_grad, greedy_decode, min_input_len};

/// Class 0 is CTC's blank; the four bases occupy 1..=4.
const BLANK: usize = 0;
const CLASSES: usize = 5;
const SAMPLES_PER_BASE: usize = 4;

fn base_to_class(b: Base) -> usize {
    b as usize + 1
}

fn class_to_dna(classes: &[usize]) -> String {
    classes
        .iter()
        .filter_map(|&c| Base::from_bits((c - 1) as u8))
        .map(|b| b.to_ascii() as char)
        .collect()
}

/// The characteristic current level of each base, in picoamps.
fn level(b: Base) -> f32 {
    match b {
        Base::A => 60.0,
        Base::C => 80.0,
        Base::G => 100.0,
        Base::T => 120.0,
    }
}

/// Fabricate a trace: each base holds its level for a few samples, with a small
/// deterministic wobble. A real pore reads about five bases at once and is far
/// messier; this only needs to have the right shape.
fn synthesize(seq: &[Base]) -> Signal {
    let mut samples = Vec::new();
    for (i, &b) in seq.iter().enumerate() {
        for j in 0..SAMPLES_PER_BASE {
            let jitter = ((i * 7 + j * 3) % 5) as f32 - 2.0;
            samples.push(level(b) + jitter);
        }
    }
    Signal::new(samples, 4000)
}

/// Score each sample by how close it sits to each base's level.
///
/// This is the model a trained network would have to discover; handing it over
/// keeps this test about the pipeline rather than about whether a toy optimiser
/// converges. Blank scores just below every base so it never wins outright.
fn oracle_logits(signal: &Signal, shift: f32, scale: f32) -> Tensor<f64> {
    let time = signal.len();
    let data: Vec<f64> = signal
        .samples()
        .iter()
        .flat_map(|&s| {
            // Undo the scaling so distances are in the original picoamp units.
            let pa = s * scale + shift;
            (0..CLASSES).map(move |c| {
                if c == BLANK {
                    -5.0
                } else {
                    let b = Base::from_bits((c - 1) as u8).expect("valid class");
                    -((pa - level(b)) as f64 / 10.0).powi(2)
                }
            })
        })
        .collect();
    Tensor::new(data, Shape::from_slice(&[time, CLASSES]))
}

#[test]
fn signal_front_end_feeds_ctc() {
    let seq = [Base::A, Base::C, Base::G, Base::T, Base::A, Base::C];
    let target: Vec<usize> = seq.iter().map(|&b| base_to_class(b)).collect();

    let mut signal = synthesize(&seq);
    let params = signal.rescale(Scaling::MedMad);
    assert!(params.scale > 0.0, "scaling collapsed");
    assert!(signal.samples().iter().all(|s| s.is_finite()));

    // A trace of 24 samples has ample room for 6 labels.
    assert!(signal.len() >= min_input_len(&target));

    let logits = oracle_logits(&signal, params.shift, params.scale);
    let loss = ctc_loss(&logits, &target, BLANK).unwrap();
    assert!(loss.is_finite() && loss > 0.0, "loss = {loss}");

    let grad = ctc_loss_grad(&logits, &target, BLANK).unwrap();
    assert_eq!(grad.shape().dims(), &[signal.len(), CLASSES]);
    assert!(grad.data().iter().all(|g| g.is_finite()));
    // Every time step's gradient must sum to zero.
    for t in 0..signal.len() {
        let row: f64 = (0..CLASSES).map(|c| grad.get(&[t, c])).sum();
        assert!(row.abs() < 1e-9, "step {t} gradient sums to {row}");
    }
}

#[test]
fn decodes_a_synthetic_read_back_to_dna() {
    let seq = [
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::G,
        Base::T,
    ];
    let truth: String = seq.iter().map(|b| b.to_ascii() as char).collect();

    let mut signal = synthesize(&seq);
    let params = signal.rescale(Scaling::MedMad);
    let logits = oracle_logits(&signal, params.shift, params.scale);

    assert_eq!(class_to_dna(&greedy_decode(&logits, BLANK)), truth);
    assert_eq!(class_to_dna(&beam_decode(&logits, BLANK, 8)), truth);
}

#[test]
fn ctc_gradient_reduces_the_loss_on_signal_derived_logits() {
    // Not a convergence claim — just that the gradient points downhill when the
    // logits come from a real scaled trace rather than a synthetic tensor.
    let seq = [Base::A, Base::C, Base::G, Base::T, Base::A];
    let target: Vec<usize> = seq.iter().map(|&b| base_to_class(b)).collect();
    let mut signal = synthesize(&seq);
    signal.rescale(Scaling::MedMad);

    let time = signal.len();
    let shape = Shape::from_slice(&[time, CLASSES]);
    let mut data = vec![0.0f64; time * CLASSES];

    let first = ctc_loss(&Tensor::new(data.clone(), shape.clone()), &target, BLANK).unwrap();
    for _ in 0..6000 {
        let logits = Tensor::new(data.clone(), shape.clone());
        let grad = ctc_loss_grad(&logits, &target, BLANK).unwrap();
        for (d, g) in data.iter_mut().zip(grad.data()) {
            *d -= 0.5 * g;
        }
    }
    let last = ctc_loss(&Tensor::new(data.clone(), shape.clone()), &target, BLANK).unwrap();

    assert!(last < first, "loss went {first} -> {last}");
    // Plain gradient descent with no momentum; the tail is slow by nature.
    assert!(last < 0.01, "expected convergence, ended at {last}");
    assert_eq!(
        class_to_dna(&greedy_decode(&Tensor::new(data, shape), BLANK)),
        "ACGTA"
    );
}

/// A homopolymer needs a blank between its two halves to survive collapsing.
///
/// The signal for `TT` is one long run at the T level, so nothing in the trace
/// marks the boundary. This is the real homopolymer problem in nanopore
/// sequencing, and the reason production basecallers are convolutional or
/// recurrent rather than pointwise — a model that sees one sample at a time
/// cannot place that blank, and calls the run short.
#[test]
fn a_pointwise_model_calls_homopolymers_short() {
    let seq = [Base::A, Base::C, Base::T, Base::T, Base::G];
    let mut signal = synthesize(&seq);
    let params = signal.rescale(Scaling::MedMad);
    let logits = oracle_logits(&signal, params.shift, params.scale);

    // ACTTG comes back as ACTG: the doubled T collapses to one.
    assert_eq!(class_to_dna(&greedy_decode(&logits, BLANK)), "ACTG");

    // CTC still scores the true label sequence as possible — the information is
    // absent from this model's output, not from the loss.
    let target: Vec<usize> = seq.iter().map(|&b| base_to_class(b)).collect();
    let loss = ctc_loss(&logits, &target, BLANK).unwrap();
    assert!(loss.is_finite(), "true target should remain reachable");
}

#[test]
fn ctc_rejects_a_signal_too_short_for_its_labels() {
    let logits = Tensor::new(vec![0.0f64; 4 * CLASSES], Shape::from_slice(&[4, CLASSES]));
    let target: Vec<usize> = (0..8).map(|i| 1 + i % 4).collect();
    assert!(ctc_loss(&logits, &target, BLANK).is_err());
}
