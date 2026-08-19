//! Putting traces from different pores on a comparable scale.

use alloc::vec::Vec;

/// How to estimate the shift and scale for a trace.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub enum Scaling {
    /// Median and median absolute deviation.
    ///
    /// The default, and what production basecallers use. A nanopore trace is
    /// full of genuine outliers — stalls, blockages, adapter — and mean/stddev
    /// chases them. The median does not.
    #[default]
    MedMad,
    /// Mean and standard deviation.
    ///
    /// Cheaper and single-pass, but sensitive to exactly the artifacts that
    /// nanopore traces are full of. Provided for comparison, not for production.
    MeanStd,
}

/// The shift and scale actually applied to a trace.
///
/// Kept as a value so a caller can record it, reuse it across reads from the
/// same pore, or invert it.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct ScalingParams {
    pub shift: f32,
    pub scale: f32,
}

impl ScalingParams {
    /// Standardise samples in place: `(x - shift) / scale`.
    pub fn apply(self, samples: &mut [f32]) {
        for s in samples.iter_mut() {
            *s = (*s - self.shift) / self.scale;
        }
    }

    /// Undo [`apply`](Self::apply).
    pub fn invert(self, samples: &mut [f32]) {
        for s in samples.iter_mut() {
            *s = *s * self.scale + self.shift;
        }
    }
}

/// Scale of 1.0 is the identity, and is what degenerate inputs fall back to.
///
/// A constant trace has zero deviation; dividing by it would produce NaN and
/// poison every downstream tensor, so we leave the samples alone instead.
const IDENTITY_SCALE: f32 = 1.0;

impl Scaling {
    /// Estimate parameters without modifying the samples.
    pub fn estimate(self, samples: &[f32]) -> ScalingParams {
        if samples.is_empty() {
            return ScalingParams {
                shift: 0.0,
                scale: IDENTITY_SCALE,
            };
        }
        match self {
            Scaling::MedMad => {
                let med = median(samples);
                // MAD scaled to be a consistent estimator of sigma under a
                // normal distribution, which is what the networks were
                // trained against.
                const MAD_TO_SIGMA: f32 = 1.4826;
                let devs: Vec<f32> = samples.iter().map(|&s| (s - med).abs()).collect();
                let mad = median(&devs) * MAD_TO_SIGMA;
                ScalingParams {
                    shift: med,
                    scale: if mad > 0.0 { mad } else { IDENTITY_SCALE },
                }
            }
            Scaling::MeanStd => {
                let n = samples.len() as f32;
                let mean = samples.iter().sum::<f32>() / n;
                let var = samples
                    .iter()
                    .map(|&s| (s - mean) * (s - mean))
                    .sum::<f32>()
                    / n;
                let sd = sqrt(var);
                ScalingParams {
                    shift: mean,
                    scale: if sd > 0.0 { sd } else { IDENTITY_SCALE },
                }
            }
        }
    }
}

/// Median by sorting a copy. Allocates, but keeps the code obvious; the hot
/// path here is dwarfed by the network that consumes the result.
fn median(xs: &[f32]) -> f32 {
    let mut v: Vec<f32> = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        // Average the two middle elements rather than picking one, so the
        // estimate is stable under reversal of the input.
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// Square root without `std`.
fn sqrt(x: f32) -> f32 {
    if x <= 0.0 {
        return 0.0;
    }
    // Newton-Raphson, seeded generously. Converges to f32 precision well
    // within this iteration count for the magnitudes seen in pA traces.
    let mut guess = x;
    for _ in 0..24 {
        guess = 0.5 * (guess + x / guess);
    }
    guess
}
