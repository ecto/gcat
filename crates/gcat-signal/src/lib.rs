//! Raw nanopore signal — scaling, windowing, tensor handoff.
//!
//! A nanopore reports picoamp current samples, not bases. Everything upstream of
//! a basecaller network is this crate's job: put the signal on a comparable
//! scale across pores and flowcells, then cut it into the fixed-size,
//! overlapping windows the network expects.
//!
//! `#![no_std]` with `alloc`, so the same code runs on an embedded box next to
//! the sequencer and in a datacenter.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

use alloc::vec::Vec;

mod scale;
mod window;

pub use scale::{Scaling, ScalingParams};
pub use window::{Window, Windows};

/// A raw current trace from one pore, in picoamps.
///
/// This is the signal as the machine reports it: no scaling applied, sample
/// rate attached so downstream code never has to guess it.
#[derive(Clone, PartialEq, Debug)]
pub struct Signal {
    samples: Vec<f32>,
    sample_rate: u32,
}

impl Signal {
    /// Wrap raw samples with their sample rate in Hz.
    pub fn new(samples: Vec<f32>, sample_rate: u32) -> Signal {
        Signal {
            samples,
            sample_rate,
        }
    }

    /// Convert from the integer ADC values a device actually emits.
    ///
    /// `scale` and `offset` come from the run metadata:
    /// `picoamps = (adc + offset) * scale`.
    pub fn from_adc(adc: &[i16], scale: f32, offset: f32, sample_rate: u32) -> Signal {
        let samples = adc.iter().map(|&a| (a as f32 + offset) * scale).collect();
        Signal {
            samples,
            sample_rate,
        }
    }

    /// The samples.
    #[inline]
    pub fn samples(&self) -> &[f32] {
        &self.samples
    }

    /// Number of samples.
    #[inline]
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// True if there are no samples.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Sample rate in Hz.
    #[inline]
    pub fn sample_rate(&self) -> u32 {
        self.sample_rate
    }

    /// Duration in seconds, or `None` if the sample rate is zero.
    pub fn duration_secs(&self) -> Option<f32> {
        if self.sample_rate == 0 {
            return None;
        }
        Some(self.samples.len() as f32 / self.sample_rate as f32)
    }

    /// Apply a scaling in place and return the parameters that were used.
    pub fn rescale(&mut self, scaling: Scaling) -> ScalingParams {
        let params = scaling.estimate(&self.samples);
        params.apply(&mut self.samples);
        params
    }

    /// Cut into overlapping fixed-size windows for a basecaller.
    pub fn windows(&self, size: usize, overlap: usize) -> Windows<'_> {
        Windows::new(&self.samples, size, overlap)
    }
}

#[cfg(feature = "tang")]
mod tensor {
    use super::Signal;
    use tang_tensor::{Shape, Tensor};

    impl Signal {
        /// Stack windows into a `[batch, 1, size]` tensor for the compute stack.
        ///
        /// The channel dimension is 1 because a pore reports a single current
        /// trace; it is kept explicit so the shape matches what `Conv1d` wants
        /// without a reshape at the call site.
        ///
        /// A tensor must be rectangular but the final window of a trace is
        /// usually short, so the tail is zero-padded to `size`. Zero is the
        /// mean of a scaled trace, so the padding reads as "no signal" rather
        /// than as a spurious current level. Callers that need to distinguish
        /// real samples from padding should use [`Signal::windows`] directly,
        /// where the final window keeps its true length.
        pub fn to_tensor(&self, size: usize, overlap: usize) -> Tensor<f32> {
            let windows: alloc::vec::Vec<_> = self.windows(size, overlap).collect();
            let batch = windows.len();
            let mut data = alloc::vec::Vec::with_capacity(batch * size);
            for w in &windows {
                data.extend_from_slice(w.samples());
                data.resize(data.len() + (size - w.len()), 0.0);
            }
            debug_assert_eq!(data.len(), batch * size);
            Tensor::new(data, Shape::from_slice(&[batch, 1, size]))
        }
    }
}
