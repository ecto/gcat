//! A 2-D histogram, for read length against quality.
//!
//! Length and quality are each uninformative alone — the diagnostic question is
//! whether the long reads are also the good ones. Binning both at once answers
//! it, and unlike a scatter plot it does not fall over at ten million reads,
//! because cost is bounded by the bin count rather than the read count.

use alloc::{vec, vec::Vec};

/// A rectangular histogram with a log-scaled x axis.
pub struct Hist2d {
    cols: usize,
    rows: usize,
    counts: Vec<u32>,
    x_min: f32,
    x_max: f32,
    y_min: f32,
    y_max: f32,
    peak: u32,
    total: u64,
}

impl Hist2d {
    /// Bins spanning `x_min..x_max` (log scale) and `y_min..y_max` (linear).
    pub fn new(cols: usize, rows: usize, x_min: f32, x_max: f32, y_min: f32, y_max: f32) -> Hist2d {
        let cols = cols.max(1);
        let rows = rows.max(1);
        Hist2d {
            cols,
            rows,
            counts: vec![0u32; cols * rows],
            // Read lengths span orders of magnitude, so the x axis is
            // logarithmic; guard the lower bound so ln stays finite.
            x_min: x_min.max(1.0),
            x_max: x_max.max(2.0),
            y_min,
            y_max,
            peak: 0,
            total: 0,
        }
    }

    /// A histogram sized for nanopore reads: 100 bp to 100 kb, Q0 to Q30.
    pub fn read_length_quality(cols: usize, rows: usize) -> Hist2d {
        Hist2d::new(cols, rows, 100.0, 100_000.0, 0.0, 30.0)
    }

    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    pub fn counts(&self) -> &[u32] {
        &self.counts
    }

    #[inline]
    pub fn peak(&self) -> u32 {
        self.peak
    }

    #[inline]
    pub fn total(&self) -> u64 {
        self.total
    }

    /// Add one observation, clamping to the edge bins rather than dropping it.
    ///
    /// Dropping out-of-range points silently understates the total, which is
    /// worse than a slightly heavy edge bin that the axis label explains.
    pub fn push(&mut self, x: f32, y: f32) {
        let lx = log(x.max(self.x_min));
        let l0 = log(self.x_min);
        let l1 = log(self.x_max);
        let fx = ((lx - l0) / (l1 - l0)).clamp(0.0, 1.0);
        let fy = ((y - self.y_min) / (self.y_max - self.y_min)).clamp(0.0, 1.0);

        let col = ((fx * self.cols as f32) as usize).min(self.cols - 1);
        let row_from_bottom = ((fy * self.rows as f32) as usize).min(self.rows - 1);
        // Row 0 is the top of the image, so high quality sits high.
        let row = self.rows - 1 - row_from_bottom;

        let idx = row * self.cols + col;
        let c = self.counts[idx] + 1;
        self.counts[idx] = c;
        if c > self.peak {
            self.peak = c;
        }
        self.total += 1;
    }

    /// Axis bounds as `(x_min, x_max, y_min, y_max)`.
    pub fn bounds(&self) -> (f32, f32, f32, f32) {
        (self.x_min, self.x_max, self.y_min, self.y_max)
    }
}

/// Natural log without `std`.
fn log(x: f32) -> f32 {
    // Split into exponent and mantissa, then use atanh's series on the
    // mantissa, which converges fast over the half-open range it lands in.
    if x <= 0.0 {
        return 0.0;
    }
    let bits = x.to_bits();
    let exp = ((bits >> 23) & 0xff) as i32 - 127;
    let mantissa = f32::from_bits((bits & 0x007f_ffff) | 0x3f80_0000);
    let z = (mantissa - 1.0) / (mantissa + 1.0);
    let z2 = z * z;
    let mut series = 0.0f32;
    let mut term = z;
    let mut k = 1.0f32;
    for _ in 0..12 {
        series += term / k;
        term *= z2;
        k += 2.0;
    }
    core::f32::consts::LN_2 * exp as f32 + 2.0 * series
}
