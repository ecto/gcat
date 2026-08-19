//! Chaos game representation.
//!
//! The classic description is a random walk: start at the centre of a square
//! with a base at each corner, and for every base in the sequence jump halfway
//! toward its corner, marking where you land.
//!
//! That description hides the useful fact. After `k` jumps, the position
//! depends *only on the last `k` bases* — earlier bases have been halved into
//! irrelevance. So a `2^k` by `2^k` chaos game image is exactly the histogram
//! of every k-mer in the sequence, laid out so that k-mers sharing a suffix land
//! near each other. Zooming in by one level is looking at k-mers one base
//! longer.
//!
//! Two consequences. The cell index rolls in `O(1)` per base like any k-mer,
//! and because each base's cell is a pure function of its window rather than of
//! the walk so far, the whole thing is a parallel histogram rather than a
//! sequential simulation — which is what makes it portable to a GPU.

use alloc::{vec, vec::Vec};
use gcat::Base;

/// Largest supported order. `k = 15` is a 32768² grid, far past what any
/// display or sane memory budget wants.
pub const MAX_ORDER: usize = 15;

/// A chaos game accumulator over a `2^order` square grid.
pub struct Cgr {
    order: usize,
    size: usize,
    counts: Vec<u32>,
    /// Column bits of the current window, newest base in the high bit.
    col: usize,
    /// Row bits of the current window.
    row: usize,
    /// Bases seen since the last discontinuity, so partial windows are skipped.
    filled: usize,
    total: u64,
    peak: u32,
}

/// Column bit for a base: 1 for the corners on the right edge (G, T).
#[inline]
const fn col_bit(b: Base) -> usize {
    (b as usize >> 1) & 1
}

/// Row bit for a base: 1 for the corners on the top edge (C, G).
#[inline]
const fn row_bit(b: Base) -> usize {
    ((b as usize >> 1) ^ b as usize) & 1
}

impl Cgr {
    /// A grid of `2^order` on a side. Orders above [`MAX_ORDER`] are clamped.
    pub fn new(order: usize) -> Cgr {
        let order = order.clamp(1, MAX_ORDER);
        let size = 1usize << order;
        Cgr {
            order,
            size,
            counts: vec![0u32; size * size],
            col: 0,
            row: 0,
            filled: 0,
            total: 0,
            peak: 0,
        }
    }

    /// Grid width and height in cells.
    #[inline]
    pub fn size(&self) -> usize {
        self.size
    }

    /// The k of the k-mers this grid counts.
    #[inline]
    pub fn order(&self) -> usize {
        self.order
    }

    /// Raw counts, row-major, `size * size` long. Row 0 is the top of the image.
    #[inline]
    pub fn counts(&self) -> &[u32] {
        &self.counts
    }

    /// Number of k-mers counted.
    #[inline]
    pub fn total(&self) -> u64 {
        self.total
    }

    /// Highest count in any single cell, for scaling a colour ramp.
    #[inline]
    pub fn peak(&self) -> u32 {
        self.peak
    }

    /// Break the rolling window.
    ///
    /// Call between records. Without it, the last bases of one read and the
    /// first of the next would form k-mers that exist in neither.
    #[inline]
    pub fn discontinuity(&mut self) {
        self.filled = 0;
    }

    /// Feed one base.
    #[inline]
    pub fn push(&mut self, b: Base) {
        let mask = self.size - 1;
        let high = self.order - 1;
        // Every existing base halves in weight; the new base takes the top bit.
        self.col = ((self.col >> 1) | (col_bit(b) << high)) & mask;
        self.row = ((self.row >> 1) | (row_bit(b) << high)) & mask;

        if self.filled < self.order {
            self.filled += 1;
            if self.filled < self.order {
                return;
            }
        }
        // Row 0 is the top of the image, but y grows upward in the chaos game.
        let img_row = mask - self.row;
        let idx = img_row * self.size + self.col;
        let c = self.counts[idx] + 1;
        self.counts[idx] = c;
        if c > self.peak {
            self.peak = c;
        }
        self.total += 1;
    }

    /// Feed a run of bases.
    pub fn extend(&mut self, bases: impl IntoIterator<Item = Base>) {
        for b in bases {
            self.push(b);
        }
    }

    /// Feed raw ASCII, breaking the window at anything that is not a base.
    ///
    /// Ambiguity codes and line noise end the current window rather than being
    /// substituted, so no k-mer is invented that the sequence did not contain.
    pub fn push_ascii(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            match Base::from_ascii(byte) {
                Some(b) => self.push(b),
                None => self.discontinuity(),
            }
        }
    }

    /// Fraction of cells that were never hit.
    ///
    /// High emptiness at a given order means the sequence avoids many k-mers —
    /// the signature that makes chaos game images of different genomes look
    /// different.
    pub fn sparsity(&self) -> f32 {
        let empty = self.counts.iter().filter(|&&c| c == 0).count();
        empty as f32 / self.counts.len() as f32
    }
}
