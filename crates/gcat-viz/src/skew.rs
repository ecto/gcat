//! Cumulative GC skew along a sequence.
//!
//! Within each window, `(G - C) / (G + C)`. Accumulated along a bacterial
//! chromosome, the running total traces a V: replication runs in opposite
//! directions on the two halves of the circle, and the leading strand
//! accumulates G over C, so the skew flips sign where replication starts and
//! again where it ends.
//!
//! The practical upshot is that the minimum of the cumulative curve marks the
//! replication origin and the maximum marks the terminus, found by counting
//! letters — no alignment and no reference.

use alloc::vec::Vec;
use gcat::Base;

/// One window of the cumulative curve.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SkewPoint {
    /// Offset of the window's first base.
    pub offset: usize,
    /// Running total of per-window skew up to and including this window.
    pub cumulative: f32,
    /// This window's own skew, in `-1..=1`.
    pub window: f32,
    /// Fraction of this window that is G or C.
    pub gc: f32,
}

/// Accumulator for the cumulative GC skew curve.
pub struct GcSkew {
    window: usize,
    points: Vec<SkewPoint>,
    g: u32,
    c: u32,
    at: u32,
    seen: usize,
    offset: usize,
    cumulative: f32,
}

impl GcSkew {
    /// Accumulate over windows of `window` bases. A window of 0 is treated as 1.
    pub fn new(window: usize) -> GcSkew {
        GcSkew {
            window: window.max(1),
            points: Vec::new(),
            g: 0,
            c: 0,
            at: 0,
            seen: 0,
            offset: 0,
            cumulative: 0.0,
        }
    }

    /// Window size in bases.
    #[inline]
    pub fn window(&self) -> usize {
        self.window
    }

    /// The curve so far.
    #[inline]
    pub fn points(&self) -> &[SkewPoint] {
        &self.points
    }

    #[inline]
    pub fn push(&mut self, b: Base) {
        match b {
            Base::G => self.g += 1,
            Base::C => self.c += 1,
            _ => self.at += 1,
        }
        self.seen += 1;
        if self.seen == self.window {
            self.close_window();
        }
    }

    /// Feed raw ASCII, ignoring anything that is not a base.
    pub fn push_ascii(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            if let Some(b) = Base::from_ascii(byte) {
                self.push(b);
            }
        }
    }

    fn close_window(&mut self) {
        let gc = self.g + self.c;
        // A window with no G or C has no defined skew; treat it as flat rather
        // than dividing by zero and poisoning the running total.
        let window_skew = if gc > 0 {
            (self.g as f32 - self.c as f32) / gc as f32
        } else {
            0.0
        };
        self.cumulative += window_skew;
        let total = gc + self.at;
        self.points.push(SkewPoint {
            offset: self.offset,
            cumulative: self.cumulative,
            window: window_skew,
            gc: if total > 0 {
                gc as f32 / total as f32
            } else {
                0.0
            },
        });
        self.offset += self.seen;
        self.g = 0;
        self.c = 0;
        self.at = 0;
        self.seen = 0;
    }

    /// Close any partial trailing window so the curve reaches the end.
    pub fn finish(&mut self) {
        if self.seen > 0 {
            self.close_window();
        }
    }

    /// Offset of the cumulative minimum — the predicted replication origin.
    pub fn origin(&self) -> Option<usize> {
        self.extreme(|a, b| a < b)
    }

    /// Offset of the cumulative maximum — the predicted replication terminus.
    pub fn terminus(&self) -> Option<usize> {
        self.extreme(|a, b| a > b)
    }

    fn extreme(&self, better: impl Fn(f32, f32) -> bool) -> Option<usize> {
        let mut it = self.points.iter();
        let first = it.next()?;
        let mut best = first;
        for p in it {
            if better(p.cumulative, best.cumulative) {
                best = p;
            }
        }
        Some(best.offset)
    }

    /// Lowest and highest cumulative values, for scaling an axis.
    pub fn range(&self) -> Option<(f32, f32)> {
        let mut it = self.points.iter().map(|p| p.cumulative);
        let first = it.next()?;
        let mut lo = first;
        let mut hi = first;
        for v in it {
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
        Some((lo, hi))
    }
}
