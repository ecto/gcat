//! Visualization kernels — chaos game, GC skew, 2-D histograms.
//!
//! Every kernel here is a histogram: bounded integer accumulators, updated by a
//! rule that depends only on local input. That is deliberate. Each one maps to
//! a GPU storage buffer plus an atomic increment without changing shape, so the
//! same code can move to a compute shader later. Nothing allocates inside a hot
//! loop and nothing needs the whole input at once.
//!
//! `#![no_std]` with `alloc`, so these run in a browser, on a laptop, or on an
//! embedded box.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

mod cgr;
mod color;
mod hist;
mod skew;

pub use cgr::{Cgr, MAX_ORDER};
pub use color::{colorize, Palette};
pub use hist::Hist2d;
pub use skew::{GcSkew, SkewPoint};
