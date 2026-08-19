//! Cutting a trace into the overlapping windows a basecaller consumes.

/// One fixed-size window of a trace.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Window<'a> {
    samples: &'a [f32],
    start: usize,
}

impl<'a> Window<'a> {
    /// The samples in this window.
    #[inline]
    pub fn samples(&self) -> &'a [f32] {
        self.samples
    }

    /// Offset of this window's first sample within the full trace.
    #[inline]
    pub fn start(&self) -> usize {
        self.start
    }

    /// Number of samples.
    #[inline]
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// True if the window is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }
}

/// Iterator over overlapping windows.
///
/// Windows overlap because a basecaller has no context at a chunk boundary and
/// its calls there are unreliable; downstream stitching keeps the confident
/// middle of each window and discards the seams. A final short window is
/// emitted rather than dropped, so no signal is silently lost.
pub struct Windows<'a> {
    samples: &'a [f32],
    size: usize,
    stride: usize,
    pos: usize,
    done: bool,
}

impl<'a> Windows<'a> {
    pub(crate) fn new(samples: &'a [f32], size: usize, overlap: usize) -> Windows<'a> {
        // An overlap at or above the window size would make no forward
        // progress; clamp it so the iterator always terminates.
        let stride = size.saturating_sub(overlap).max(1);
        Windows {
            samples,
            size,
            stride,
            pos: 0,
            done: size == 0 || samples.is_empty(),
        }
    }
}

impl<'a> Iterator for Windows<'a> {
    type Item = Window<'a>;

    fn next(&mut self) -> Option<Window<'a>> {
        if self.done {
            return None;
        }
        let start = self.pos;
        let end = (start + self.size).min(self.samples.len());
        // Reaching the end of the trace ends the iteration, after emitting the
        // partial window we are standing on.
        if end >= self.samples.len() {
            self.done = true;
        }
        self.pos += self.stride;
        Some(Window {
            samples: &self.samples[start..end],
            start,
        })
    }
}
