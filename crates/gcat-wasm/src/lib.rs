//! Browser bindings — streaming analysis entirely client-side.
//!
//! Nothing here uploads anything. A file is read in chunks by JavaScript and
//! pushed through [`Analyzer::push`]; the accumulators live in WebAssembly
//! memory and the page reads them back as typed arrays over that same memory,
//! so a megapixel image crosses the boundary as a pointer rather than as JSON.
//!
//! The analyzer is a state machine over the byte stream rather than a parser
//! over a whole file, because the point is to handle a FASTQ larger than the
//! tab's memory and to show results while it is still reading.

use gcat::Base;
use gcat_viz::{colorize, Cgr, GcSkew, Hist2d, Palette};
use wasm_bindgen::prelude::*;

/// Which line of a FASTQ record the stream is currently inside.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Section {
    Header,
    Sequence,
    Plus,
    Quality,
}

/// What the file turned out to be.
#[derive(Copy, Clone, PartialEq, Eq)]
enum Kind {
    Unknown,
    Fasta,
    Fastq,
}

/// Streaming analyzer over FASTA or FASTQ bytes.
#[wasm_bindgen]
pub struct Analyzer {
    kind: Kind,
    section: Section,
    line: Vec<u8>,

    cgr: Cgr,
    skew: GcSkew,
    hist: Hist2d,
    rgba: Vec<u8>,

    records: u64,
    bases: u64,
    gc: u64,
    ambiguous: u64,
    /// Length of the record currently being read.
    current_len: usize,
    /// Sum of error probabilities for the current record's quality string.
    current_err: f64,
    /// Read lengths, capped so a huge file cannot exhaust memory.
    lengths: Vec<u32>,
    truncated: bool,
    qual_sum: f64,
    qual_n: u64,
}

/// Beyond this many reads, stop retaining individual lengths. Summary
/// statistics stay exact via the running counters; only the exact N50 and
/// median degrade, and they degrade to an estimate over a large sample rather
/// than to nonsense.
const MAX_RETAINED_LENGTHS: usize = 2_000_000;

#[wasm_bindgen]
impl Analyzer {
    /// `cgr_order` gives a `2^order` square; `skew_window` is in bases.
    #[wasm_bindgen(constructor)]
    pub fn new(cgr_order: usize, skew_window: usize) -> Analyzer {
        let cgr = Cgr::new(cgr_order);
        let cells = cgr.size() * cgr.size();
        let hist = Hist2d::read_length_quality(160, 120);
        let pixels = cells.max(hist.cols() * hist.rows());
        Analyzer {
            kind: Kind::Unknown,
            section: Section::Header,
            line: Vec::with_capacity(256),
            cgr,
            skew: GcSkew::new(skew_window),
            hist,
            rgba: vec![0u8; pixels * 4],
            records: 0,
            bases: 0,
            gc: 0,
            ambiguous: 0,
            current_len: 0,
            current_err: 0.0,
            lengths: Vec::new(),
            truncated: false,
            qual_sum: 0.0,
            qual_n: 0,
        }
    }

    /// Feed the next chunk of file bytes. Any split point is safe — the
    /// analyzer buffers a partial line across calls.
    pub fn push(&mut self, chunk: &[u8]) {
        for &byte in chunk {
            if byte == b'\n' || byte == b'\r' {
                if !self.line.is_empty() {
                    self.consume_line();
                    self.line.clear();
                }
            } else {
                self.line.push(byte);
            }
        }
    }

    /// Flush any trailing partial line and close open windows.
    pub fn finish(&mut self) {
        if !self.line.is_empty() {
            let line = core::mem::take(&mut self.line);
            self.line = line;
            self.consume_line();
            self.line.clear();
        }
        self.close_record();
        self.skew.finish();
    }

    fn consume_line(&mut self) {
        let first = self.line[0];
        if self.kind == Kind::Unknown {
            self.kind = match first {
                b'>' => Kind::Fasta,
                b'@' => Kind::Fastq,
                // Not a recognisable header: treat the whole file as bare
                // sequence rather than refusing, so pasted sequence works.
                _ => Kind::Fasta,
            };
            if first == b'>' || first == b'@' {
                self.section = Section::Sequence;
                self.records += 1;
                return;
            }
        }

        match self.kind {
            Kind::Fasta => {
                if first == b'>' {
                    self.close_record();
                    self.records += 1;
                    self.cgr.discontinuity();
                } else {
                    self.take_sequence();
                }
            }
            Kind::Fastq => match self.section {
                Section::Header => {
                    self.close_record();
                    self.records += 1;
                    self.cgr.discontinuity();
                    self.section = Section::Sequence;
                }
                Section::Sequence => {
                    self.take_sequence();
                    self.section = Section::Plus;
                }
                Section::Plus => self.section = Section::Quality,
                Section::Quality => {
                    self.take_quality();
                    self.section = Section::Header;
                }
            },
            Kind::Unknown => {}
        }
    }

    fn take_sequence(&mut self) {
        let line = core::mem::take(&mut self.line);
        for &byte in &line {
            match Base::from_ascii(byte) {
                Some(b) => {
                    self.cgr.push(b);
                    self.skew.push(b);
                    self.bases += 1;
                    self.current_len += 1;
                    if b.is_gc() {
                        self.gc += 1;
                    }
                }
                None => {
                    // Ambiguity breaks the k-mer window but still counts toward
                    // read length, because the base is really there.
                    self.cgr.discontinuity();
                    self.ambiguous += 1;
                    self.bases += 1;
                    self.current_len += 1;
                }
            }
        }
        self.line = line;
    }

    fn take_quality(&mut self) {
        let line = core::mem::take(&mut self.line);
        for &byte in &line {
            let phred = byte.saturating_sub(33) as f64;
            // Accumulate error probability, not the Phred score: scores are
            // logarithmic and averaging them directly overstates quality.
            self.current_err += pow10(-phred / 10.0);
        }
        self.line = line;
    }

    fn close_record(&mut self) {
        if self.current_len == 0 {
            return;
        }
        if self.lengths.len() < MAX_RETAINED_LENGTHS {
            self.lengths.push(self.current_len as u32);
        } else {
            self.truncated = true;
        }

        if self.kind == Kind::Fastq && self.current_err > 0.0 {
            let mean_p = self.current_err / self.current_len as f64;
            let q = if mean_p > 0.0 {
                -10.0 * log10(mean_p)
            } else {
                60.0
            };
            self.hist.push(self.current_len as f32, q as f32);
            self.qual_sum += q;
            self.qual_n += 1;
        }
        self.current_len = 0;
        self.current_err = 0.0;
    }

    /// Render the chaos game grid into the shared pixel buffer.
    pub fn render_cgr(&mut self, palette: u8, gamma: f32) {
        colorize(
            self.cgr.counts(),
            self.cgr.peak(),
            palette_of(palette),
            gamma,
            &mut self.rgba,
        );
    }

    /// Render the length-versus-quality histogram into the shared pixel buffer.
    pub fn render_hist(&mut self, palette: u8, gamma: f32) {
        colorize(
            self.hist.counts(),
            self.hist.peak(),
            palette_of(palette),
            gamma,
            &mut self.rgba,
        );
    }

    /// Pointer to the RGBA buffer, for a zero-copy `ImageData` on the JS side.
    pub fn pixels_ptr(&self) -> *const u8 {
        self.rgba.as_ptr()
    }

    pub fn cgr_size(&self) -> usize {
        self.cgr.size()
    }

    pub fn cgr_order(&self) -> usize {
        self.cgr.order()
    }

    pub fn cgr_sparsity(&self) -> f32 {
        self.cgr.sparsity()
    }

    pub fn hist_cols(&self) -> usize {
        self.hist.cols()
    }

    pub fn hist_rows(&self) -> usize {
        self.hist.rows()
    }

    pub fn hist_total(&self) -> f64 {
        self.hist.total() as f64
    }

    /// The cumulative GC skew curve, as a flat `[offset, cumulative, ...]` array.
    pub fn skew_curve(&self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.skew.points().len() * 2);
        for p in self.skew.points() {
            out.push(p.offset as f32);
            out.push(p.cumulative);
        }
        out
    }

    /// Predicted replication origin in bases, or -1 when there is no curve.
    pub fn origin(&self) -> f64 {
        self.skew.origin().map_or(-1.0, |v| v as f64)
    }

    /// Predicted replication terminus in bases, or -1 when there is no curve.
    pub fn terminus(&self) -> f64 {
        self.skew.terminus().map_or(-1.0, |v| v as f64)
    }

    pub fn records(&self) -> f64 {
        self.records as f64
    }

    pub fn bases(&self) -> f64 {
        self.bases as f64
    }

    pub fn gc_fraction(&self) -> f64 {
        let known = self.bases - self.ambiguous;
        if known == 0 {
            0.0
        } else {
            self.gc as f64 / known as f64
        }
    }

    pub fn ambiguous(&self) -> f64 {
        self.ambiguous as f64
    }

    pub fn mean_quality(&self) -> f64 {
        if self.qual_n == 0 {
            -1.0
        } else {
            self.qual_sum / self.qual_n as f64
        }
    }

    pub fn is_fastq(&self) -> bool {
        self.kind == Kind::Fastq
    }

    pub fn lengths_truncated(&self) -> bool {
        self.truncated
    }

    /// `[n, min, median, mean, max, n50]` over retained read lengths.
    pub fn length_stats(&self) -> Vec<f64> {
        if self.lengths.is_empty() {
            return vec![0.0; 6];
        }
        let mut sorted = self.lengths.clone();
        sorted.sort_unstable();
        let n = sorted.len();
        let total: u64 = sorted.iter().map(|&l| l as u64).sum();
        let median = if n % 2 == 1 {
            sorted[n / 2] as f64
        } else {
            (sorted[n / 2 - 1] as f64 + sorted[n / 2] as f64) / 2.0
        };

        let half = total / 2;
        let mut acc = 0u64;
        let mut n50 = 0f64;
        for &l in sorted.iter().rev() {
            acc += l as u64;
            if acc >= half {
                n50 = l as f64;
                break;
            }
        }

        vec![
            n as f64,
            sorted[0] as f64,
            median,
            total as f64 / n as f64,
            sorted[n - 1] as f64,
            n50,
        ]
    }
}

fn palette_of(v: u8) -> Palette {
    match v {
        1 => Palette::Ice,
        2 => Palette::Mono,
        _ => Palette::Ember,
    }
}

fn pow10(x: f64) -> f64 {
    libm_exp(x * core::f64::consts::LN_10)
}

fn log10(x: f64) -> f64 {
    libm_ln(x) / core::f64::consts::LN_10
}

fn libm_exp(x: f64) -> f64 {
    if x < -700.0 {
        return 0.0;
    }
    let n = (x / core::f64::consts::LN_2).round();
    let r = x - n * core::f64::consts::LN_2;
    let mut term = 1.0f64;
    let mut sum = 1.0f64;
    for i in 1..18 {
        term *= r / i as f64;
        sum += term;
    }
    sum * pow2(n as i32)
}

fn libm_ln(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let bits = x.to_bits();
    let exp = (((bits >> 52) & 0x7ff) as i32) - 1023;
    let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    let z = (mantissa - 1.0) / (mantissa + 1.0);
    let z2 = z * z;
    let mut series = 0.0f64;
    let mut term = z;
    let mut k = 1.0f64;
    for _ in 0..20 {
        series += term / k;
        term *= z2;
        k += 2.0;
    }
    core::f64::consts::LN_2 * exp as f64 + 2.0 * series
}

fn pow2(n: i32) -> f64 {
    let clamped = n.clamp(-1022, 1023);
    f64::from_bits(((clamped + 1023) as u64) << 52)
}
