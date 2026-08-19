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

    /// Retained sequence for the editor, uppercased ASCII.
    seq: Vec<u8>,
    /// Retained Phred bytes, aligned to `seq` while both are within cap.
    qual: Vec<u8>,
    seq_truncated: bool,
    /// Scratch the editor reads through, so rows cross the boundary as a
    /// pointer rather than as a fresh string per scroll frame.
    scratch: Vec<u8>,
}

/// How much sequence to keep for the editor.
///
/// Stored as raw ASCII rather than 2-bit packed, because the editor must show
/// `N` where the genome says `N` — a packed alphabet would have to render
/// ambiguity as some real base, and a sequence viewer that quietly lies about
/// its input is worse than one that stops early. At this cap a bacterial genome
/// fits many times over and the truncation is surfaced in the UI.
const SEQ_CAP: usize = 32_000_000;

/// Quality costs a byte per base with no packing, so it is capped lower.
const QUAL_CAP: usize = 8_000_000;

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
            seq: Vec::new(),
            qual: Vec::new(),
            seq_truncated: false,
            scratch: vec![0u8; 4096],
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
        if self.seq.len() < SEQ_CAP {
            let room = SEQ_CAP - self.seq.len();
            let take = room.min(line.len());
            self.seq
                .extend(line[..take].iter().map(|b| b.to_ascii_uppercase()));
            if take < line.len() {
                self.seq_truncated = true;
            }
        } else {
            self.seq_truncated = true;
        }
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
        if self.qual.len() < QUAL_CAP {
            let room = QUAL_CAP - self.qual.len();
            let take = room.min(line.len());
            self.qual.extend_from_slice(&line[..take]);
        }
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

#[wasm_bindgen]
impl Analyzer {
    /// Bases retained for the editor.
    pub fn seq_len(&self) -> f64 {
        self.seq.len() as f64
    }

    /// True when the sequence outgrew the retention cap.
    pub fn seq_truncated(&self) -> bool {
        self.seq_truncated
    }

    /// True when Phred bytes are available and aligned with the sequence.
    pub fn has_quality(&self) -> bool {
        !self.qual.is_empty()
    }

    /// Copy `len` bases from `start` into the scratch buffer.
    ///
    /// Returns how many bytes were written, which is short at the end of the
    /// sequence. Read them through [`Analyzer::scratch_ptr`].
    pub fn read_bases(&mut self, start: usize, len: usize) -> usize {
        Self::copy_out(&self.seq, &mut self.scratch, start, len)
    }

    /// Copy `len` Phred bytes from `start` into the scratch buffer.
    pub fn read_quality(&mut self, start: usize, len: usize) -> usize {
        Self::copy_out(&self.qual, &mut self.scratch, start, len)
    }

    fn copy_out(src: &[u8], scratch: &mut Vec<u8>, start: usize, len: usize) -> usize {
        if start >= src.len() || len == 0 {
            return 0;
        }
        let end = (start + len).min(src.len());
        let n = end - start;
        if scratch.len() < n {
            scratch.resize(n, 0);
        }
        scratch[..n].copy_from_slice(&src[start..end]);
        n
    }

    /// Pointer to the scratch buffer.
    pub fn scratch_ptr(&self) -> *const u8 {
        self.scratch.as_ptr()
    }

    /// GC fraction over a span, or -1 when the span holds no unambiguous base.
    pub fn gc_of(&self, start: usize, len: usize) -> f32 {
        if start >= self.seq.len() {
            return -1.0;
        }
        let end = (start + len).min(self.seq.len());
        let mut gc = 0u32;
        let mut known = 0u32;
        for &b in &self.seq[start..end] {
            match b {
                b'G' | b'C' => {
                    gc += 1;
                    known += 1;
                }
                b'A' | b'T' | b'U' => known += 1,
                _ => {}
            }
        }
        if known == 0 {
            -1.0
        } else {
            gc as f32 / known as f32
        }
    }

    /// GC fraction in each of `bins` equal spans, for the minimap.
    pub fn gc_profile(&self, bins: usize) -> Vec<f32> {
        let bins = bins.max(1);
        let mut out = Vec::with_capacity(bins);
        if self.seq.is_empty() {
            return out;
        }
        let step = self.seq.len() as f64 / bins as f64;
        for i in 0..bins {
            let start = (i as f64 * step) as usize;
            let end = (((i + 1) as f64 * step) as usize).min(self.seq.len());
            out.push(self.gc_of(start, end.saturating_sub(start)));
        }
        out
    }

    /// Offset of the next occurrence of `motif` at or after `from`, or -1.
    ///
    /// Case-insensitive. `N` in the motif matches any base, so a caller can
    /// search a degenerate pattern without a regex engine.
    pub fn find_motif(&mut self, motif: &str, from: usize) -> f64 {
        let needle: Vec<u8> = motif.bytes().map(|b| b.to_ascii_uppercase()).collect();
        if needle.is_empty() || needle.len() > self.seq.len() {
            return -1.0;
        }
        let last = self.seq.len() - needle.len();
        for i in from.min(last + 1)..=last {
            if Self::matches(&self.seq[i..i + needle.len()], &needle) {
                return i as f64;
            }
        }
        -1.0
    }

    /// Count every occurrence of `motif`, overlaps included.
    pub fn count_motif(&self, motif: &str) -> f64 {
        let needle: Vec<u8> = motif.bytes().map(|b| b.to_ascii_uppercase()).collect();
        if needle.is_empty() || needle.len() > self.seq.len() {
            return 0.0;
        }
        let mut n = 0f64;
        for i in 0..=(self.seq.len() - needle.len()) {
            if Self::matches(&self.seq[i..i + needle.len()], &needle) {
                n += 1.0;
            }
        }
        n
    }

    fn matches(window: &[u8], needle: &[u8]) -> bool {
        window
            .iter()
            .zip(needle)
            .all(|(&a, &b)| b == b'N' || a == b || (b == b'T' && a == b'U'))
    }

    /// Reverse complement of a span, as a string for the selection readout.
    pub fn revcomp_of(&self, start: usize, len: usize) -> String {
        if start >= self.seq.len() {
            return String::new();
        }
        let end = (start + len).min(self.seq.len());
        self.seq[start..end]
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => 'T',
                b'T' | b'U' => 'A',
                b'C' => 'G',
                b'G' => 'C',
                _ => 'N',
            })
            .collect()
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
