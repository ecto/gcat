//! Streaming FASTA/FASTQ readers and writers.
//!
//! Records stream one at a time and reuse their buffers, so a 100 GB FASTQ
//! costs the same memory as a small one. Sequences stay as raw bytes here —
//! parsing into [`gcat::Seq`] is the caller's choice, because plenty of work
//! (grepping, subsetting, counting) never needs it.

pub mod annot;

use gcat::{Base, Seq, SeqError};
use std::io::{self, BufRead, Write};

/// A sequence record, borrowed from the reader's internal buffers.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct Record<'a> {
    /// The header line, without the leading `>` or `@`.
    pub id: &'a str,
    /// Raw sequence bytes as they appeared in the file.
    pub seq: &'a [u8],
    /// Phred quality bytes, present for FASTQ and absent for FASTA.
    pub qual: Option<&'a [u8]>,
}

impl Record<'_> {
    /// The identifier up to the first whitespace.
    pub fn name(&self) -> &str {
        self.id.split_whitespace().next().unwrap_or("")
    }

    /// Sequence length in bases.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// True if the record has no sequence.
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Parse into a packed [`Seq`], rejecting ambiguity codes.
    pub fn to_seq(&self) -> Result<Seq, SeqError> {
        Seq::from_ascii(self.seq)
    }

    /// Parse into a packed [`Seq`], substituting `fill` for ambiguity codes.
    pub fn to_seq_lossy(&self, fill: Base) -> Result<(Seq, usize), SeqError> {
        Seq::from_ascii_lossy(self.seq, fill)
    }

    /// Mean Phred quality, or `None` for FASTA records and empty reads.
    ///
    /// Averaged in the probability domain, not the log domain: Phred scores are
    /// log-scaled, and averaging them directly overstates read quality.
    pub fn mean_quality(&self, offset: u8) -> Option<f64> {
        let q = self.qual?;
        if q.is_empty() {
            return None;
        }
        let sum_p: f64 = q
            .iter()
            .map(|&c| {
                let phred = c.saturating_sub(offset) as f64;
                libm_pow10(-phred / 10.0)
            })
            .sum();
        let mean_p = sum_p / q.len() as f64;
        Some(if mean_p > 0.0 {
            -10.0 * mean_p.log10()
        } else {
            0.0
        })
    }
}

fn libm_pow10(x: f64) -> f64 {
    10f64.powf(x)
}

/// The standard Phred+33 offset used by essentially all modern FASTQ.
pub const PHRED33: u8 = 33;

/// What went wrong while reading.
#[derive(Debug)]
pub enum ReadError {
    Io(io::Error),
    /// A record did not start with `>` or `@` where one was expected.
    BadHeader {
        line: usize,
    },
    /// A FASTQ record's quality string length did not match its sequence.
    LengthMismatch {
        line: usize,
        seq: usize,
        qual: usize,
    },
    /// The file ended in the middle of a record.
    Truncated {
        line: usize,
    },
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadError::Io(e) => write!(f, "io error: {e}"),
            ReadError::BadHeader { line } => write!(f, "malformed header at line {line}"),
            ReadError::LengthMismatch { line, seq, qual } => write!(
                f,
                "line {line}: sequence is {seq} bases but quality is {qual}"
            ),
            ReadError::Truncated { line } => write!(f, "file ends mid-record at line {line}"),
        }
    }
}

impl std::error::Error for ReadError {}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> ReadError {
        ReadError::Io(e)
    }
}

/// Which format a reader is parsing.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Format {
    Fasta,
    Fastq,
}

impl Format {
    /// Guess the format from the first byte of a stream.
    pub fn sniff(first: u8) -> Option<Format> {
        match first {
            b'>' => Some(Format::Fasta),
            b'@' => Some(Format::Fastq),
            _ => None,
        }
    }
}

/// A streaming reader for FASTA and FASTQ.
///
/// Detects the format from the first record, so callers do not have to care
/// which they were handed.
pub struct Reader<R: BufRead> {
    inner: R,
    format: Option<Format>,
    id: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    peeked: Option<String>,
    line_no: usize,
    finished: bool,
}

impl<R: BufRead> Reader<R> {
    /// Wrap a buffered reader.
    pub fn new(inner: R) -> Reader<R> {
        Reader {
            inner,
            format: None,
            id: String::new(),
            seq: Vec::new(),
            qual: Vec::new(),
            peeked: None,
            line_no: 0,
            finished: false,
        }
    }

    /// The detected format, available once the first record has been read.
    pub fn format(&self) -> Option<Format> {
        self.format
    }

    /// Read one line, returning `None` at end of input.
    fn read_line(&mut self) -> Result<Option<String>, ReadError> {
        if let Some(line) = self.peeked.take() {
            return Ok(Some(line));
        }
        let mut buf = String::new();
        let n = self.inner.read_line(&mut buf)?;
        if n == 0 {
            return Ok(None);
        }
        self.line_no += 1;
        while buf.ends_with('\n') || buf.ends_with('\r') {
            buf.pop();
        }
        Ok(Some(buf))
    }

    /// Advance to the next record.
    ///
    /// Returns `Ok(None)` at clean end of input. The returned [`Record`] borrows
    /// from the reader, so it must be dropped before the next call.
    pub fn next_record(&mut self) -> Result<Option<Record<'_>>, ReadError> {
        if self.finished {
            return Ok(None);
        }
        // Skip blank lines between records.
        let header = loop {
            match self.read_line()? {
                None => {
                    self.finished = true;
                    return Ok(None);
                }
                Some(l) if l.trim().is_empty() => continue,
                Some(l) => break l,
            }
        };

        let first = header.as_bytes()[0];
        let format = match Format::sniff(first) {
            Some(f) => f,
            None => return Err(ReadError::BadHeader { line: self.line_no }),
        };
        // Lock the format to whatever the first record was.
        let format = *self.format.get_or_insert(format);

        self.id.clear();
        self.id.push_str(&header[1..]);
        self.seq.clear();
        self.qual.clear();

        match format {
            Format::Fasta => {
                // Sequence continues over as many lines as it likes, until the
                // next header or end of file.
                loop {
                    match self.read_line()? {
                        None => {
                            self.finished = true;
                            break;
                        }
                        Some(l) if l.starts_with('>') => {
                            self.peeked = Some(l);
                            break;
                        }
                        Some(l) => self.seq.extend_from_slice(l.trim().as_bytes()),
                    }
                }
            }
            Format::Fastq => {
                let seq_line = self
                    .read_line()?
                    .ok_or(ReadError::Truncated { line: self.line_no })?;
                self.seq.extend_from_slice(seq_line.trim().as_bytes());

                let plus = self
                    .read_line()?
                    .ok_or(ReadError::Truncated { line: self.line_no })?;
                if !plus.starts_with('+') {
                    return Err(ReadError::BadHeader { line: self.line_no });
                }

                let qual_line = self
                    .read_line()?
                    .ok_or(ReadError::Truncated { line: self.line_no })?;
                self.qual.extend_from_slice(qual_line.trim().as_bytes());

                if self.qual.len() != self.seq.len() {
                    return Err(ReadError::LengthMismatch {
                        line: self.line_no,
                        seq: self.seq.len(),
                        qual: self.qual.len(),
                    });
                }
            }
        }

        Ok(Some(Record {
            id: &self.id,
            seq: &self.seq,
            qual: if format == Format::Fastq {
                Some(&self.qual)
            } else {
                None
            },
        }))
    }
}

/// Write a FASTA record, wrapping the sequence at `width` columns.
///
/// A `width` of 0 writes the sequence on a single line.
pub fn write_fasta<W: Write>(w: &mut W, id: &str, seq: &[u8], width: usize) -> io::Result<()> {
    writeln!(w, ">{id}")?;
    if width == 0 || seq.len() <= width {
        w.write_all(seq)?;
        writeln!(w)?;
    } else {
        for chunk in seq.chunks(width) {
            w.write_all(chunk)?;
            writeln!(w)?;
        }
    }
    Ok(())
}

/// Write a FASTQ record.
pub fn write_fastq<W: Write>(w: &mut W, id: &str, seq: &[u8], qual: &[u8]) -> io::Result<()> {
    writeln!(w, "@{id}")?;
    w.write_all(seq)?;
    writeln!(w)?;
    writeln!(w, "+")?;
    w.write_all(qual)?;
    writeln!(w)?;
    Ok(())
}
