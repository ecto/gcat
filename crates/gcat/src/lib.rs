//! Core genomic types — alphabets, 2-bit packed sequences, k-mers.
//!
//! `gcat` is `#![no_std]` with `alloc`. Everything here is plain data: no I/O,
//! no numerics, no allocation beyond what a sequence needs. The heavier crates
//! (`gcat-io`, `gcat-signal`) build on top.
//!
//! ```
//! use gcat::{Base, Seq};
//!
//! let s: Seq = "ACGTACGT".parse().unwrap();
//! assert_eq!(s.len(), 8);
//! assert_eq!(s.get(0), Some(Base::A));
//! assert_eq!(s.revcomp().to_string(), "ACGTACGT");
//! ```

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

use alloc::{string::String, vec::Vec};
use core::{fmt, str::FromStr};

mod kmer;
pub use kmer::{Kmer, Kmers};

/// A single DNA base.
///
/// Encoded as 2 bits so that [`Seq`] can pack four bases per byte. Ambiguity
/// codes (`N`, `R`, `Y`, …) are deliberately not representable — see
/// [`SeqError::Ambiguous`] for how they are handled on parse.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
#[repr(u8)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    /// All four bases in encoding order.
    pub const ALL: [Base; 4] = [Base::A, Base::C, Base::G, Base::T];

    /// Decode from the 2-bit representation. Returns `None` for values > 3.
    #[inline]
    pub const fn from_bits(bits: u8) -> Option<Base> {
        match bits {
            0 => Some(Base::A),
            1 => Some(Base::C),
            2 => Some(Base::G),
            3 => Some(Base::T),
            _ => None,
        }
    }

    /// Parse a single ASCII byte. Accepts both cases; `U` maps to `T`.
    #[inline]
    pub const fn from_ascii(b: u8) -> Option<Base> {
        match b {
            b'A' | b'a' => Some(Base::A),
            b'C' | b'c' => Some(Base::C),
            b'G' | b'g' => Some(Base::G),
            b'T' | b't' | b'U' | b'u' => Some(Base::T),
            _ => None,
        }
    }

    /// The uppercase ASCII byte for this base.
    #[inline]
    pub const fn to_ascii(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::C => b'C',
            Base::G => b'G',
            Base::T => b'T',
        }
    }

    /// Watson-Crick complement.
    #[inline]
    pub const fn complement(self) -> Base {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
        }
    }

    /// True for G or C — the bases held by three hydrogen bonds.
    #[inline]
    pub const fn is_gc(self) -> bool {
        matches!(self, Base::C | Base::G)
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Base::A => "A",
            Base::C => "C",
            Base::G => "G",
            Base::T => "T",
        })
    }
}

/// Errors from parsing a [`Seq`].
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum SeqError {
    /// An IUPAC ambiguity code (`N`, `R`, `Y`, …) at the given offset.
    ///
    /// Ambiguity is a real signal, not a parse failure — callers that want to
    /// tolerate it should use [`Seq::from_ascii_lossy`], which substitutes a
    /// concrete base and reports how many positions were affected.
    Ambiguous { offset: usize, byte: u8 },
    /// A byte that is not a base or an ambiguity code at all.
    Invalid { offset: usize, byte: u8 },
}

impl core::error::Error for SeqError {}

impl fmt::Display for SeqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SeqError::Ambiguous { offset, byte } => {
                write!(f, "ambiguity code {:?} at offset {offset}", *byte as char)
            }
            SeqError::Invalid { offset, byte } => {
                write!(f, "invalid byte {:?} at offset {offset}", *byte as char)
            }
        }
    }
}

/// Is this byte a recognised IUPAC ambiguity code?
#[inline]
const fn is_ambiguity(b: u8) -> bool {
    matches!(
        b.to_ascii_uppercase(),
        b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V'
    )
}

/// A 2-bit packed DNA sequence.
///
/// Four bases per byte. A 3 Gb genome costs 750 MB packed rather than 3 GB as
/// ASCII, which is the difference between fitting in memory and not.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct Seq {
    /// Packed bases, 4 per byte, low bits first.
    packed: Vec<u8>,
    /// Base count. The final byte may be partially used.
    len: usize,
}

impl Seq {
    /// An empty sequence.
    #[inline]
    pub const fn new() -> Seq {
        Seq {
            packed: Vec::new(),
            len: 0,
        }
    }

    /// Allocate space for `cap` bases up front.
    pub fn with_capacity(cap: usize) -> Seq {
        Seq {
            packed: Vec::with_capacity(cap.div_ceil(4)),
            len: 0,
        }
    }

    /// Number of bases.
    #[inline]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// True if there are no bases.
    #[inline]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Bytes of packed storage.
    #[inline]
    pub fn packed_bytes(&self) -> usize {
        self.packed.len()
    }

    /// Append one base.
    #[inline]
    pub fn push(&mut self, base: Base) {
        let slot = self.len % 4;
        if slot == 0 {
            self.packed.push(0);
        }
        let byte = self.packed.last_mut().expect("just pushed if needed");
        *byte |= (base as u8) << (slot * 2);
        self.len += 1;
    }

    /// The base at `i`, or `None` if out of bounds.
    #[inline]
    pub fn get(&self, i: usize) -> Option<Base> {
        if i >= self.len {
            return None;
        }
        let bits = (self.packed[i / 4] >> ((i % 4) * 2)) & 0b11;
        Base::from_bits(bits)
    }

    /// Iterate over bases.
    pub fn iter(&self) -> impl Iterator<Item = Base> + '_ {
        (0..self.len).map(move |i| self.get(i).expect("i < len"))
    }

    /// Parse ASCII, rejecting anything that is not an unambiguous base.
    pub fn from_ascii(bytes: &[u8]) -> Result<Seq, SeqError> {
        let mut s = Seq::with_capacity(bytes.len());
        for (offset, &b) in bytes.iter().enumerate() {
            match Base::from_ascii(b) {
                Some(base) => s.push(base),
                None if is_ambiguity(b) => return Err(SeqError::Ambiguous { offset, byte: b }),
                None => return Err(SeqError::Invalid { offset, byte: b }),
            }
        }
        Ok(s)
    }

    /// Parse ASCII, substituting `fill` for ambiguity codes.
    ///
    /// Returns the sequence and the number of positions substituted, so the
    /// caller can decide whether the read is worth keeping. Still errors on
    /// bytes that are not bases or ambiguity codes at all.
    pub fn from_ascii_lossy(bytes: &[u8], fill: Base) -> Result<(Seq, usize), SeqError> {
        let mut s = Seq::with_capacity(bytes.len());
        let mut substituted = 0usize;
        for (offset, &b) in bytes.iter().enumerate() {
            match Base::from_ascii(b) {
                Some(base) => s.push(base),
                None if is_ambiguity(b) => {
                    s.push(fill);
                    substituted += 1;
                }
                None => return Err(SeqError::Invalid { offset, byte: b }),
            }
        }
        Ok((s, substituted))
    }

    /// Reverse complement.
    pub fn revcomp(&self) -> Seq {
        let mut out = Seq::with_capacity(self.len);
        for i in (0..self.len).rev() {
            out.push(self.get(i).expect("i < len").complement());
        }
        out
    }

    /// Fraction of bases that are G or C, or `None` for an empty sequence.
    pub fn gc_content(&self) -> Option<f64> {
        if self.len == 0 {
            return None;
        }
        let gc = self.iter().filter(|b| b.is_gc()).count();
        Some(gc as f64 / self.len as f64)
    }

    /// Iterate over overlapping k-mers. Empty if `k` is 0 or exceeds the length.
    pub fn kmers(&self, k: usize) -> Kmers<'_> {
        Kmers::new(self, k)
    }

    /// Render as uppercase ASCII.
    pub fn to_ascii(&self) -> Vec<u8> {
        self.iter().map(Base::to_ascii).collect()
    }
}

impl FromStr for Seq {
    type Err = SeqError;
    fn from_str(s: &str) -> Result<Seq, SeqError> {
        Seq::from_ascii(s.as_bytes())
    }
}

impl FromIterator<Base> for Seq {
    fn from_iter<I: IntoIterator<Item = Base>>(iter: I) -> Seq {
        let mut s = Seq::new();
        for b in iter {
            s.push(b);
        }
        s
    }
}

impl fmt::Display for Seq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for b in self.iter() {
            fmt::Display::fmt(&b, f)?;
        }
        Ok(())
    }
}

impl fmt::Debug for Seq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Seq({} bp, ", self.len)?;
        // Long sequences elide the middle so debug output stays readable.
        const HEAD: usize = 20;
        if self.len <= HEAD * 2 {
            for b in self.iter() {
                fmt::Display::fmt(&b, f)?;
            }
        } else {
            for i in 0..HEAD {
                fmt::Display::fmt(&self.get(i).expect("i < len"), f)?;
            }
            write!(f, "…")?;
            for i in (self.len - HEAD)..self.len {
                fmt::Display::fmt(&self.get(i).expect("i < len"), f)?;
            }
        }
        f.write_str(")")
    }
}

/// Convenience: collect a [`Seq`] into an owned `String`.
impl From<&Seq> for String {
    fn from(s: &Seq) -> String {
        s.iter().map(|b| b.to_ascii() as char).collect()
    }
}
