//! Canonical k-mers over a packed [`Seq`].

use crate::{Base, Seq};
use core::fmt;

/// A k-mer encoded as 2 bits per base in a `u64`.
///
/// `k` is capped at 32 because that is what fits. The first base of the k-mer
/// occupies the most significant used bits, so lexicographic order on the
/// encoding matches lexicographic order on the sequence.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Kmer {
    bits: u64,
    k: u8,
}

/// The largest `k` representable in a [`Kmer`].
pub const MAX_K: usize = 32;

impl Kmer {
    /// Build from bases, or `None` if `k` is 0 or above [`MAX_K`].
    pub fn from_bases(bases: &[Base]) -> Option<Kmer> {
        let k = bases.len();
        if k == 0 || k > MAX_K {
            return None;
        }
        let mut bits = 0u64;
        for &b in bases {
            bits = (bits << 2) | (b as u64);
        }
        Some(Kmer { bits, k: k as u8 })
    }

    /// The k of this k-mer.
    #[inline]
    pub const fn k(self) -> usize {
        self.k as usize
    }

    /// The raw 2-bit encoding.
    #[inline]
    pub const fn bits(self) -> u64 {
        self.bits
    }

    /// Slide one base in from the right, dropping the leftmost base.
    #[inline]
    pub const fn roll(self, next: Base) -> Kmer {
        // Mask to k*2 bits. k <= 32 so k*2 <= 64; shifting a u64 by 64 is UB in
        // release, so build the mask by shifting a u128 down instead.
        let mask = ((1u128 << (self.k as u32 * 2)) - 1) as u64;
        Kmer {
            bits: ((self.bits << 2) | next as u64) & mask,
            k: self.k,
        }
    }

    /// Decode base `i` of the k-mer, counting from the left.
    #[inline]
    pub const fn base(self, i: usize) -> Option<Base> {
        if i >= self.k as usize {
            return None;
        }
        let shift = (self.k as usize - 1 - i) * 2;
        Base::from_bits(((self.bits >> shift) & 0b11) as u8)
    }

    /// Reverse complement of this k-mer.
    pub fn revcomp(self) -> Kmer {
        let mut bits = 0u64;
        for i in 0..self.k as usize {
            let b = self.base(i).expect("i < k").complement();
            // Writing forward while reading forward reverses the order.
            bits |= (b as u64) << (i * 2);
        }
        Kmer { bits, k: self.k }
    }

    /// The lexicographically smaller of this k-mer and its reverse complement.
    ///
    /// Strand is usually an artifact of which way the molecule went through the
    /// machine, so most counting wants the canonical form.
    pub fn canonical(self) -> Kmer {
        let rc = self.revcomp();
        if rc.bits < self.bits {
            rc
        } else {
            self
        }
    }
}

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.k as usize {
            fmt::Display::fmt(&self.base(i).expect("i < k"), f)?;
        }
        Ok(())
    }
}

impl fmt::Debug for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Kmer({self})")
    }
}

/// Iterator over the overlapping k-mers of a [`Seq`].
pub struct Kmers<'a> {
    seq: &'a Seq,
    k: usize,
    pos: usize,
    current: Option<Kmer>,
}

impl<'a> Kmers<'a> {
    pub(crate) fn new(seq: &'a Seq, k: usize) -> Kmers<'a> {
        Kmers {
            seq,
            k,
            pos: 0,
            current: None,
        }
    }
}

impl Iterator for Kmers<'_> {
    type Item = Kmer;

    fn next(&mut self) -> Option<Kmer> {
        if self.k == 0 || self.k > MAX_K || self.k > self.seq.len() {
            return None;
        }
        match self.current {
            // Seed the first window.
            None => {
                let bases: heapless_vec::Buf = (0..self.k)
                    .map(|i| self.seq.get(i).expect("k <= len"))
                    .collect();
                let km = Kmer::from_bases(bases.as_slice())?;
                self.current = Some(km);
                self.pos = self.k;
                Some(km)
            }
            // Roll forward one base at a time — O(1) per k-mer.
            Some(prev) => {
                if self.pos >= self.seq.len() {
                    return None;
                }
                let next = self.seq.get(self.pos).expect("pos < len");
                let km = prev.roll(next);
                self.current = Some(km);
                self.pos += 1;
                Some(km)
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.k == 0 || self.k > MAX_K || self.k > self.seq.len() {
            return (0, Some(0));
        }
        let total = self.seq.len() - self.k + 1;
        let done = if self.current.is_none() {
            0
        } else {
            self.pos - self.k + 1
        };
        let left = total - done;
        (left, Some(left))
    }
}

/// A fixed-capacity base buffer so seeding a window needs no allocation.
mod heapless_vec {
    use super::{Base, MAX_K};

    pub struct Buf {
        data: [Base; MAX_K],
        len: usize,
    }

    impl Buf {
        pub fn as_slice(&self) -> &[Base] {
            &self.data[..self.len]
        }
    }

    impl FromIterator<Base> for Buf {
        fn from_iter<I: IntoIterator<Item = Base>>(iter: I) -> Buf {
            let mut data = [Base::A; MAX_K];
            let mut len = 0;
            for b in iter.into_iter().take(MAX_K) {
                data[len] = b;
                len += 1;
            }
            Buf { data, len }
        }
    }
}
