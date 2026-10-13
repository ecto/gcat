//! Pairwise alignment for comparing a region across species.
//!
//! The question it answers is narrow: for each base of a reference, did a
//! second sequence keep it, change it, or lose it? That is what a
//! conservation track draws, and it does not need a general aligner.
//!
//! Seed and chain, then fill:
//!
//! 1. **Seed.** Every 12-mer of the query that occurs at most a few times in
//!    the reference is an anchor. Repeats are dropped, since an anchor that
//!    could sit in ten places says nothing about where it belongs.
//! 2. **Chain.** Anchors on the same diagonal merge into runs, and the best
//!    colinear chain of runs is kept: the backbone of the alignment.
//! 3. **Fill.** Between consecutive runs, a small global alignment decides
//!    base by base what matched, mismatched or was lost.
//!
//! Both strands of the query are tried and the better chain wins, so a
//! sequence fetched in the other orientation still lines up.
//!
//! `#![no_std]` with `alloc`, so it runs in the browser on files the user drops.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

use alloc::vec;
use alloc::vec::Vec;

const K: usize = 12;
/// A seed occurring more often than this in the reference is a repeat.
const MAX_OCC: usize = 4;
/// Gaps larger than this on either side are left unaligned rather than
/// filled. Global alignment of two unrelated stretches still lines up about
/// half their bases by chance, so filling a long unanchored gap invents
/// conservation; on cat against mouse, raising this past 500 aligns kilobases
/// that no seed supports.
const MAX_FILL: usize = 500;
/// How many earlier runs a run may chain onto.
const LOOKBACK: usize = 400;

/// What happened to one reference base.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Fate {
    /// The query has the same base here.
    Match,
    /// The query has a different base here.
    Mismatch,
    /// The query skips this base: a deletion relative to the reference.
    Gap,
    /// Outside the aligned stretch, or between runs too far apart to fill.
    Unaligned,
}

/// The fate of every reference base, plus which strand of the query matched.
#[derive(Clone, Debug)]
pub struct Alignment {
    pub fates: Vec<Fate>,
    /// True when the query aligned as its reverse complement.
    pub reverse: bool,
}

impl Alignment {
    /// Percent identity in a centred window around each base: matches over
    /// aligned bases, gaps counting as differences. `None` where fewer than
    /// half the window's bases are aligned at all, so the track shows nothing
    /// rather than a confident number from three bases.
    pub fn identity(&self, window: usize) -> Vec<Option<f32>> {
        let n = self.fates.len();
        let half = window.max(1) / 2;
        // Prefix sums of matches and of aligned bases.
        let mut m = vec![0u32; n + 1];
        let mut a = vec![0u32; n + 1];
        for (i, f) in self.fates.iter().enumerate() {
            m[i + 1] = m[i] + (*f == Fate::Match) as u32;
            a[i + 1] = a[i] + (*f != Fate::Unaligned) as u32;
        }
        (0..n)
            .map(|i| {
                let lo = i.saturating_sub(half);
                let hi = (i + half + 1).min(n);
                let aligned = a[hi] - a[lo];
                if self.fates[i] == Fate::Unaligned || (aligned as usize) * 2 < hi - lo {
                    return None;
                }
                Some((m[hi] - m[lo]) as f32 / aligned as f32 * 100.0)
            })
            .collect()
    }
}

fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

// Every sort below is over plain `u64`s, with the fields packed high to low in
// sort order. One element type means one copy of the sort routine in the
// browser bundle instead of one per tuple; it is also faster.

/// `kmer << 32 | position` for every clean k-mer, sorted. A 12-mer is 24 bits.
fn kmers(s: &[u8]) -> Vec<u64> {
    let mask = (1u64 << (2 * K)) - 1;
    let mut out = Vec::with_capacity(s.len());
    let (mut v, mut run) = (0u64, 0usize);
    for (i, &b) in s.iter().enumerate() {
        match code(b) {
            Some(c) => {
                v = ((v << 2) | c) & mask;
                run += 1;
                if run >= K {
                    out.push(v << 32 | (i + 1 - K) as u64);
                }
            }
            None => run = 0,
        }
    }
    out.sort_unstable();
    out
}

/// A stretch where reference and query agree exactly, on one diagonal.
#[derive(Copy, Clone, Debug)]
struct Run {
    r: usize,
    q: usize,
    len: usize,
}

fn runs(reference: &[u8], query: &[u8]) -> Vec<Run> {
    const LOW: u64 = 0xffff_ffff;
    let index = kmers(reference);
    // Diagonals run from -reference.len() to query.len(); shift them positive.
    let shift = reference.len() as u64;
    let mut anchors: Vec<u64> = Vec::new();
    for entry in kmers(query) {
        let (km, qp) = (entry >> 32, entry & LOW);
        let lo = index.partition_point(|&e| e >> 32 < km);
        let hi = index.partition_point(|&e| e >> 32 <= km);
        if hi > lo && hi - lo <= MAX_OCC {
            anchors.extend(index[lo..hi].iter().map(|&e| (qp + shift - (e & LOW)) << 32 | (e & LOW)));
        }
    }
    // Sorted by diagonal, then position, so overlapping anchors on one
    // diagonal are adjacent and merge in a single pass.
    anchors.sort_unstable();
    let mut out: Vec<Run> = Vec::new();
    for a in anchors {
        let r = (a & LOW) as usize;
        let q = ((a >> 32) + r as u64 - shift) as usize;
        if let Some(last) = out.last_mut() {
            if q as i64 - r as i64 == last.q as i64 - last.r as i64 && r <= last.r + last.len {
                last.len = last.len.max(r + K - last.r);
                continue;
            }
        }
        out.push(Run { r, q, len: K });
    }
    let mut order: Vec<u64> = out.iter().enumerate().map(|(i, run)| (run.r as u64) << 32 | i as u64).collect();
    order.sort_unstable();
    order.into_iter().map(|k| out[(k & LOW) as usize]).collect()
}

/// The highest-scoring colinear chain: each run starts after the previous one
/// ends in both sequences. Score is bases covered, less a small charge for
/// changing diagonal so the chain prefers one consistent alignment.
///
/// Runs may overlap their predecessor by a few bases, where sequence beside a
/// deletion matches by chance; the later run is trimmed rather than refused.
fn chain(runs: &[Run]) -> Vec<Run> {
    if runs.is_empty() {
        return Vec::new();
    }
    let mut score = vec![0f64; runs.len()];
    let mut prev = vec![usize::MAX; runs.len()];
    let mut trim = vec![0usize; runs.len()];
    for i in 0..runs.len() {
        score[i] = runs[i].len as f64;
        let ri = runs[i];
        for j in (i.saturating_sub(LOOKBACK)..i).rev() {
            let rj = runs[j];
            let over = (rj.r + rj.len).saturating_sub(ri.r).max((rj.q + rj.len).saturating_sub(ri.q));
            if over >= ri.len || rj.r >= ri.r || rj.q >= ri.q {
                continue;
            }
            let dr = ri.r + over - (rj.r + rj.len);
            let dq = ri.q + over - (rj.q + rj.len);
            let cost = 0.05 * dr.abs_diff(dq) as f64 + 0.01 * dr.min(dq) as f64;
            let s = score[j] + (ri.len - over) as f64 - cost;
            if s > score[i] {
                score[i] = s;
                prev[i] = j;
                trim[i] = over;
            }
        }
    }
    let mut best = 0;
    for i in 1..runs.len() {
        if score[i] > score[best] {
            best = i;
        }
    }
    let mut out = Vec::new();
    let mut i = best;
    loop {
        let t = trim[i];
        out.push(Run { r: runs[i].r + t, q: runs[i].q + t, len: runs[i].len - t });
        if prev[i] == usize::MAX {
            break;
        }
        i = prev[i];
    }
    out.reverse();
    out
}

/// Global alignment of a reference gap against a query gap, writing each
/// reference base's fate. Linear gap cost; small enough cells that a byte of
/// traceback each is fine.
fn fill(r: &[u8], q: &[u8], out: &mut [Fate]) {
    let (n, m) = (r.len(), q.len());
    if n == 0 {
        return;
    }
    if m == 0 {
        out.fill(Fate::Gap);
        return;
    }
    const MATCH: i32 = 1;
    const MISMATCH: i32 = -1;
    const GAP: i32 = -2;
    let w = m + 1;
    let mut dir = vec![0u8; (n + 1) * w]; // 0 diag, 1 up (ref gap→query skips), 2 left
    let mut row: Vec<i32> = (0..=m as i32).map(|j| j * GAP).collect();
    dir[1..=m].fill(2);
    for i in 1..=n {
        let mut diag = row[0];
        row[0] = i as i32 * GAP;
        dir[i * w] = 1;
        for j in 1..=m {
            let same = r[i - 1].eq_ignore_ascii_case(&q[j - 1]);
            let d = diag + if same { MATCH } else { MISMATCH };
            let u = row[j] + GAP;
            let l = row[j - 1] + GAP;
            diag = row[j];
            let (best, how) = if d >= u && d >= l { (d, 0) } else if u >= l { (u, 1) } else { (l, 2) };
            row[j] = best;
            dir[i * w + j] = how;
        }
    }
    let (mut i, mut j) = (n, m);
    while i > 0 {
        match dir[i * w + j] {
            0 if j > 0 => {
                out[i - 1] = if r[i - 1].eq_ignore_ascii_case(&q[j - 1]) { Fate::Match } else { Fate::Mismatch };
                i -= 1;
                j -= 1;
            }
            2 if j > 0 => j -= 1,
            _ => {
                out[i - 1] = Fate::Gap;
                i -= 1;
            }
        }
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn align_strand(reference: &[u8], query: &[u8]) -> (Vec<Fate>, usize) {
    let backbone = chain(&runs(reference, query));
    let mut fates = vec![Fate::Unaligned; reference.len()];
    let mut covered = 0;
    for (i, run) in backbone.iter().enumerate() {
        for f in &mut fates[run.r..run.r + run.len] {
            *f = Fate::Match;
        }
        covered += run.len;
        let Some(next) = backbone.get(i + 1) else { break };
        let (r0, r1) = (run.r + run.len, next.r);
        let (q0, q1) = (run.q + run.len, next.q);
        if r1 - r0 <= MAX_FILL && q1 - q0 <= MAX_FILL {
            fill(&reference[r0..r1], &query[q0..q1], &mut fates[r0..r1]);
        }
    }
    (fates, covered)
}

/// Align `query` onto `reference`, trying both strands.
pub fn align(reference: &[u8], query: &[u8]) -> Alignment {
    let (fwd, fwd_cov) = align_strand(reference, query);
    let rc = revcomp(query);
    let (rev, rev_cov) = align_strand(reference, &rc);
    if rev_cov > fwd_cov {
        Alignment { fates: rev, reverse: true }
    } else {
        Alignment { fates: fwd, reverse: false }
    }
}

