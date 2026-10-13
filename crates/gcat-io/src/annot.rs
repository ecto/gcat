//! Annotation: GFF3 features, BED regions, and the gene models built from them.
//!
//! Coordinates are 0-based half-open throughout, whatever the file said. GFF3
//! is 1-based inclusive and BED is already 0-based half-open; converting once
//! at the edge means nothing downstream has to remember which it came from.

use std::collections::HashMap;
use std::io::{self, BufRead};

/// Which strand a feature is on.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    fn parse(s: &str) -> Strand {
        match s {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unknown,
        }
    }
}

/// One annotated interval.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Feature {
    /// The sequence it sits on, as named in the FASTA header (`NC_018741.3`).
    pub seqid: String,
    /// GFF3 type column (`gene`, `mRNA`, `exon`, `CDS`), or `region` for BED.
    pub kind: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub id: Option<String>,
    /// GFF3 `Parent`. Only the first parent is kept: a feature shared by
    /// several transcripts contributes the same interval to each, and the gene
    /// model collapses transcripts anyway.
    pub parent: Option<String>,
    /// `Name`, falling back to `gene`, for display.
    pub name: Option<String>,
}

/// What went wrong while reading annotation.
#[derive(Debug)]
pub enum AnnotError {
    Io(io::Error),
    /// A data line had too few columns or an unparseable coordinate.
    BadLine { line: usize },
}

impl std::fmt::Display for AnnotError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AnnotError::Io(e) => write!(f, "io error: {e}"),
            AnnotError::BadLine { line } => write!(f, "malformed annotation at line {line}"),
        }
    }
}

impl std::error::Error for AnnotError {}

impl From<io::Error> for AnnotError {
    fn from(e: io::Error) -> AnnotError {
        AnnotError::Io(e)
    }
}

/// Undo GFF3's percent-encoding of reserved characters (`%2C`, `%3B`, …).
fn unescape(s: &str) -> String {
    if !s.contains('%') {
        return s.to_string();
    }
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        let hex = |c: u8| (c as char).to_digit(16);
        if b[i] == b'%' && i + 2 < b.len() {
            if let (Some(hi), Some(lo)) = (hex(b[i + 1]), hex(b[i + 2])) {
                out.push((hi * 16 + lo) as u8);
                i += 3;
                continue;
            }
        }
        out.push(b[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// Read GFF3. Stops at a `##FASTA` directive; the sequence after it is the
/// FASTA reader's job.
pub fn read_gff3<R: BufRead>(r: R) -> Result<Vec<Feature>, AnnotError> {
    let mut out = Vec::new();
    for (i, line) in r.lines().enumerate() {
        let line = line?;
        if line.starts_with("##FASTA") {
            break;
        }
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            return Err(AnnotError::BadLine { line: i + 1 });
        }
        let bad = || AnnotError::BadLine { line: i + 1 };
        let start: u64 = cols[3].parse().map_err(|_| bad())?;
        let end: u64 = cols[4].parse().map_err(|_| bad())?;
        let (mut id, mut parent, mut name, mut gene) = (None, None, None, None);
        for kv in cols[8].split(';') {
            let Some((k, v)) = kv.split_once('=') else { continue };
            match k {
                "ID" => id = Some(unescape(v)),
                "Parent" => parent = v.split(',').next().map(unescape),
                "Name" => name = Some(unescape(v)),
                "gene" => gene = Some(unescape(v)),
                _ => {}
            }
        }
        out.push(Feature {
            seqid: cols[0].to_string(),
            kind: cols[2].to_string(),
            start: start.saturating_sub(1),
            end,
            strand: Strand::parse(cols[6]),
            id,
            parent,
            name: name.or(gene),
        });
    }
    Ok(out)
}

/// Read BED: `chrom start end [name [score [strand]]]`. Header lines
/// (`track`, `browser`, `#`) are skipped.
pub fn read_bed<R: BufRead>(r: R) -> Result<Vec<Feature>, AnnotError> {
    let mut out = Vec::new();
    for (i, line) in r.lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") {
            continue;
        }
        let cols: Vec<&str> = t.split_whitespace().collect();
        let bad = || AnnotError::BadLine { line: i + 1 };
        if cols.len() < 3 {
            return Err(bad());
        }
        out.push(Feature {
            seqid: cols[0].to_string(),
            kind: "region".to_string(),
            start: cols[1].parse().map_err(|_| bad())?,
            end: cols[2].parse().map_err(|_| bad())?,
            strand: cols.get(5).map_or(Strand::Unknown, |s| Strand::parse(s)),
            id: None,
            parent: None,
            name: cols.get(3).map(|s| s.to_string()),
        });
    }
    Ok(out)
}

/// A gene with its transcripts collapsed into one row: every exon any
/// transcript uses, merged where they overlap, and the span of coding sequence
/// across all of them. That loses which exons go together, and keeps what a
/// track at gene scale can actually show.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Gene {
    pub seqid: String,
    pub name: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    /// Sorted, non-overlapping.
    pub exons: Vec<(u64, u64)>,
    pub cds: Option<(u64, u64)>,
}

/// Build gene models from GFF3 features by following `Parent` links from
/// exons and CDS up through transcripts to genes.
pub fn genes(features: &[Feature]) -> Vec<Gene> {
    let mut by_id: HashMap<&str, &Feature> = HashMap::new();
    for f in features {
        if let Some(id) = &f.id {
            by_id.insert(id, f);
        }
    }
    let mut order: Vec<&str> = Vec::new();
    let mut models: HashMap<&str, Gene> = HashMap::new();
    for f in features {
        let is_gene = f.parent.is_none() && matches!(f.kind.as_str(), "gene" | "pseudogene");
        if !is_gene {
            continue;
        }
        let Some(id) = f.id.as_deref() else { continue };
        order.push(id);
        models.insert(
            id,
            Gene {
                seqid: f.seqid.clone(),
                name: f.name.clone().unwrap_or_else(|| id.to_string()),
                start: f.start,
                end: f.end,
                strand: f.strand,
                exons: Vec::new(),
                cds: None,
            },
        );
    }
    for f in features {
        let is_exon = f.kind == "exon";
        if !is_exon && f.kind != "CDS" {
            continue;
        }
        let Some(g) = root(&by_id, f).and_then(|r| models.get_mut(r)) else { continue };
        if is_exon {
            g.exons.push((f.start, f.end));
        } else {
            g.cds = Some(match g.cds {
                Some((a, b)) => (a.min(f.start), b.max(f.end)),
                None => (f.start, f.end),
            });
        }
    }
    order
        .into_iter()
        .filter_map(|id| models.remove(id))
        .map(|mut g| {
            g.exons.sort_unstable();
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(g.exons.len());
            for (a, b) in g.exons.drain(..) {
                match merged.last_mut() {
                    Some(last) if a <= last.1 => last.1 = last.1.max(b),
                    _ => merged.push((a, b)),
                }
            }
            g.exons = merged;
            g
        })
        .collect()
}

/// Walk up `Parent` links to the outermost ancestor; for an exon, that is its
/// gene. Bounded so a cycle in a broken file cannot hang the page.
fn root<'a>(by_id: &HashMap<&'a str, &'a Feature>, f: &'a Feature) -> Option<&'a str> {
    let mut cur = f.parent.as_deref()?;
    for _ in 0..8 {
        match by_id.get(cur).and_then(|p| p.parent.as_deref()) {
            Some(up) if by_id.contains_key(up) => cur = up,
            _ => break,
        }
    }
    Some(cur)
}

/// Where a FASTA record sits in a larger sequence, when its header says so.
///
/// NCBI names a subrange fetch `NC_018741.3:110430079-110439152` (1-based,
/// inclusive), and samtools faidx writes the same form. Returns the parent
/// sequence id and the 0-based offset of the record's first base, so
/// annotation in chromosome coordinates lines up with it. A header with no
/// range is its own coordinate system: offset 0.
pub fn header_locus(id: &str) -> (&str, u64) {
    let name = id.split_whitespace().next().unwrap_or("");
    let parsed = name.rsplit_once(':').and_then(|(seq, range)| {
        let (a, b) = range.split_once('-')?;
        let a: u64 = a.replace(',', "").parse().ok()?;
        let b: u64 = b.replace(',', "").parse().ok()?;
        (a >= 1 && b >= a).then_some((seq, a - 1))
    });
    parsed.unwrap_or((name, 0))
}

/// A bedGraph track: a value over each interval, such as read depth.
#[derive(Clone, PartialEq, Debug, Default)]
pub struct Signal {
    /// From the `track name="…"` line, if the file had one.
    pub name: Option<String>,
    /// Values are percentages on a fixed 0–100 scale (identity to another
    /// species) rather than counts scaled to their own maximum (read depth).
    pub percent: bool,
    /// `(seqid, start, end, value)`, in file order.
    pub spans: Vec<(String, u64, u64, f32)>,
}

impl Signal {
    /// The largest value in each of `bins` equal slices of `[start, end)` on
    /// `seqid`, NaN where nothing is recorded — which is not the same as a
    /// recorded zero. Max rather than mean so a peak is
    /// never averaged away; a bin reads zero only when nothing in it is
    /// covered, which is exactly what a deletion looks like.
    pub fn binned(&self, seqid: &str, start: f64, end: f64, bins: usize) -> Vec<f32> {
        let mut out = vec![f32::NAN; bins];
        let width = (end - start) / bins.max(1) as f64;
        if width <= 0.0 {
            return out;
        }
        for (id, a, b, v) in &self.spans {
            if id != seqid || (*b as f64) <= start || (*a as f64) >= end {
                continue;
            }
            let lo = (((*a as f64).max(start) - start) / width).floor() as usize;
            let hi = ((((*b as f64).min(end) - start) / width).ceil() as usize).min(bins);
            for slot in &mut out[lo.min(bins)..hi] {
                // f32::max ignores a NaN operand, so the first value wins.
                *slot = slot.max(*v);
            }
        }
        out
    }
}

impl Signal {
    /// A percentage track from one value per base starting at `offset`,
    /// skipping NaN (no data). Values are rounded to whole percent so runs of
    /// similar identity collapse into one span.
    pub fn from_percent(name: &str, seqid: &str, offset: u64, values: &[f32]) -> Signal {
        let mut spans: Vec<(String, u64, u64, f32)> = Vec::new();
        for (i, v) in values.iter().enumerate() {
            if v.is_nan() {
                continue;
            }
            let (pos, v) = (offset + i as u64, v.round());
            match spans.last_mut() {
                Some(last) if last.2 == pos && last.3 == v => last.2 += 1,
                _ => spans.push((seqid.to_string(), pos, pos + 1, v)),
            }
        }
        Signal { name: Some(name.to_string()), percent: true, spans }
    }
}

/// Read bedGraph: `chrom start end value`, 0-based half-open like BED.
pub fn read_bedgraph<R: BufRead>(r: R) -> Result<Signal, AnnotError> {
    let mut sig = Signal::default();
    for (i, line) in r.lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if let Some(rest) = t.strip_prefix("track") {
            sig.name = rest
                .split_once("name=")
                .map(|(_, n)| n.trim_start_matches('"').split('"').next().unwrap_or("").to_string());
            continue;
        }
        if t.is_empty() || t.starts_with('#') || t.starts_with("browser") {
            continue;
        }
        let cols: Vec<&str> = t.split_whitespace().collect();
        let bad = || AnnotError::BadLine { line: i + 1 };
        if cols.len() < 4 {
            return Err(bad());
        }
        sig.spans.push((
            cols[0].to_string(),
            cols[1].parse().map_err(|_| bad())?,
            cols[2].parse().map_err(|_| bad())?,
            cols[3].parse().map_err(|_| bad())?,
        ));
    }
    Ok(sig)
}
