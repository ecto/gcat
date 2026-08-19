//! Subcommand implementations.

use crate::{echo, inputs, Fallible};
use gcat::Base;
use gcat_io::{Reader, PHRED33};
use std::collections::HashMap;
use std::io::Write;

const COMMANDS: &[&str] = &["cat", "stats", "head", "grep", "rc", "kmer"];

/// Does this argument name a subcommand?
pub fn is_command(s: &str) -> bool {
    COMMANDS.contains(&s)
}

/// Split arguments into flags (`-x`, `--long`, with values) and file paths.
///
/// `-` is a file (stdin), not a flag.
pub fn split_args(args: &[String]) -> (HashMap<String, String>, Vec<String>) {
    let mut flags = HashMap::new();
    let mut files = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if a == "-" || !a.starts_with('-') {
            files.push(a.clone());
            i += 1;
            continue;
        }
        let key = a.trim_start_matches('-').to_string();
        // A flag consumes the next argument as its value unless that argument
        // is itself a flag or a file that follows another flag.
        if let Some((k, v)) = key.split_once('=') {
            flags.insert(k.to_string(), v.to_string());
            i += 1;
        } else if i + 1 < args.len() && !args[i + 1].starts_with('-') {
            flags.insert(key, args[i + 1].clone());
            i += 2;
        } else {
            flags.insert(key, "true".to_string());
            i += 1;
        }
    }
    (flags, files)
}

fn flag_usize(flags: &HashMap<String, String>, names: &[&str], default: usize) -> usize {
    for n in names {
        if let Some(v) = flags.get(*n) {
            if let Ok(x) = v.parse() {
                return x;
            }
        }
    }
    default
}

const FASTA_WIDTH: usize = 60;

/// Concatenate records to stdout.
pub fn cat<W: Write>(files: &[String], out: &mut W) -> Fallible {
    for (_, r) in inputs(files) {
        let mut reader = Reader::new(r);
        while let Some(rec) = reader.next_record()? {
            echo(out, rec.id, rec.seq, rec.qual, FASTA_WIDTH)?;
        }
    }
    Ok(())
}

/// First N records.
pub fn head<W: Write>(files: &[String], flags: &HashMap<String, String>, out: &mut W) -> Fallible {
    let n = flag_usize(flags, &["n", "num"], 10);
    let mut seen = 0;
    for (_, r) in inputs(files) {
        let mut reader = Reader::new(r);
        while let Some(rec) = reader.next_record()? {
            if seen >= n {
                return Ok(());
            }
            echo(out, rec.id, rec.seq, rec.qual, FASTA_WIDTH)?;
            seen += 1;
        }
    }
    Ok(())
}

/// Records whose id or sequence contains a pattern.
pub fn grep<W: Write>(files: &[String], flags: &HashMap<String, String>, out: &mut W) -> Fallible {
    let pattern = flags
        .get("p")
        .or_else(|| flags.get("pattern"))
        .or_else(|| flags.get("e"))
        .ok_or("grep needs a pattern: gcat grep -p ACGT reads.fq")?;
    let ids_only = flags.contains_key("id");
    let pat_bytes = pattern.as_bytes();

    for (_, r) in inputs(files) {
        let mut reader = Reader::new(r);
        while let Some(rec) = reader.next_record()? {
            let hit = if ids_only {
                rec.id.contains(pattern.as_str())
            } else {
                rec.id.contains(pattern.as_str()) || contains(rec.seq, pat_bytes)
            };
            if hit {
                echo(out, rec.id, rec.seq, rec.qual, FASTA_WIDTH)?;
            }
        }
    }
    Ok(())
}

/// Naive substring search. Fine for the read lengths this sees; a real index
/// belongs in a different crate, not in the `cat` equivalent.
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() || needle.len() > haystack.len() {
        return false;
    }
    haystack
        .windows(needle.len())
        .any(|w| w.eq_ignore_ascii_case(needle))
}

/// Reverse complement every record.
pub fn rc<W: Write>(files: &[String], out: &mut W) -> Fallible {
    for (_, r) in inputs(files) {
        let mut reader = Reader::new(r);
        while let Some(rec) = reader.next_record()? {
            let (seq, _) = rec.to_seq_lossy(Base::A)?;
            let rc = seq.revcomp();
            // Quality is per-base, so it reverses with the sequence.
            let qual: Option<Vec<u8>> = rec.qual.map(|q| q.iter().rev().copied().collect());
            echo(out, rec.id, &rc.to_ascii(), qual.as_deref(), FASTA_WIDTH)?;
        }
    }
    Ok(())
}

/// Count canonical k-mers across all records.
pub fn kmer<W: Write>(files: &[String], flags: &HashMap<String, String>, out: &mut W) -> Fallible {
    let k = flag_usize(flags, &["k"], 21);
    if k == 0 || k > 32 {
        return Err(format!("k must be between 1 and 32, got {k}").into());
    }
    let top = flag_usize(flags, &["top", "n"], 20);

    let mut counts: HashMap<u64, u32> = HashMap::new();
    let mut skipped = 0usize;
    for (_, r) in inputs(files) {
        let mut reader = Reader::new(r);
        while let Some(rec) = reader.next_record()? {
            // Reads with ambiguity codes are skipped rather than filled: a
            // substituted base would create k-mers that were never observed.
            let seq = match rec.to_seq() {
                Ok(s) => s,
                Err(_) => {
                    skipped += 1;
                    continue;
                }
            };
            for km in seq.kmers(k) {
                *counts.entry(km.canonical().bits()).or_insert(0) += 1;
            }
        }
    }

    let mut ranked: Vec<(u64, u32)> = counts.into_iter().collect();
    // Sort by count descending, then by encoding so output is deterministic.
    ranked.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    writeln!(out, "# k={k}  distinct={}", ranked.len())?;
    if skipped > 0 {
        writeln!(
            out,
            "# skipped {skipped} record(s) containing ambiguity codes"
        )?;
    }
    for (bits, count) in ranked.iter().take(top) {
        writeln!(out, "{}\t{}", decode(*bits, k), count)?;
    }
    Ok(())
}

/// Render a 2-bit k-mer encoding back to bases.
fn decode(bits: u64, k: usize) -> String {
    (0..k)
        .map(|i| {
            let shift = (k - 1 - i) * 2;
            let b = Base::from_bits(((bits >> shift) & 0b11) as u8).expect("2 bits");
            b.to_ascii() as char
        })
        .collect()
}

/// Summarise records: count, bases, N50, GC.
pub fn stats<W: Write>(files: &[String], out: &mut W) -> Fallible {
    let mut lengths: Vec<usize> = Vec::new();
    let mut total_bases = 0usize;
    let mut gc = 0usize;
    let mut acgt = 0usize;
    let mut ambiguous = 0usize;
    let mut qual_sum = 0f64;
    let mut qual_n = 0usize;

    for (_, r) in inputs(files) {
        let mut reader = Reader::new(r);
        while let Some(rec) = reader.next_record()? {
            lengths.push(rec.len());
            total_bases += rec.len();
            for &b in rec.seq {
                match Base::from_ascii(b) {
                    Some(base) => {
                        acgt += 1;
                        if base.is_gc() {
                            gc += 1;
                        }
                    }
                    None => ambiguous += 1,
                }
            }
            if let Some(q) = rec.mean_quality(PHRED33) {
                qual_sum += q;
                qual_n += 1;
            }
        }
    }

    if lengths.is_empty() {
        writeln!(out, "no records")?;
        return Ok(());
    }

    lengths.sort_unstable();
    let n = lengths.len();
    let min = lengths[0];
    let max = lengths[n - 1];
    let mean = total_bases as f64 / n as f64;
    let median = if n % 2 == 1 {
        lengths[n / 2]
    } else {
        (lengths[n / 2 - 1] + lengths[n / 2]) / 2
    };

    writeln!(out, "records      {n}")?;
    writeln!(out, "bases        {total_bases}")?;
    writeln!(out, "min          {min}")?;
    writeln!(out, "median       {median}")?;
    writeln!(out, "mean         {mean:.1}")?;
    writeln!(out, "max          {max}")?;
    writeln!(out, "N50          {}", n50(&lengths, total_bases))?;
    if acgt > 0 {
        writeln!(out, "GC           {:.2}%", 100.0 * gc as f64 / acgt as f64)?;
    }
    if ambiguous > 0 {
        writeln!(
            out,
            "ambiguous    {ambiguous} ({:.3}%)",
            100.0 * ambiguous as f64 / total_bases as f64
        )?;
    }
    if qual_n > 0 {
        writeln!(out, "mean Q       {:.1}", qual_sum / qual_n as f64)?;
    }
    Ok(())
}

/// N50: the length such that half of all bases sit in records at least that
/// long. Expects `lengths` sorted ascending.
fn n50(lengths: &[usize], total: usize) -> usize {
    let half = total / 2;
    let mut acc = 0usize;
    for &l in lengths.iter().rev() {
        acc += l;
        if acc >= half {
            return l;
        }
    }
    0
}
