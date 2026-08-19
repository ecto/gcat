//! `gcat` — cat for genomes.
//!
//! Reads FASTA/FASTQ from files or stdin, writes to stdout, and composes with
//! pipes. Subcommands stay small and single-purpose on the theory that the
//! shell is a perfectly good way to combine them.

use gcat_io::{write_fasta, write_fastq};
use std::io::{self, BufRead, BufWriter, Write};
use std::process::ExitCode;

mod cmd;

const USAGE: &str = "\
gcat — cat for genomes

USAGE:
    gcat [COMMAND] [OPTIONS] [FILE...]

COMMANDS:
    cat      Concatenate records to stdout (default)
    stats    Summarise records: count, bases, N50, GC
    head     First N records
    grep     Records whose id or sequence matches a pattern
    rc       Reverse complement every record
    kmer     Count canonical k-mers

OPTIONS:
    -h, --help       Show this help
    -V, --version    Show version

With no FILE, or when FILE is -, reads standard input.
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return ExitCode::SUCCESS;
    }
    if args.iter().any(|a| a == "-V" || a == "--version") {
        println!("gcat {}", env!("CARGO_PKG_VERSION"));
        return ExitCode::SUCCESS;
    }

    // The first argument is a command only if it names one; otherwise it is a
    // file and the command defaults to `cat`, so `gcat reads.fq` does the
    // obvious thing.
    let (command, rest) = match args.split_first() {
        None => ("cat", &[][..]),
        Some((first, rest)) if cmd::is_command(first) => (first.as_str(), rest),
        _ => ("cat", &args[..]),
    };

    match run(command, rest) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            // A closed pipe is how `gcat ... | head` ends. That is success.
            if let Some(io_err) = e.downcast_ref::<io::Error>() {
                if io_err.kind() == io::ErrorKind::BrokenPipe {
                    return ExitCode::SUCCESS;
                }
            }
            eprintln!("gcat: {e}");
            ExitCode::FAILURE
        }
    }
}

pub type Fallible = Result<(), Box<dyn std::error::Error>>;

fn run(command: &str, args: &[String]) -> Fallible {
    let (flags, files) = cmd::split_args(args);
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    let result = match command {
        "cat" => cmd::cat(&files, &mut out),
        "stats" => cmd::stats(&files, &mut out),
        "head" => cmd::head(&files, &flags, &mut out),
        "grep" => cmd::grep(&files, &flags, &mut out),
        "rc" => cmd::rc(&files, &mut out),
        "kmer" => cmd::kmer(&files, &flags, &mut out),
        other => Err(format!("unknown command {other:?}\n\n{USAGE}").into()),
    };
    result?;
    out.flush()?;
    Ok(())
}

/// A named input stream.
pub type Input = (String, Box<dyn BufRead>);

/// Open each input as a buffered reader, treating no files or `-` as stdin.
///
/// A file that cannot be opened is an error, not an empty stream: silently
/// reading nothing would make `gcat stats typo.fa` report zero records and
/// exit 0, which is worse than useless.
pub fn inputs(files: &[String]) -> Result<Vec<Input>, Box<dyn std::error::Error>> {
    if files.is_empty() {
        return Ok(vec![(
            "-".into(),
            Box::new(io::BufReader::new(io::stdin())) as Box<dyn BufRead>,
        )]);
    }
    let mut out = Vec::with_capacity(files.len());
    for f in files {
        let r: Box<dyn BufRead> = if f == "-" {
            Box::new(io::BufReader::new(io::stdin()))
        } else {
            let file = std::fs::File::open(f).map_err(|e| format!("{f}: {e}"))?;
            Box::new(io::BufReader::new(file))
        };
        out.push((f.clone(), r));
    }
    Ok(out)
}

/// Echo a record in whatever format it arrived as.
pub fn echo<W: Write>(
    out: &mut W,
    id: &str,
    seq: &[u8],
    qual: Option<&[u8]>,
    width: usize,
) -> io::Result<()> {
    match qual {
        Some(q) => write_fastq(out, id, seq, q),
        None => write_fasta(out, id, seq, width),
    }
}
