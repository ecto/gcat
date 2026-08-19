use gcat_io::{write_fasta, write_fastq, Format, Reader, PHRED33};
use std::io::Cursor;

/// id, sequence, optional quality — flattened so records outlive the reader.
type Owned = (String, String, Option<String>);

fn read_all(src: &str) -> (Option<Format>, Vec<Owned>) {
    let mut r = Reader::new(Cursor::new(src.as_bytes().to_vec()));
    let mut out = Vec::new();
    while let Some(rec) = r.next_record().unwrap() {
        out.push((
            rec.id.to_string(),
            String::from_utf8(rec.seq.to_vec()).unwrap(),
            rec.qual.map(|q| String::from_utf8(q.to_vec()).unwrap()),
        ));
    }
    (r.format(), out)
}

#[test]
fn reads_single_line_fasta() {
    let (fmt, recs) = read_all(">read1 desc here\nACGTACGT\n>read2\nTTTT\n");
    assert_eq!(fmt, Some(Format::Fasta));
    assert_eq!(recs.len(), 2);
    assert_eq!(recs[0].0, "read1 desc here");
    assert_eq!(recs[0].1, "ACGTACGT");
    assert_eq!(recs[1].1, "TTTT");
    assert!(recs[0].2.is_none());
}

#[test]
fn joins_multiline_fasta() {
    let (_, recs) = read_all(">r\nACGT\nACGT\nAC\n>s\nTTTT\n");
    assert_eq!(recs[0].1, "ACGTACGTAC");
    assert_eq!(recs[1].1, "TTTT");
}

#[test]
fn reads_fastq_with_quality() {
    let (fmt, recs) = read_all("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\n!!!!\n");
    assert_eq!(fmt, Some(Format::Fastq));
    assert_eq!(recs.len(), 2);
    assert_eq!(recs[0].1, "ACGT");
    assert_eq!(recs[0].2.as_deref(), Some("IIII"));
    assert_eq!(recs[1].2.as_deref(), Some("!!!!"));
}

#[test]
fn record_name_stops_at_whitespace() {
    let mut r = Reader::new(Cursor::new(b">read1 some description\nACGT\n".to_vec()));
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(rec.name(), "read1");
    assert_eq!(rec.id, "read1 some description");
}

#[test]
fn skips_blank_lines_between_records() {
    let (_, recs) = read_all(">a\nACGT\n\n\n>b\nTTTT\n\n");
    assert_eq!(recs.len(), 2);
}

#[test]
fn empty_input_yields_no_records() {
    let (_, recs) = read_all("");
    assert!(recs.is_empty());
}

#[test]
fn quality_length_mismatch_is_an_error() {
    let mut r = Reader::new(Cursor::new(b"@r\nACGT\n+\nII\n".to_vec()));
    assert!(r.next_record().is_err());
}

#[test]
fn truncated_fastq_is_an_error() {
    let mut r = Reader::new(Cursor::new(b"@r\nACGT\n+\n".to_vec()));
    assert!(r.next_record().is_err());
}

#[test]
fn bad_header_is_an_error() {
    let mut r = Reader::new(Cursor::new(b"ACGT\nACGT\n".to_vec()));
    assert!(r.next_record().is_err());
}

#[test]
fn handles_crlf_line_endings() {
    let (_, recs) = read_all("@r\r\nACGT\r\n+\r\nIIII\r\n");
    assert_eq!(recs[0].1, "ACGT");
    assert_eq!(recs[0].2.as_deref(), Some("IIII"));
}

#[test]
fn mean_quality_averages_in_probability_space() {
    // Half the bases at Q40 (p=1e-4), half at Q10 (p=0.1).
    // Mean p = 0.05005 -> Q ~= 13.0, NOT the naive (40+10)/2 = 25.
    let src = format!("@r\nACGT\n+\n{}{}{}{}\n", 'I', 'I', '+', '+');
    let mut r = Reader::new(Cursor::new(src.into_bytes()));
    let rec = r.next_record().unwrap().unwrap();
    let q = rec.mean_quality(PHRED33).unwrap();
    assert!(q > 12.0 && q < 14.0, "expected ~13, got {q}");
}

#[test]
fn fasta_records_have_no_mean_quality() {
    let mut r = Reader::new(Cursor::new(b">r\nACGT\n".to_vec()));
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(rec.mean_quality(PHRED33), None);
}

#[test]
fn records_convert_to_packed_seq() {
    let mut r = Reader::new(Cursor::new(b">r\nACGTNNGT\n".to_vec()));
    let rec = r.next_record().unwrap().unwrap();
    assert!(
        rec.to_seq().is_err(),
        "ambiguity should be rejected by default"
    );
    let (seq, n) = rec.to_seq_lossy(gcat::Base::A).unwrap();
    assert_eq!(n, 2);
    assert_eq!(seq.to_string(), "ACGTAAGT");
}

#[test]
fn fasta_writer_wraps_at_width() {
    let mut buf = Vec::new();
    write_fasta(&mut buf, "r", b"ACGTACGTAC", 4).unwrap();
    assert_eq!(String::from_utf8(buf).unwrap(), ">r\nACGT\nACGT\nAC\n");
}

#[test]
fn fasta_writer_width_zero_is_one_line() {
    let mut buf = Vec::new();
    write_fasta(&mut buf, "r", b"ACGTACGTAC", 0).unwrap();
    assert_eq!(String::from_utf8(buf).unwrap(), ">r\nACGTACGTAC\n");
}

#[test]
fn written_records_read_back_identically() {
    let mut buf = Vec::new();
    write_fastq(&mut buf, "r1 desc", b"ACGTACGT", b"IIIIIIII").unwrap();
    write_fastq(&mut buf, "r2", b"TTTT", b"!!!!").unwrap();

    let mut r = Reader::new(Cursor::new(buf));
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(rec.id, "r1 desc");
    assert_eq!(rec.seq, b"ACGTACGT");
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(rec.id, "r2");
    assert_eq!(rec.qual, Some(&b"!!!!"[..]));
}
