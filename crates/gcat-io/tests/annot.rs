use gcat_io::annot::{genes, header_locus, read_bed, read_bedgraph, read_gff3, Strand};
use std::io::Cursor;

/// Two transcripts of one gene sharing their last exon, with CDS on each.
const GFF: &str = "##gff-version 3
chrX\tRefSeq\tregion\t1\t1000\t.\t+\t.\tID=chrX:1..1000;Name=X
chrX\tGnomon\tgene\t100\t900\t.\t+\t.\tID=gene-A;Name=GENEA;gene=GENEA
chrX\tGnomon\tmRNA\t100\t900\t.\t+\t.\tID=rna-1;Parent=gene-A
chrX\tGnomon\texon\t100\t150\t.\t+\t.\tID=e1;Parent=rna-1
chrX\tGnomon\texon\t800\t900\t.\t+\t.\tID=e2;Parent=rna-1
chrX\tGnomon\tCDS\t820\t880\t.\t+\t0\tID=c1;Parent=rna-1
chrX\tGnomon\tmRNA\t300\t900\t.\t+\t.\tID=rna-2;Parent=gene-A
chrX\tGnomon\texon\t300\t320\t.\t+\t.\tID=e3;Parent=rna-2
chrX\tGnomon\texon\t790\t900\t.\t+\t.\tID=e4;Parent=rna-2
chrX\tGnomon\tCDS\t810\t880\t.\t+\t0\tID=c2;Parent=rna-2
chrX\tGnomon\tgene\t950\t990\t.\t-\t.\tID=gene-B;gene=GENEB;product=a%2C b
##FASTA
>chrX
ACGT
";

#[test]
fn gff3_is_converted_to_half_open() {
    let f = read_gff3(Cursor::new(GFF)).unwrap();
    let gene = f.iter().find(|f| f.id.as_deref() == Some("gene-A")).unwrap();
    assert_eq!((gene.start, gene.end), (99, 900));
    assert_eq!(gene.strand, Strand::Forward);
    assert_eq!(gene.name.as_deref(), Some("GENEA"));
}

#[test]
fn gff3_stops_at_fasta_and_falls_back_to_gene_name() {
    let f = read_gff3(Cursor::new(GFF)).unwrap();
    assert_eq!(f.len(), 11);
    assert_eq!(f.last().unwrap().name.as_deref(), Some("GENEB"));
    assert_eq!(f.last().unwrap().strand, Strand::Reverse);
}

#[test]
fn gff3_rejects_short_lines() {
    assert!(read_gff3(Cursor::new("chrX\tgene\t1\t2\n")).is_err());
}

#[test]
fn genes_collapse_transcripts() {
    let g = genes(&read_gff3(Cursor::new(GFF)).unwrap());
    assert_eq!(g.len(), 2);
    assert_eq!(g[0].name, "GENEA");
    // 800–900 and 790–900 overlap and merge; the two first exons stay apart.
    assert_eq!(g[0].exons, vec![(99, 150), (299, 320), (789, 900)]);
    assert_eq!(g[0].cds, Some((809, 880)));
    assert!(g[1].exons.is_empty() && g[1].cds.is_none());
}

#[test]
fn bed_is_already_half_open() {
    let f = read_bed(Cursor::new("track name=x\n# c\nchrX 2000 7074 orange_del 0 +\nchrX\t5\t6\n")).unwrap();
    assert_eq!(f.len(), 2);
    assert_eq!((f[0].start, f[0].end), (2000, 7074));
    assert_eq!(f[0].name.as_deref(), Some("orange_del"));
    assert_eq!(f[0].strand, Strand::Forward);
    assert_eq!(f[1].name, None);
}

#[test]
fn header_locus_reads_ncbi_ranges() {
    assert_eq!(
        header_locus("NC_018741.3:110430079-110439152 Felis catus chromosome X"),
        ("NC_018741.3", 110430078)
    );
    assert_eq!(header_locus("chr1:1,001-2,000"), ("chr1", 1000));
    assert_eq!(header_locus("NC_001422.1 phiX174"), ("NC_001422.1", 0));
    // A colon that is not a range leaves the name alone.
    assert_eq!(header_locus("weird:name"), ("weird:name", 0));
}

#[test]
fn bedgraph_bins_by_max_and_keeps_zero_runs() {
    let sig = read_bedgraph(Cursor::new(
        "track type=bedGraph name=\"orange (Fca126)\"\nchrX 0 40 30\nchrX 40 60 0\nchrX 60 100 20\nchr1 0 100 99\n",
    ))
    .unwrap();
    assert_eq!(sig.name.as_deref(), Some("orange (Fca126)"));
    assert_eq!(sig.binned("chrX", 0.0, 100.0, 5), vec![30.0, 30.0, 0.0, 20.0, 20.0]);
    // A window past the data has no values, which is different from zero.
    assert!(sig.binned("chrX", 200.0, 300.0, 2).iter().all(|v| v.is_nan()));
    assert!(read_bedgraph(Cursor::new("chrX 0 10\n")).is_err());
}

#[test]
fn percent_signal_collapses_runs_and_skips_gaps() {
    use gcat_io::annot::Signal;
    let s = Signal::from_percent("tiger", "chrX", 100, &[99.6, 100.0, f32::NAN, 50.2, 50.0]);
    assert!(s.percent);
    assert_eq!(s.spans, vec![("chrX".into(), 100, 102, 100.0), ("chrX".into(), 103, 105, 50.0)]);
}
