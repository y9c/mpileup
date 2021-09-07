use csv;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use serde::Deserialize;
use std::char;
use std::io;
use std::str;

#[derive(Deserialize, Debug)]
struct PosRecord {
    chrom: String,
    start: u32,
    end: u32,
}

pub fn do_count(
    region_path: std::string::String,
    cram_path: std::string::String,
    fasta_path: std::string::String,
) {
    let mut bam_reader = bam::IndexedReader::from_path(&cram_path).unwrap();
    let bam_header = bam_reader.header().clone();

    bam_reader.set_reference(fasta_path).unwrap();

    println!("{:?}", region_path);
    let mut pos_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(&region_path)
        .unwrap();

    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::BufWriter::new(io::stdout()));

    for (_i, record) in pos_reader.deserialize().enumerate() {
        let record: PosRecord = record.unwrap();

        // pileup over all covered ites
        let tid = bam_header.tid(record.chrom.as_bytes()).unwrap();
        let start = record.start as u32;
        let end = record.end as u32;
        bam_reader.fetch((tid, start, end)).unwrap();

        let mut alt_count = 0;

        for p in bam_reader.pileup() {
            let pileup = p.unwrap();

            for alignment in pileup.alignments() {
                if !alignment.is_del() && !alignment.is_refskip() {
                    /* println!("Base {}", alignment.record().seq()[alignment.qpos().unwrap()]); */
                    let base = alignment.record().seq()[alignment.qpos().unwrap()] as char;
                    if base.to_string().to_uppercase() != "A" {
                        alt_count += 1;
                    }
                }
            }

            csv_writer
                .serialize((&record.chrom, pileup.pos(), pileup.depth(), alt_count))
                .unwrap();
        }
    }
}
