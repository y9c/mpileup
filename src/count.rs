use csv;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::faidx;
use serde::Deserialize;
use std::char;
use std::io;
use std::path::PathBuf;
use std::str;

#[derive(Deserialize, Debug)]
struct PosRecord {
    chrom: String,
    start: u32,
    end: u32,
}

pub fn run(region_path: PathBuf, cram_path: PathBuf, fasta_path: PathBuf) {
    let fa_reader = faidx::Reader::from_path(&fasta_path).unwrap();

    let mut bam_reader = bam::IndexedReader::from_path(&cram_path).unwrap();
    let bam_header = bam_reader.header().clone();

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

        for p in bam_reader.pileup() {
            let mut ref_count = 0;
            let mut alt_count = 0;
            let mut indel_count = 0;

            let pileup = p.unwrap();
            let ref_pos = pileup.pos() as usize;
            let r = fa_reader
                .fetch_seq_string(&record.chrom, ref_pos, ref_pos)
                .unwrap();

            for alignment in pileup.alignments() {
                if !alignment.is_del() && !alignment.is_refskip() {
                    /* println!("Base {}", alignment.record().seq()[alignment.qpos().unwrap()]); */
                    let base = alignment.record().seq()[alignment.qpos().unwrap()] as char;
                    if base.to_string().to_uppercase() != r {
                        alt_count += 1;
                    } else {
                        ref_count += 1;
                    }
                } else {
                    indel_count += 1;
                }
            }

            csv_writer
                .serialize((
                    &record.chrom,
                    ref_pos,
                    r,
                    pileup.depth(),
                    ref_count,
                    alt_count,
                    indel_count,
                ))
                .unwrap();
        }
    }
}
