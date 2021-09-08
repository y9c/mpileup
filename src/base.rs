use bio_types::strand::ReqStrand;
use csv;
use rust_htslib::bam::{self, Read};
use serde::Deserialize;
use std::path::PathBuf;

use std::collections::HashMap;

#[derive(Deserialize, Debug)]
struct PosRecord {
    chrom: String,
    start: u32,
    end: u32,
}

pub fn run(region_path: PathBuf, bam_path_list: Vec<PathBuf>) {
    // read region file
    let mut pos_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(&region_path)
        .unwrap();

    // Read through all records in region.
    for (_i, record) in pos_reader.deserialize().enumerate() {
        let record: PosRecord = record.unwrap();
        let start = record.start as u32;
        let end = record.end as u32;

        // let mut p2s: <u32, String> = HashMap::new();
        let mut p2s: HashMap<(u32, usize), String> = HashMap::new();

        for (i, bam_path) in bam_path_list.iter().enumerate() {
            // read bam file
            let mut bam_reader = bam::IndexedReader::from_path(&bam_path).unwrap();
            let bam_header = bam_reader.header().clone();

            let tid = bam_header.tid(record.chrom.as_bytes()).unwrap();

            bam_reader.fetch((tid, start, end)).unwrap();
            // pileup over all covered sites
            for p in bam_reader.pileup() {
                let pileup = p.unwrap();

                let mut base_list: Vec<u8> = Vec::new();
                if (start <= pileup.pos()) && (pileup.pos() < end) {
                    for alignment in pileup.alignments() {
                        if !alignment.is_del() && !alignment.is_refskip() {
                            if alignment.record().strand() == ReqStrand::Forward {
                                let read_base = alignment.record().seq()[alignment.qpos().unwrap()];
                                base_list.push(read_base);
                            } else if alignment.record().strand() == ReqStrand::Reverse {
                                let read_base = alignment.record().seq()[alignment.qpos().unwrap()];
                                base_list.push(read_base);
                            }
                        }
                    }

                    let base_string = String::from_utf8(base_list.clone()).unwrap();
                    p2s.insert((pileup.pos(), i), base_string);

                    // println!( "{}\t{}\t{}\t{}", pileup.tid(), pileup.pos(), pileup.depth(), base_string);
                }
            }
        }

        for p in start..end {
            print!("{}\t", p);
            let val = (0..bam_path_list.len())
                .map(|x| match p2s.get(&(p, x)) {
                    Some(val) => val,
                    None => "",
                })
                .collect::<Vec<_>>()
                .join("\t");

            println!("{}", val);
        }
    }
}
