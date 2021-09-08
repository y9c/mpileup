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

pub fn run(region_path: PathBuf, bam_path_list: Vec<PathBuf>, min_mean_depth: u32) {
    // A, C, G, T
    let dna_bases: Vec<u8> = vec![65, 67, 71, 84];

    // read region file
    let mut pos_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(&region_path)
        .unwrap();

    // Read through all records in region.
    for (_i, record) in pos_reader.deserialize().enumerate() {
        let record: PosRecord = record.unwrap();
        let start = record.start as u32;
        let end = record.end as u32;

        let mut p2base: HashMap<(u32, usize), String> = HashMap::new();
        let mut p2ins: HashMap<(u32, usize), String> = HashMap::new();
        let mut p2del: HashMap<(u32, usize), String> = HashMap::new();
        let mut p2depth: HashMap<(u32, usize), u32> = HashMap::new();

        for (i, bam_path) in bam_path_list.iter().enumerate() {
            // read bam file
            let mut bam_reader = bam::IndexedReader::from_path(&bam_path).unwrap();
            let bam_header = bam_reader.header().clone();

            let tid = bam_header.tid(record.chrom.as_bytes()).unwrap();

            bam_reader.fetch((tid, start, end)).unwrap();
            // pileup over all covered sites
            for p in bam_reader.pileup() {
                let pileup = p.unwrap();

                p2depth.insert((pileup.pos(), i), pileup.depth());

                let mut base_list: Vec<u8> = Vec::new();
                let mut insertion_list: Vec<u32> = Vec::new();
                let mut deletion_list: Vec<u32> = Vec::new();
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
                        match alignment.indel() {
                            bam::pileup::Indel::Ins(len) => {
                                insertion_list.push(len);
                            }
                            bam::pileup::Indel::Del(len) => {
                                deletion_list.push(len);
                            }
                            _ => {}
                        }
                    }

                    // let base_string = String::from_utf8(base_list.clone()).unwrap();
                    let base_counter = dna_bases
                        .clone()
                        .into_iter()
                        .map(|b| base_list.iter().filter(|&x| *x == b).count().to_string())
                        .collect::<Vec<_>>()
                        .join(",");
                    p2base.insert((pileup.pos(), i), base_counter);

                    let insertion_counter = insertion_list
                        .into_iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join("|");
                    p2ins.insert((pileup.pos(), i), insertion_counter);

                    let deletion_counter = deletion_list
                        .into_iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join("|");
                    p2del.insert((pileup.pos(), i), deletion_counter);
                }
            }
        }

        for p in start..end {
            // filter depth
            let dep: u32 = (0..bam_path_list.len())
                .map(|x| match p2depth.get(&(p, x)) {
                    Some(val) => *val,
                    None => 0,
                })
                .sum();
            if (dep / bam_path_list.len() as u32) < min_mean_depth {
                break;
            }

            print!("{}\t{}\t", record.chrom, p);
            let val = (0..bam_path_list.len())
                .map(|x| {
                    format!(
                        "{},{},{}",
                        match p2base.get(&(p, x)) {
                            Some(val) => val,
                            None => "0,0,0,0",
                        },
                        match p2ins.get(&(p, x)) {
                            Some(val) => val,
                            None => "",
                        },
                        match p2del.get(&(p, x)) {
                            Some(val) => val,
                            None => "",
                        }
                    )
                })
                .collect::<Vec<_>>()
                .join("\t");

            println!("{}", val);
        }
    }
}
