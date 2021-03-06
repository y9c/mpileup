use csv;
use itertools::Itertools;
use rust_htslib::bam::{self, Read};
use rust_htslib::faidx;
use serde::Deserialize;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Deserialize, Debug)]
struct PosRecord {
    chrom: String,
    start: u32,
    end: u32,
}

// fn complement_base_code(c: u8) -> u8 {
//     // A = 65, a = 97
//     // C = 67, c = 99
//     // G = 71, g = 103
//     // T = 84, t = 116
//     // U = 85, u = 117
//     match c {
//         65 | 97 | 85 | 117 => 84,
//         67 | 99 => 71,
//         71 | 103 => 67,
//         84 | 116 => 65,
//         _ => 78,
//     }
// }

pub fn run(
    region_path: PathBuf,
    fasta_path: PathBuf,
    bam_path_list: Vec<PathBuf>,
    min_mean_depth: u32,
    count_indel: bool,
) {
    // A, C, G, T
    let dna_bases: Vec<u8> = vec![65, 67, 71, 84];

    // read fasta file
    let fa_reader = faidx::Reader::from_path(&fasta_path).unwrap();

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

                let mut base_list_fwd: Vec<u8> = Vec::new();
                let mut base_list_rev: Vec<u8> = Vec::new();
                let mut insertion_list: Vec<u32> = Vec::new();
                let mut deletion_list: Vec<u32> = Vec::new();
                let ref_pos = pileup.pos();
                if (start <= ref_pos) && (ref_pos < end) {
                    let mut total_reads = 0;

                    // for alignment in pileup.alignments() {
                    // TODO: pick by quality
                    // START: group by qname
                    let grouped_by_qname = pileup
                        .alignments()
                        .map(|aln| {
                            let record = aln.record();
                            (aln, record)
                        })
                        .sorted_by(|a, b| Ord::cmp(a.1.qname(), b.1.qname()))
                        .group_by(|a| a.1.qname().to_owned());

                    for (_qname, reads) in grouped_by_qname.into_iter() {
                        let (alignment, _record) = reads
                            .into_iter()
                            .max_by(|a, b| match a.1.mapq().cmp(&b.1.mapq()) {
                                Ordering::Greater => Ordering::Greater,
                                Ordering::Less => Ordering::Less,
                                Ordering::Equal => {
                                    if a.1.flags() & 64 == 0 {
                                        Ordering::Greater
                                    } else if b.1.flags() & 64 == 0 {
                                        Ordering::Less
                                    } else {
                                        Ordering::Greater
                                    }
                                }
                            })
                            .unwrap();
                        // END: group by qname

                        if !alignment.is_del() && !alignment.is_refskip() {
                            total_reads += 1;
                            let read_base = alignment.record().seq()[alignment.qpos().unwrap()];
                            if alignment.record().flags() & 128 == 128 {
                                if alignment.record().flags() & 16 == 16 {
                                    base_list_fwd.push(read_base);
                                } else {
                                    base_list_rev.push(read_base);
                                }
                            } else {
                                if alignment.record().flags() & 16 == 16 {
                                    base_list_rev.push(read_base);
                                } else {
                                    base_list_fwd.push(read_base);
                                }
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
                    p2depth.insert((ref_pos, i), total_reads);

                    // let base_string = String::from_utf8(base_list.clone()).unwrap();

                    let mut base_counter = dna_bases
                        .clone()
                        .into_iter()
                        .map(|b| {
                            base_list_fwd
                                .iter()
                                .filter(|&x| *x == b)
                                .count()
                                .to_string()
                        })
                        .collect::<Vec<_>>();
                    base_counter.append(
                        &mut dna_bases
                            .clone()
                            .into_iter()
                            .map(|b| {
                                base_list_rev
                                    .iter()
                                    .filter(|&x| *x == b)
                                    .count()
                                    .to_string()
                            })
                            .collect::<Vec<_>>(),
                    );
                    p2base.insert((ref_pos, i), base_counter.join(","));

                    if count_indel {
                        let insertion_counter = insertion_list
                            .into_iter()
                            .map(|x| x.to_string())
                            .collect::<Vec<_>>()
                            .join("|");
                        p2ins.insert((ref_pos, i), insertion_counter);

                        let deletion_counter = deletion_list
                            .into_iter()
                            .map(|x| x.to_string())
                            .collect::<Vec<_>>()
                            .join("|");
                        p2del.insert((ref_pos, i), deletion_counter);
                    }
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

            let r = fa_reader
                .fetch_seq_string(&record.chrom, p as usize, p as usize)
                .unwrap();
            if (dep / bam_path_list.len() as u32) >= min_mean_depth {
                print!("{}\t{}\t{}\t", record.chrom, p, r);
                let val = (0..bam_path_list.len())
                    .map(|x| {
                        if count_indel {
                            format!(
                                "{},{},{}",
                                match p2base.get(&(p, x)) {
                                    Some(val) => val,
                                    None => "0,0,0,0,0,0,0,0",
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
                        } else {
                            match p2base.get(&(p, x)) {
                                Some(val) => val,
                                None => "0,0,0,0,0,0,0,0",
                            }
                            .to_string()
                        }
                    })
                    .collect::<Vec<_>>()
                    .join("\t");

                println!("{}", val);
            }
        }
    }
}
