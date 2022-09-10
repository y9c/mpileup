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

fn complement_base_code(c: u8) -> u8 {
    // A = 65, a = 97
    // C = 67, c = 99
    // G = 71, g = 103
    // T = 84, t = 116
    // U = 85, u = 117
    match c {
        65 | 97 | 85 | 117 => 84,
        67 | 99 => 71,
        71 | 103 => 67,
        84 | 116 => 65,
        _ => 78,
    }
}

fn parse_region(
    chrom: &String,
    start: u32,
    end: u32,
    dna_bases: &Vec<u8>,
    fa_reader: &faidx::Reader,
    bam_path_list: &Vec<PathBuf>,
    min_depth: u32,
    count_indel: bool,
    ignore_strand: bool,
    by_strand: bool,
) {
    let mut p2depth: HashMap<(u32, usize), (u32, u32)> = HashMap::new();
    let mut p2base: HashMap<(u32, usize), (Vec<usize>, Vec<usize>)> = HashMap::new();
    let mut p2ins: HashMap<(u32, usize), (Vec<u32>, Vec<u32>)> = HashMap::new();
    let mut p2del: HashMap<(u32, usize), (Vec<u32>, Vec<u32>)> = HashMap::new();

    for (i, bam_path) in bam_path_list.iter().enumerate() {
        // read bam file
        let mut bam_reader = bam::IndexedReader::from_path(&bam_path).unwrap();
        let bam_header = bam_reader.header().clone();

        let tid = bam_header.tid(chrom.as_bytes()).unwrap();

        bam_reader.fetch((tid, start, end)).unwrap();
        // pileup over all covered sites
        for p in bam_reader.pileup() {
            let pileup = p.unwrap();

            let mut base_list_fwd: Vec<u8> = Vec::new();
            let mut base_list_rev: Vec<u8> = Vec::new();
            let mut insertion_list_fwd: Vec<u32> = Vec::new();
            let mut insertion_list_rev: Vec<u32> = Vec::new();
            let mut deletion_list_fwd: Vec<u32> = Vec::new();
            let mut deletion_list_rev: Vec<u32> = Vec::new();
            let mut total_reads_fwd = 0;
            let mut total_reads_rev = 0;
            let ref_pos = pileup.pos();
            if (start <= ref_pos) && (ref_pos < end) {
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

                    let strand: char;
                    if alignment.record().flags() & 128 == 128 {
                        if alignment.record().flags() & 16 == 16 {
                            strand = '+'
                        } else {
                            strand = '-'
                        }
                    } else {
                        if alignment.record().flags() & 16 == 16 {
                            strand = '-'
                        } else {
                            strand = '+'
                        }
                    }

                    if !alignment.is_del() && !alignment.is_refskip() {
                        let read_base = alignment.record().seq()[alignment.qpos().unwrap()];
                        if strand == '+' {
                            total_reads_fwd += 1;
                            base_list_fwd.push(read_base);
                        } else if strand == '-' {
                            total_reads_rev += 1;
                            base_list_rev.push(read_base);
                        }
                    }
                    if count_indel {
                        match alignment.indel() {
                            bam::pileup::Indel::Ins(len) => {
                                if strand == '+' {
                                    insertion_list_fwd.push(len);
                                } else if strand == '-' {
                                    insertion_list_rev.push(len);
                                }
                            }
                            bam::pileup::Indel::Del(len) => {
                                if strand == '+' {
                                    deletion_list_fwd.push(len);
                                } else if strand == '-' {
                                    deletion_list_rev.push(len);
                                }
                            }
                            _ => {}
                        }
                    }
                }
                p2depth.insert((ref_pos, i), (total_reads_fwd, total_reads_rev));

                let base_counter_fwd: Vec<usize>;
                let base_counter_rev: Vec<usize>;
                // count forward bases
                base_counter_fwd = dna_bases
                    .clone()
                    .into_iter()
                    .map(|b| base_list_fwd.iter().filter(|&x| *x == b).count())
                    .collect::<Vec<_>>();
                // count reverse bases
                base_counter_rev = dna_bases
                    .clone()
                    .into_iter()
                    .map(|b| base_list_rev.iter().filter(|&x| *x == b).count())
                    .collect::<Vec<_>>();

                p2base.insert((ref_pos, i), (base_counter_fwd, base_counter_rev));

                if count_indel {
                    p2ins.insert((ref_pos, i), (insertion_list_fwd, insertion_list_rev));
                    p2del.insert((ref_pos, i), (deletion_list_fwd, deletion_list_rev));
                }
            }
        }
    }

    for p in start..end {
        let rec_list = (0..bam_path_list.len())
            .map(|x| {
                if ignore_strand {
                    let mut rec = vec![match p2base.get(&(p, x)) {
                        Some((v1, v2)) => format!(
                            "{}",
                            (0..4)
                                .map(|i| v1[i] + v2[i])
                                .collect::<Vec<usize>>()
                                .iter()
                                .join(",")
                        ),
                        None => "0,0,0,0".to_string(),
                    }];
                    if count_indel {
                        rec.push(match p2ins.get(&(p, x)) {
                            Some((v1, v2)) => {
                                format!("{}", v1.iter().chain(v2.iter()).join("|"))
                            }
                            None => "".to_string(),
                        });
                        rec.push(match p2del.get(&(p, x)) {
                            Some((v1, v2)) => {
                                format!("{}", v1.iter().chain(v2.iter()).join("|"))
                            }
                            None => "".to_string(),
                        });
                    }
                    vec![rec.join(",")]
                } else if by_strand {
                    let mut rec = match p2base.get(&(p, x)) {
                        Some((v1, v2)) => vec![v1.iter().join(","), v2.iter().join(",")],
                        None => vec!["0,0,0,0".to_string(), "0,0,0,0".to_string()],
                    };
                    if count_indel {
                        rec.append(&mut match p2ins.get(&(p, x)) {
                            Some((v1, v2)) => vec![v1.iter().join("|"), v2.iter().join("|")],
                            None => vec!["".to_string(), "".to_string()],
                        });
                        rec.append(&mut match p2del.get(&(p, x)) {
                            Some((v1, v2)) => vec![v1.iter().join("|"), v2.iter().join("|")],
                            None => vec!["".to_string(), "".to_string()],
                        });
                    }
                    vec![
                        rec.iter().step_by(2).join(","),
                        rec.iter().skip(1).step_by(2).join(","),
                    ]
                } else {
                    let mut rec = vec![match p2base.get(&(p, x)) {
                        Some((v1, v2)) => {
                            format!("{},{}", v1.iter().join(","), v2.iter().join(","))
                        }
                        None => "0,0,0,0,0,0,0,0".to_string(),
                    }];
                    if count_indel {
                        rec.push(match p2ins.get(&(p, x)) {
                            Some((v1, v2)) => {
                                format!("{}|{}", v1.iter().join("|"), v2.iter().join("|"))
                            }
                            None => "".to_string(),
                        });
                        rec.push(match p2del.get(&(p, x)) {
                            Some((v1, v2)) => {
                                format!("{}|{}", v1.iter().join("|"), v2.iter().join("|"))
                            }
                            None => "".to_string(),
                        });
                    }
                    vec![rec.join(",")]
                }
            })
            .collect::<Vec<_>>();

        let r = fa_reader
            .fetch_seq_string(&chrom, p as usize, p as usize)
            .unwrap();
        if ignore_strand {
            // filter depth
            let passed_filter = (0..bam_path_list.len())
                .map(|x| match p2depth.get(&(p, x)) {
                    Some(val) => (*val).0 + (*val).1,
                    None => 0,
                })
                .max()
                .unwrap()
                >= min_depth;
            if passed_filter {
                let val = rec_list.iter().map(|x| &x[0]).join("\t");
                println!("{}\t{}\t{}\t{}\t{}", chrom, p, '.', r, val);
            }
        } else if by_strand {
            let passed_filter = (0..bam_path_list.len())
                .map(|x| match p2depth.get(&(p, x)) {
                    Some(val) => (*val).0,
                    None => 0,
                })
                .max()
                .unwrap()
                >= min_depth;
            if passed_filter {
                let val = rec_list.iter().map(|x| &x[0]).join("\t");
                println!("{}\t{}\t{}\t{}\t{}", chrom, p, "+", r, val);
            }
            let passed_filter = (0..bam_path_list.len())
                .map(|x| match p2depth.get(&(p, x)) {
                    Some(val) => (*val).1,
                    None => 0,
                })
                .max()
                .unwrap()
                >= min_depth;
            if passed_filter {
                let val = rec_list.iter().map(|x| &x[1]).join("\t");
                println!(
                    "{}\t{}\t{}\t{}\t{}",
                    chrom,
                    p,
                    "-",
                    complement_base_code(r.as_bytes()[0]) as char,
                    val
                );
            }
        } else {
            let passed_filter = (0..bam_path_list.len())
                .map(|x| match p2depth.get(&(p, x)) {
                    Some(val) => (*val).0 + (*val).1,
                    None => 0,
                })
                .max()
                .unwrap()
                >= min_depth;
            if passed_filter {
                let val = rec_list.iter().map(|x| &x[0]).join("\t");
                println!("{}\t{}\t{}\t{}\t{}", chrom, p, "+/-", r, val);
            }
        }
    }
}
pub fn run(
    region_path: PathBuf,
    fasta_path: PathBuf,
    bam_path_list: Vec<PathBuf>,
    min_depth: u32,
    count_indel: bool,
    with_header: bool,
    ignore_strand: bool,
    by_strand: bool,
) {
    // check parameters
    if by_strand & ignore_strand {
        eprintln!("Output records by strand, but `--ignore-strand` flag is set.");
        std::process::exit(1);
    }

    let chunk_size = 51;
    // A, C, G, T
    let dna_bases = &vec![65, 67, 71, 84];

    // read fasta file
    let fa_reader = &faidx::Reader::from_path(&fasta_path).unwrap();

    // read region file
    let mut pos_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(&region_path)
        .unwrap();

    if with_header {
        print!("Chrom\tPos\tStrand\tRef");
        for pth in &bam_path_list {
            print!("\t{}", pth.to_str().unwrap().to_string())
        }
        print!("\n");
    }
    // Read through all records in region.
    for (_i, record) in pos_reader.deserialize().enumerate() {
        let record: PosRecord = record.unwrap();
        let chrom: &String = &record.chrom;
        let start = record.start as u32;
        let end = record.end as u32;

        // split into chunks
        let mut spans: Vec<(u32, u32)> = Vec::new();
        let mut splited_start = start;
        for _ in 1..(end - start) / chunk_size {
            splited_start = splited_start + chunk_size;
            spans.push((splited_start, splited_start + chunk_size));
        }
        spans.push((splited_start, end));
        for (s, e) in spans {
            parse_region(
                &chrom,
                s,
                e,
                dna_bases,
                fa_reader,
                &bam_path_list,
                min_depth,
                count_indel,
                ignore_strand,
                by_strand,
            );
        }
    }
}
