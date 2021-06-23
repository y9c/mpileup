use bio_types::strand::ReqStrand;
use rust_htslib::{bam, bam::Read};

fn main() {
    let mut bam = bam::Reader::from_path(&"test/sample1.bam").unwrap();

    // pileup over all covered sites
    for p in bam.pileup() {
        let pileup = p.unwrap();
        println!(
            "Site = {}:{}; depth = {}",
            pileup.tid(),
            pileup.pos(),
            pileup.depth()
        );

        for alignment in pileup.alignments() {
            if !alignment.is_del()
                && !alignment.is_refskip()
                && alignment.record().strand() == ReqStrand::Forward
            {
                let read_base = alignment.record().seq()[alignment.qpos().unwrap()];
                println!("Base {}", read_base);
            }
        }
    }
}
