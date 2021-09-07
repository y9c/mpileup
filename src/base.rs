use bio_types::strand::ReqStrand;
use rust_htslib::{bam, bam::Read};
pub fn run() {
    let mut bam = bam::Reader::from_path(&"test/sample1.bam").unwrap();

    // pileup over all covered sites
    for p in bam.pileup() {
        let pileup = p.unwrap();

        let mut base_list: Vec<u8> = Vec::new();
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
        println!(
            "{}\t{}\t{}\t{}",
            pileup.tid(),
            pileup.pos(),
            pileup.depth(),
            base_string
        );
    }
}
