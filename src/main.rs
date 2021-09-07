mod count;
use bio_types::strand::ReqStrand;
use rust_htslib::{bam, bam::Read};

use clap::{AppSettings, Clap};

/// This doc string acts as a help message when the user runs '--help'
/// as do all doc strings on fields
#[derive(Clap)]
#[clap(version = "0.0.1", author = "Chang Ye <yech1990@gmail.com>")]
#[clap(setting = AppSettings::ColoredHelp)]
struct Opts {
    /// Sets a custom config file. Could have been an Option<T> with no default too
    #[clap(short, long, default_value = "default.conf")]
    config: String,
    /// Some input. Because this isn't an Option<T> it's required to be used
    /// A level of verbosity, and can be used multiple times
    #[clap(short, long, parse(from_occurrences))]
    verbose: i32,
    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Clap)]
enum SubCommand {
    #[clap(version = "0.0.1", author = "Chang Ye <yech1990@gmail.com>")]
    Mpileup,
    #[clap(version = "0.0.1", author = "Chang Ye <yech1990@gmail.com>")]
    Count(Count),
}

#[derive(Clap)]
struct Count {
    #[clap(short, long, about = "debug")]
    debug: bool,
    #[clap(short = 'r', long = "region", about = "input bed file..")]
    bed: String,
    #[clap(short = 'i', long = "input", about = "input bam file..")]
    bam: String,
}

/// A subcommand for controlling testing
impl SubCommand {
    /// Print debug info
    fn run(&self) {
        match self {
            SubCommand::Mpileup => {
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
            SubCommand::Count(..) => {}
        }
    }
}

fn main() {
    let opts: Opts = Opts::parse();

    // Gets a value for config if supplied by user, or defaults to "default.conf"
    println!("Value for config: {}", opts.config);

    // Vary the output based on how many times the user used the "verbose" flag
    // (i.e. 'myprog -v -v -v' or 'myprog -vvv' vs 'myprog -v'
    match opts.verbose {
        0 => println!("No verbose info"),
        1 => println!("Some verbose info"),
        2 => println!("Tons of verbose info"),
        _ => println!("Don't be ridiculous"),
    }

    // You can handle information about subcommands by requesting their matches by name
    // (as below), requesting just the name used, or both at the same time
    match opts.subcmd {
        SubCommand::Mpileup => {
            println!("Printing debug info of mpileup...");
            SubCommand::Mpileup.run();
        }
        SubCommand::Count(c) => {
            println!("Printing debug info of count...");
            count::do_count(c.bed, c.bam, "tes2".to_string());
        }
    }
}
