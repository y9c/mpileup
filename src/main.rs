mod base;
mod count;

use clap::Parser;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[clap(about, version, author)]
struct Opts {
    /// A level of verbosity, and can be used multiple times
    #[clap(short, long, parse(from_occurrences))]
    verbose: i32,
    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Parser)]
enum SubCommand {
    Base(Base),
    Count(Count),
}

#[derive(Parser)]
struct Base {
    #[clap(short = 'g', long = "count-indel", help = "Count indel?")]
    indel: bool,
    #[clap(short = 't', long = "target", help = "input bed file..", validator = file_path_validation)]
    bed: PathBuf,
    #[clap(short = 'r', long = "reference", help = "input fa file..", validator = file_path_validation)]
    fa: PathBuf,
    #[clap(
        short = 'i',
        long = "input",
        help = "input bam files..",
        required = true,
        parse(from_os_str),
        takes_value = true, multiple_values = true, validator = file_path_validation,
    )]
    bam: Vec<PathBuf>,
    #[clap(
        short = 'd',
        long = "min-depth",
        help = "Set min depth for output. Note: anyone of the samples passing the cutoff is ok.",
        default_value = "0"
    )]
    min_depth: u32,
    #[clap(
        short = 'm',
        long = "mean-depth",
        help = "Set min cutoff of mean depth for output",
        default_value = "0"
    )]
    mean_depth: u32,
    #[clap(
        short = 'q',
        long = "min-qual",
        help = "Set min quality for base. (greater or equal to)",
        default_value = "0"
    )]
    qual: u8,
    #[clap(
        short = 'H',
        long = "headless",
        help = "Write without header in the output"
    )]
    headless: bool,
    #[clap(
        short = 'S',
        long = "strandless",
        help = "Calcualte counts ignore strand info"
    )]
    strandless: bool,
    #[clap(
        short = 's',
        long = "split-strand",
        help = "Split counts into different rows by strand"
    )]
    bystrand: bool,
    #[clap(
        short = 'c',
        long = "--chunk-size",
        default_value = "8",
        help = "Start the job in multiple threads"
    )]
    chunk: u32,
    #[clap(
        short = 'j',
        long = "--threads",
        default_value = "8",
        help = "Start the job in multiple threads"
    )]
    njobs: usize,
    #[clap(
        short = 'l',
        long = "--log-type",
        default_value = "0",
        help = "Log type. 0: no log; 1: spans; 2: progress bar"
    )]
    logtype: u8,
}

#[derive(Parser)]
struct Count {
    #[clap(short, long, help = "debug")]
    debug: bool,
    #[clap(short = 't', long = "target", help = "input bed file..")]
    bed: PathBuf,
    #[clap(short = 'r', long = "reference", help = "input fa file..")]
    fa: PathBuf,
    #[clap(short = 'i', long = "input", help = "input bam file..")]
    bam: PathBuf,
}

impl SubCommand {}

fn main() {
    let opts: Opts = Opts::parse();

    // Vary the output based on how many times the user used the "verbose" flag
    // (i.e. 'myprog -v -v -v' or 'myprog -vvv' vs 'myprog -v'
    match opts.verbose {
        0 => print!(""),
        1 => println!("Some verbose info"),
        2 => println!("Tons of verbose info"),
        _ => println!("Don't be ridiculous"),
    }

    // You can handle information about subcommands by requesting their matches by name
    // (as below), requesting just the name used, or both at the same time
    match opts.subcmd {
        SubCommand::Base(o) => {
            base::run(
                o.bed,
                o.fa,
                o.bam,
                o.min_depth,
                o.mean_depth,
                o.qual,
                o.indel,
                o.headless,
                o.strandless,
                o.bystrand,
                o.chunk,
                o.njobs,
                o.logtype,
            );
        }
        SubCommand::Count(o) => {
            count::run(o.bed, o.fa, o.bam);
        }
    }
}

fn file_path_validation(path: &str) -> Result<(), String> {
    let path = Path::new(path);
    if !path.exists() {
        Err(format!("{path:?} file doesn't exists"))
    } else if !path.is_file() {
        Err(format!("{path:?} is not a file"))
    } else {
        Ok(())
    }
}
