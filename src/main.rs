mod count;
mod mpu;

use clap::{AppSettings, Clap};
use std::path::PathBuf;

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
#[clap(setting = AppSettings::ColoredHelp)]
enum SubCommand {
    #[clap(version = "0.0.1", author = "Chang Ye <yech1990@gmail.com>")]
    Mpileup,
    #[clap(version = "0.0.1", author = "Chang Ye <yech1990@gmail.com>")]
    Count(Count),
}

#[derive(Clap)]
#[clap(setting = AppSettings::ColoredHelp)]
struct Count {
    #[clap(short, long, about = "debug")]
    debug: bool,
    #[clap(short = 'i', long = "input", about = "input bam file..")]
    bam: PathBuf,
    #[clap(short = 't', long = "target", about = "input bed file..")]
    bed: PathBuf,
    #[clap(short = 'r', long = "reference", about = "input fa file..")]
    fa: PathBuf,
}

impl SubCommand {}

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
            mpu::run();
        }
        SubCommand::Count(c) => {
            println!("Printing debug info of count...");
            count::run(c.bed, c.bam, c.fa);
        }
    }
}
