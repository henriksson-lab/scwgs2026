use std::{num::NonZeroUsize, path::PathBuf};

use clap::Parser;
use sag2tirp::{Options, convert};

#[derive(Parser)]
#[command(
    version,
    about = "Convert SAG assemblies to sorted, indexed Bascet TIRP (synthetic qualities)"
)]
struct Args {
    /// GORG-Dark ZIP, directory of per-cell FASTAs, or one FASTA (.gz supported)
    #[arg(short, long)]
    input: PathBuf,
    /// Output .tirp.gz; also creates .tbi and .records.tsv; refuses existing files
    #[arg(short, long)]
    output: PathBuf,
    /// Maximum fragment length; overlap is k-1 and the final tail is retained
    #[arg(long, default_value = "150")]
    fragment_length: NonZeroUsize,
    /// K-mer length; defaults to 31, giving 30-base overlaps
    #[arg(long, default_value = "31")]
    kmer_length: NonZeroUsize,
    /// Synthetic Phred+33 score (0..=93); default Q0, NOT measured quality
    #[arg(long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=93))]
    phred: u8,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    let stats = convert(
        &args.input,
        &args.output,
        &Options {
            fragment_length: args.fragment_length,
            kmer_length: args.kmer_length,
            phred: args.phred,
        },
    )?;
    eprintln!(
        "Wrote {} cells, {} contigs, {} records, {} source bases, {} emitted bases (synthetic Q{}).",
        stats.cells, stats.contigs, stats.records, stats.source_bases, stats.bases, args.phred
    );
    Ok(())
}
