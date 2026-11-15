use anyhow::{Context, Result, ensure};
use clap::Parser;
use rand::{Rng, SeedableRng, rngs::SmallRng, seq::index};
use std::{
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    num::NonZeroUsize,
    path::{Path, PathBuf},
    time::Instant,
};

const READ_LEN: usize = 150;
const PREFIX_LEN: usize = 45;
const QUALITY: [u8; READ_LEN] = [b'F'; READ_LEN];
const BARCODE_TABLE: &str = include_str!("../atrandi_barcodes.tsv");
type Pools = [Vec<[u8; 8]>; 4];

#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    /// Number of distinct cell barcodes (1..=331776)
    #[arg(long)]
    cells: usize,

    /// Total read PAIRS; each output FASTQ contains this many 150 bp reads
    #[arg(long)]
    reads: u64,

    /// Output prefix; writes <prefix>_R1.fastq.gz and <prefix>_R2.fastq.gz (BGZF)
    #[arg(long, default_value = "simulated")]
    out_prefix: PathBuf,

    /// Random seed (repeatable for the same binary and arguments)
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Write plain .fastq files instead of BGZF (Bascet requires BGZF input)
    #[arg(long)]
    uncompressed: bool,

    /// BGZF compression workers per mate (two independent writer pools)
    #[arg(long, default_value = "8")]
    compression_threads: NonZeroUsize,
}

enum FastqWriter {
    Plain(BufWriter<File>),
    Bgzf(BufWriter<noodles_bgzf::io::MultithreadedWriter<BufWriter<File>>>),
}

impl FastqWriter {
    fn new(file: File, uncompressed: bool, workers: NonZeroUsize) -> Self {
        let buffered = BufWriter::with_capacity(1024 * 1024, file);
        if uncompressed {
            Self::Plain(buffered)
        } else {
            let compressed =
                noodles_bgzf::io::MultithreadedWriter::with_worker_count(workers, buffered);
            Self::Bgzf(BufWriter::with_capacity(1024 * 1024, compressed))
        }
    }

    fn finish(self) -> std::io::Result<()> {
        match self {
            Self::Plain(mut writer) => writer.flush(),
            Self::Bgzf(writer) => writer
                .into_inner()
                .map_err(|e| e.into_error())?
                .finish()?
                .flush(),
        }
    }
}

impl Write for FastqWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(bytes),
            Self::Bgzf(writer) => writer.write(bytes),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Bgzf(writer) => writer.flush(),
        }
    }
}

fn barcode_pools() -> Result<Pools> {
    let mut pools: Pools = std::array::from_fn(|_| Vec::new());
    for line in BARCODE_TABLE.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        ensure!(fields.len() == 3, "Invalid embedded barcode row");
        let round = fields[0].parse::<usize>()?;
        ensure!((1..=4).contains(&round), "Invalid barcode round");
        let barcode: [u8; 8] = fields[2].as_bytes().try_into()?;
        ensure!(
            barcode.iter().all(|b| b"ACGT".contains(b)),
            "Invalid barcode bases"
        );
        ensure!(!pools[round - 1].contains(&barcode), "Duplicate barcode");
        pools[round - 1].push(barcode);
    }
    ensure!(
        pools.iter().all(|p| p.len() == 24),
        "Expected 24 barcodes per round"
    );
    Ok(pools)
}

// Mixed-radix decoding gives every full barcode combination a unique index.
// Matches Bascet's atrandi_wgs_barcode_illumina.rs: round 4 at R2 offset 0,
// round 3 at 12, round 2 at 24, round 1 at 36; trim 45 bases from R2.
fn cell_prefix(mut combination: usize, pools: &Pools) -> [u8; PREFIX_LEN] {
    let mut prefix = [b'T'; PREFIX_LEN];
    for (position, pool) in pools.iter().rev().enumerate() {
        let choice = combination % pool.len();
        combination /= pool.len();
        prefix[position * 12..position * 12 + 8].copy_from_slice(&pool[choice]);
    }
    prefix[8..12].copy_from_slice(b"AGGA");
    prefix[20..24].copy_from_slice(b"ACTC");
    prefix[32..36].copy_from_slice(b"AAGG");
    prefix // terminal T is the ligation base
}

// One RNG call per 32 genomic bases, with no per-read heap allocation.
fn random_dna(rng: &mut SmallRng, sequence: &mut [u8]) {
    for chunk in sequence.chunks_mut(32) {
        let mut bits = rng.random::<u64>();
        for base in chunk {
            *base = b"ACGT"[(bits & 3) as usize];
            bits >>= 2;
        }
    }
}

fn output_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name = prefix.as_os_str().to_os_string();
    name.push(suffix);
    PathBuf::from(name)
}

fn create_output(path: &Path) -> Result<File> {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .with_context(|| {
            format!(
                "Cannot create {} (existing files are not overwritten)",
                path.display()
            )
        })
}

fn write_record(
    writer: &mut impl Write,
    id: u64,
    cell: usize,
    mate: u8,
    sequence: &[u8; READ_LEN],
) -> std::io::Result<()> {
    writeln!(writer, "@sim:{id}:cell{cell}/{mate}")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(&QUALITY)?;
    writer.write_all(b"\n")
}

fn run(args: Args) -> Result<()> {
    let pools = barcode_pools()?;
    let capacity: usize = pools.iter().map(Vec::len).product();
    ensure!(
        (1..=capacity).contains(&args.cells),
        "--cells must be between 1 and {capacity}"
    );
    ensure!(args.reads > 0, "--reads must be positive");
    ensure!(
        !args.out_prefix.as_os_str().is_empty(),
        "Output prefix must not be empty"
    );

    let mut rng = SmallRng::seed_from_u64(args.seed);
    let cells: Vec<_> = index::sample(&mut rng, capacity, args.cells)
        .into_iter()
        .map(|combination| cell_prefix(combination, &pools))
        .collect();

    let extension = if args.uncompressed {
        "fastq"
    } else {
        "fastq.gz"
    };
    let path_r1 = output_path(&args.out_prefix, &format!("_R1.{extension}"));
    let path_r2 = output_path(&args.out_prefix, &format!("_R2.{extension}"));
    if let Some(parent) = path_r1.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)?;
    }
    let file_r1 = create_output(&path_r1)?;
    let file_r2 = match create_output(&path_r2) {
        Ok(file) => file,
        Err(error) => {
            drop(file_r1);
            std::fs::remove_file(&path_r1)?;
            return Err(error);
        }
    };
    let mut r1_out = FastqWriter::new(file_r1, args.uncompressed, args.compression_threads);
    let mut r2_out = FastqWriter::new(file_r2, args.uncompressed, args.compression_threads);
    let start = Instant::now();
    let mut r1 = [b'A'; READ_LEN];
    let mut r2 = [b'A'; READ_LEN];
    eprintln!(
        "Generating {} read pairs from {} cells (seed {}); mean {:.2} pairs/cell",
        args.reads,
        args.cells,
        args.seed,
        args.reads as f64 / args.cells as f64
    );

    for id in 0..args.reads {
        // Independent uniform draws give random input order and approximately
        // equal coverage without constructing and shuffling all read pairs.
        let cell = rng.random_range(0..cells.len());
        random_dna(&mut rng, &mut r1);
        random_dna(&mut rng, &mut r2[PREFIX_LEN..]);
        r2[..PREFIX_LEN].copy_from_slice(&cells[cell]);
        write_record(&mut r1_out, id, cell, 1, &r1)?;
        write_record(&mut r2_out, id, cell, 2, &r2)?;
    }
    r1_out.finish()?;
    r2_out.finish()?;
    eprintln!(
        "Wrote {} and {} in {:.2}s",
        path_r1.display(),
        path_r2.display(),
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

fn main() -> Result<()> {
    run(Args::parse())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn all_cell_combinations_are_unique_and_match_chemistry() {
        let pools = barcode_pools().unwrap();
        let capacity: usize = pools.iter().map(Vec::len).product();
        assert_eq!(capacity, 331776);
        let mut seen = HashSet::new();
        for combination in 0..capacity {
            let prefix = cell_prefix(combination, &pools);
            for (position, pool) in pools.iter().rev().enumerate() {
                let barcode: [u8; 8] = prefix[position * 12..position * 12 + 8].try_into().unwrap();
                assert!(pool.contains(&barcode));
            }
            assert_eq!(&prefix[8..12], b"AGGA");
            assert_eq!(&prefix[20..24], b"ACTC");
            assert_eq!(&prefix[32..36], b"AAGG");
            assert_eq!(prefix[44], b'T');
            assert!(seen.insert(prefix));
        }
    }

    #[test]
    fn record_has_paired_name_and_150_bases_and_qualities() {
        let mut record = Vec::new();
        write_record(&mut record, 123, 42, 2, &[b'A'; READ_LEN]).unwrap();
        let record = String::from_utf8(record).unwrap();
        let lines: Vec<_> = record.lines().collect();
        assert_eq!(lines[0], "@sim:123:cell42/2");
        assert_eq!(lines[1], "A".repeat(150));
        assert_eq!(lines[2], "+");
        assert_eq!(lines[3], "F".repeat(150));
    }
}
