//! Per-file SAG identity becomes TIRP cell identity. Assemblies are represented
//! as single-end records with explicitly synthetic quality scores.
mod bgzf;

use std::{
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    num::NonZeroUsize,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail, ensure};
use flate2::read::MultiGzDecoder;
use noodles::{
    core::Position,
    csi::binning_index::index::{header::Builder, reference_sequence::bin::Chunk},
    tabix,
};
use tempfile::NamedTempFile;
use zip::ZipArchive;

#[derive(Clone, Debug)]
pub struct Options {
    /// Fragment length, k..=30,000; defaults to 150 bases.
    pub fragment_length: NonZeroUsize,
    /// Fragments overlap by k-1 bases; defaults to 31.
    pub kmer_length: NonZeroUsize,
    /// Synthetic Phred score, 0..=93. No quality is inferred from FASTA.
    pub phred: u8,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            fragment_length: NonZeroUsize::new(150).unwrap(),
            kmer_length: NonZeroUsize::new(31).unwrap(),
            phred: 0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Stats {
    pub cells: u64,
    pub contigs: u64,
    pub records: u64,
    pub bases: u64,
    pub source_bases: u64,
}

enum Source {
    File(PathBuf),
    Zip(usize),
}
struct Input {
    cell: String,
    gzip: bool,
    source: Source,
}

/// Appends a sidecar suffix without replacing any part of the original path.
pub fn sidecar(output: &Path, suffix: &str) -> PathBuf {
    let mut name = output.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn input_name(name: &str, source: Source) -> Result<Option<Input>> {
    // ZIPs created on macOS contain AppleDouble records, not real gzip files.
    if name
        .split('/')
        .any(|p| p == "__MACOSX" || p.starts_with('.'))
    {
        return Ok(None);
    }
    let base = name.rsplit('/').next().unwrap_or(name);
    let (stem, gzip) = match base.strip_suffix(".gz") {
        Some(s) => (s, true),
        None => (base, false),
    };
    let Some(stem) = [".fasta", ".fna", ".fa"]
        .iter()
        .find_map(|s| stem.strip_suffix(s))
    else {
        return Ok(None);
    };
    let cell = stem.strip_suffix("_contigs").unwrap_or(stem);
    ensure!(
        !cell.is_empty()
            && cell.len() <= 1024
            && !cell.starts_with('#')
            && cell
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b)),
        "unsafe cell ID derived from filename {base:?}; use ASCII letters, digits, _, ., -"
    );
    Ok(Some(Input {
        cell: cell.to_owned(),
        gzip,
        source,
    }))
}

fn collect_directory(dir: &Path, inputs: &mut Vec<Input>) -> Result<()> {
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        if entry.file_name().as_encoded_bytes().starts_with(b".") || entry.file_name() == "__MACOSX"
        {
            continue;
        }
        let kind = entry.file_type()?;
        if kind.is_dir() {
            collect_directory(&entry.path(), inputs)?;
        } else if kind.is_file() {
            let name = entry
                .file_name()
                .into_string()
                .map_err(|_| anyhow::anyhow!("non-UTF8 filename"))?;
            if let Some(input) = input_name(&name, Source::File(entry.path()))? {
                inputs.push(input);
            }
        }
        // Symlinks are deliberately not followed.
    }
    Ok(())
}

/// Converts a GORG ZIP, directory tree, or single FASTA(.gz). Output and
/// sidecars must not exist. Temporary files are removed on conversion errors.
/// Memory is bounded by one FASTA contig plus compression buffers and the index.
pub fn convert(input: &Path, output: &Path, options: &Options) -> Result<Stats> {
    ensure!(options.phred <= 93, "Phred must be 0..=93");
    ensure!(
        options.fragment_length.get() <= 30_000,
        "fragment length must be <=30,000 so each TIRP record fits in one BGZF block"
    );
    ensure!(
        options.kmer_length <= options.fragment_length,
        "fragment length must be at least k"
    );
    ensure!(
        output.to_string_lossy().ends_with(".tirp.gz"),
        "output must end in .tirp.gz"
    );
    let paths = [
        output.to_path_buf(),
        sidecar(output, ".tbi"),
        sidecar(output, ".records.tsv"),
    ];
    for path in &paths {
        ensure!(
            !path.try_exists()?,
            "refusing to overwrite {}",
            path.display()
        );
    }
    let mut archive = None;
    let mut inputs = Vec::new();
    if input.is_dir() {
        collect_directory(input, &mut inputs)?;
    } else if input.extension().is_some_and(|e| e == "zip") {
        let mut zip = ZipArchive::new(File::open(input)?)?;
        for i in 0..zip.len() {
            let entry = zip.by_index(i)?;
            if entry.is_file()
                && let Some(input) = input_name(entry.name(), Source::Zip(i))?
            {
                inputs.push(input);
            }
        }
        archive = Some(zip);
    } else {
        let name = input
            .file_name()
            .and_then(|s| s.to_str())
            .context("invalid input filename")?;
        inputs.push(
            input_name(name, Source::File(input.to_path_buf()))?
                .context("expected FASTA(.gz) or ZIP")?,
        );
    }
    ensure!(!inputs.is_empty(), "no SAG FASTA files found");
    inputs.sort_unstable_by(|a, b| a.cell.as_bytes().cmp(b.cell.as_bytes()));
    for pair in inputs.windows(2) {
        ensure!(
            pair[0].cell != pair[1].cell,
            "duplicate cell ID: {}",
            pair[0].cell
        );
    }

    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut data_tmp = NamedTempFile::new_in(parent)?;
    let mut index_tmp = NamedTempFile::new_in(parent)?;
    let mut map_tmp = NamedTempFile::new_in(parent)?;
    let mut writer = bgzf::Writer::new(BufWriter::new(data_tmp.as_file_mut()));
    let mut map = BufWriter::new(map_tmp.as_file_mut());
    writeln!(
        map,
        "cell_id\trecord_index\tcontig_index\tstart_0based\tend_exclusive\tfasta_header_escaped"
    )?;
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(Builder::bed().build());
    let mut stats = Stats::default();
    for entry in inputs {
        writer.set_cell()?;
        let reader: Box<dyn Read + '_> = match entry.source {
            Source::File(path) => {
                Box::new(File::open(&path).with_context(|| format!("open {}", path.display()))?)
            }
            Source::Zip(i) => Box::new(archive.as_mut().context("missing ZIP")?.by_index(i)?),
        };
        let reader: Box<dyn Read + '_> = if entry.gzip {
            Box::new(MultiGzDecoder::new(reader))
        } else {
            reader
        };
        let mut contig_index = 0u64;
        let mut record_index = 0u64;
        parse_fasta(BufReader::new(reader), |header, sequence| {
            contig_index += 1;
            stats.contigs += 1;
            stats.source_bases += sequence.len() as u64;
            let size = options.fragment_length.get();
            let step = size - options.kmer_length.get() + 1;
            let mut start: usize = 0;
            loop {
                let stop = sequence.len().min(start.saturating_add(size));
                let bases = &sequence[start..stop];
                record_index += 1;
                // id + seven tabs + two coordinate bytes + two base-length
                // fields (R1/Q1) + newline. Never split a record across blocks.
                writer.reserve_record(entry.cell.len() + 10 + 2 * bases.len())?;
                let begin = writer.position();
                writer.write_all(entry.cell.as_bytes())?;
                writer.write_all(b"\t1\t1\t")?;
                writer.write_all(bases)?;
                writer.write_all(b"\t\t")?; // empty R2
                let qualities = [options.phred + 33; 8192];
                let mut left = bases.len();
                while left > 0 {
                    let n = left.min(qualities.len());
                    writer.write_all(&qualities[..n])?;
                    left -= n;
                }
                writer.write_all(b"\t\t\n")?; // empty Q2 and UMI; preserve trailing tab
                let end = writer.position();
                // Match Bascet's BED indexer exactly: raw BED start=1,end=1
                // becomes start=2,end=1 in noodles' one-based inclusive API.
                indexer.add_record(
                    &entry.cell,
                    Position::try_from(2)?,
                    Position::MIN,
                    Chunk::new(begin, end),
                )?;
                writeln!(
                    map,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    entry.cell,
                    record_index,
                    contig_index,
                    start,
                    start + bases.len(),
                    escape_header(header)
                )?;
                stats.records += 1;
                stats.bases += bases.len() as u64;
                if stop == sequence.len() {
                    break;
                }
                start += step;
            }
            Ok(())
        })
        .with_context(|| format!("cell {}", entry.cell))?;
        ensure!(contig_index > 0, "cell {} contains no contigs", entry.cell);
        stats.cells += 1;
    }
    writer.finish()?.flush()?;
    map.flush()?;
    drop(map);
    {
        let mut writer = tabix::io::Writer::new(index_tmp.as_file_mut());
        writer.write_index(&indexer.build())?;
        writer.try_finish()?;
    }
    for file in [&data_tmp, &index_tmp, &map_tmp] {
        file.as_file().sync_all()?;
    }
    // Publish data last. persist_noclobber protects against races as well as
    // the preflight checks. Roll back only files created by this invocation.
    index_tmp.persist_noclobber(&paths[1])?;
    if let Err(e) = map_tmp.persist_noclobber(&paths[2]) {
        let _ = fs::remove_file(&paths[1]);
        return Err(e.into());
    }
    if let Err(e) = data_tmp.persist_noclobber(&paths[0]) {
        let _ = fs::remove_file(&paths[1]);
        let _ = fs::remove_file(&paths[2]);
        return Err(e.into());
    }
    Ok(stats)
}

fn escape_header(header: &str) -> String {
    header
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
}

fn parse_fasta(
    reader: impl BufRead,
    mut emit: impl FnMut(&str, &[u8]) -> Result<()>,
) -> Result<()> {
    let mut header: Option<String> = None;
    let mut sequence = Vec::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        if let Some(next) = line.strip_prefix('>') {
            if let Some(previous) = header.take() {
                ensure!(!sequence.is_empty(), "empty FASTA record {previous:?}");
                emit(&previous, &sequence)?;
                sequence.clear();
            }
            ensure!(
                !next.trim().is_empty(),
                "empty FASTA header at line {}",
                line_index + 1
            );
            header = Some(next.to_owned());
        } else {
            ensure!(
                header.is_some(),
                "sequence before FASTA header at line {}",
                line_index + 1
            );
            ensure!(
                line.bytes()
                    .all(|b| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(&b)),
                "invalid DNA sequence at line {} (expected IUPAC DNA without whitespace)",
                line_index + 1
            );
            sequence.extend_from_slice(line.as_bytes());
        }
    }
    if let Some(header) = header {
        ensure!(!sequence.is_empty(), "empty FASTA record {header:?}");
        emit(&header, &sequence)?;
    } else {
        bail!("empty FASTA input");
    }
    Ok(())
}
