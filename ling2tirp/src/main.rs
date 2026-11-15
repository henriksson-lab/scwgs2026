use anyhow::{Context, Result, ensure};
use clap::Parser;
use flate2::{Compress, Compression, FlushCompress, read::MultiGzDecoder};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::atomic::{AtomicUsize, Ordering},
};
use tempfile::NamedTempFile;

#[derive(Parser)]
#[command(args_override_self = true)]
struct Args {
    #[arg(long, default_value = "/husky/henriksson/atrandi/v9_ling")]
    root: PathBuf,
    #[arg(long)]
    manifest: PathBuf,
    #[arg(long, default_value_t = 8)]
    workers: usize,
    #[arg(long, default_value = "4G")]
    sort_memory: String,
    #[arg(long)]
    sort_temp: Option<PathBuf>,
}
#[derive(Clone, Deserialize, Serialize)]
struct Row {
    library: String,
    barcode_a: String,
    role: String,
    filename: String,
    bytes: u64,
    md5: String,
}
type Group = ((String, String), BTreeMap<String, Row>);
type Fastq = BufReader<MultiGzDecoder<BufReader<File>>>;
const EOF: &[u8] = &[
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];
fn suffix(path: &Path, end: &str) -> PathBuf {
    PathBuf::from(format!("{}{end}", path.display()))
}
fn source(root: &Path, row: &Row) -> PathBuf {
    root.join("source").join(&row.library).join(&row.filename)
}
fn digest(path: &Path) -> Result<String> {
    let mut f = File::open(path)?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1024 * 1024];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(format!("{:x}", h.finalize()))
}
fn trim(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
}
fn record(reader: &mut Fastq, fields: &mut [Vec<u8>; 4]) -> Result<bool> {
    for line in fields.iter_mut() {
        line.clear();
    }
    if reader.read_until(b'\n', &mut fields[0])? == 0 {
        return Ok(false);
    }
    for line in fields.iter_mut().skip(1) {
        ensure!(reader.read_until(b'\n', line)? > 0, "Truncated FASTQ");
    }
    for line in fields.iter_mut() {
        trim(line);
    }
    ensure!(
        fields[0].starts_with(b"@") && fields[2].starts_with(b"+"),
        "Malformed FASTQ"
    );
    ensure!(
        fields[1].len() == fields[3].len(),
        "Sequence/quality length differs"
    );
    ensure!(
        !fields[1]
            .iter()
            .chain(&fields[3])
            .any(|b| b.is_ascii_whitespace()),
        "Whitespace in sequence/quality"
    );
    Ok(true)
}
fn read_key(header: &[u8]) -> &[u8] {
    let token = header[1..]
        .split(|b| b.is_ascii_whitespace())
        .next()
        .unwrap();
    let token = token
        .windows(4)
        .position(|s| s == b"rbc:")
        .map_or(token, |i| &token[..i]);
    let pieces: Vec<&[u8]> = token.split(|b| *b == b':').collect();
    if pieces.len() >= 4 {
        let n = pieces.len();
        let coordinate = pieces[n - 4];
        if matches!(pieces[n - 3], b"Y" | b"N")
            && pieces[n - 2].iter().all(u8::is_ascii_digit)
            && matches!(coordinate.last(), Some(b'1' | b'2'))
        {
            let end =
                token.len() - pieces[n - 1].len() - pieces[n - 2].len() - pieces[n - 3].len() - 4;
            return &token[..end];
        }
    }
    token
        .strip_suffix(b"/1")
        .or_else(|| token.strip_suffix(b"/2"))
        .unwrap_or(token)
}
struct Bbgz<'a> {
    out: BufWriter<&'a mut File>,
    cell: Vec<u8>,
    raw: Vec<u8>,
}
impl<'a> Bbgz<'a> {
    fn new(file: &'a mut File) -> Self {
        Self {
            out: BufWriter::new(file),
            cell: Vec::new(),
            raw: Vec::new(),
        }
    }
    fn record(&mut self, cell: &[u8], line: &[u8]) -> Result<()> {
        ensure!(
            !cell.is_empty() && cell.len() <= 1024 && line.len() <= 60 * 1024,
            "BBGZ record/cell exceeds limit"
        );
        if cell != self.cell || self.raw.len() + line.len() > 60 * 1024 {
            self.flush_block()?;
            self.cell.clear();
            self.cell.extend_from_slice(cell);
        }
        self.raw.extend_from_slice(line);
        Ok(())
    }
    fn flush_block(&mut self) -> Result<()> {
        if self.raw.is_empty() {
            return Ok(());
        }
        let mut compressor = Compress::new(Compression::new(6), false);
        let mut packed = Vec::with_capacity(65536);
        compressor.compress_vec(&self.raw, &mut packed, FlushCompress::Sync)?;
        ensure!(
            compressor.total_in() == self.raw.len() as u64,
            "Compression input not consumed"
        );
        packed.extend_from_slice(&[3, 0]);
        let xlen = 4 + self.cell.len() + 6;
        let size = 12 + xlen + packed.len() + 8;
        ensure!(size <= 65536, "BBGZ block exceeds limit");
        self.out.write_all(&[31, 139, 8, 5, 0, 0, 0, 0, 0, 255])?;
        self.out.write_all(&(xlen as u16).to_le_bytes())?;
        self.out.write_all(b"ID")?;
        self.out
            .write_all(&(self.cell.len() as u16).to_le_bytes())?;
        self.out.write_all(&self.cell)?;
        self.out.write_all(b"BC\x02\x00")?;
        self.out.write_all(&((size - 1) as u16).to_le_bytes())?;
        self.out.write_all(&packed)?;
        self.out
            .write_all(&crc32fast::hash(&self.raw).to_le_bytes())?;
        self.out.write_all(&(self.raw.len() as u32).to_le_bytes())?;
        self.raw.clear();
        Ok(())
    }
    fn finish(mut self) -> Result<()> {
        self.flush_block()?;
        self.out.write_all(EOF)?;
        self.out.flush()?;
        Ok(())
    }
}
fn identity(group: &Group) -> Value {
    let sources: BTreeMap<_, _> = group
        .1
        .iter()
        .map(|(role, row)| {
            (
                role,
                json!({"filename":row.filename,"md5":row.md5,"bytes":row.bytes}),
            )
        })
        .collect();
    json!({"library":group.0.0,"barcode_a":group.0.1,"sources":sources})
}
fn import(root: &Path, index: usize, group: &Group, args: &Args) -> Result<Value> {
    let output = root.join(format!("debarcoded.{index}.tirp.gz"));
    let hist = suffix(&output, ".hist");
    let receipt = suffix(&output, ".json");
    let mut id = identity(group);
    if output.exists() && hist.exists() && receipt.exists() {
        let saved: Value = serde_json::from_reader(File::open(&receipt)?)?;
        ensure!(
            ["library", "barcode_a", "sources"]
                .iter()
                .all(|k| saved[k] == id[k])
                && saved["tirp_bytes"].as_u64() == Some(output.metadata()?.len())
                && saved["hist_bytes"].as_u64() == Some(hist.metadata()?.len()),
            "Existing output/receipt mismatch"
        );
        println!("Reusing {}", output.display());
        return Ok(saved);
    }
    ensure!(
        !output.exists() && !hist.exists() && !receipt.exists(),
        "Refusing incomplete output {}",
        output.display()
    );
    let sortdir = args
        .sort_temp
        .clone()
        .unwrap_or_else(|| root.join("temp/ling_import_rust"))
        .join(format!("{}_{}", group.0.0, group.0.1));
    fs::create_dir_all(&sortdir)?;
    let mut child = Command::new("sort")
        .args([
            "--stable",
            "--field-separator=\t",
            "--key=1,1",
            "--parallel=1",
        ])
        .arg(format!("--buffer-size={}", args.sort_memory))
        .arg(format!("--temporary-directory={}", sortdir.display()))
        .env("LC_ALL", "C")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    let result = (|| -> Result<Value> {
        println!("Importing shard {index}: {} {}", group.0.0, group.0.1);
        let mut readers: Vec<Fastq> = ["r1", "barcode", "r2_fragment"]
            .iter()
            .map(|role| -> Result<Fastq> {
                Ok(BufReader::with_capacity(
                    1024 * 1024,
                    MultiGzDecoder::new(BufReader::with_capacity(
                        1024 * 1024,
                        File::open(source(root, &group.1[*role]))?,
                    )),
                ))
            })
            .collect::<Result<_>>()?;
        let mut records: Vec<[Vec<u8>; 4]> = (0..3)
            .map(|_| std::array::from_fn(|_| Vec::new()))
            .collect();
        let prefix = format!("{}_{}_", group.0.0, group.0.1).into_bytes();
        let mut total = 0u64;
        let mut input = BufWriter::with_capacity(1024 * 1024, child.stdin.take().unwrap());
        let mut line = Vec::with_capacity(1024);
        loop {
            let present: Vec<bool> = readers
                .iter_mut()
                .zip(&mut records)
                .map(|(reader, fields)| record(reader, fields))
                .collect::<Result<_>>()?;
            if present.iter().all(|p| !*p) {
                break;
            }
            ensure!(present.iter().all(|p| *p), "Triplet lengths differ");
            ensure!(
                read_key(&records[0][0]) == read_key(&records[1][0])
                    && read_key(&records[0][0]) == read_key(&records[2][0]),
                "Triplet IDs differ at {}",
                total + 1
            );
            let barcode = &records[1][1];
            ensure!(
                barcode.len() == 24
                    && barcode
                        .iter()
                        .all(|b| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(b)),
                "Invalid 24bp barcode"
            );
            line.clear();
            line.extend_from_slice(&prefix);
            line.extend_from_slice(barcode);
            line.extend_from_slice(b"\t1\t1\t");
            for field in [
                &records[0][1],
                &records[2][1],
                &records[0][3],
                &records[2][3],
            ] {
                line.extend_from_slice(field);
                line.push(b'\t');
            }
            line.push(b'\n');
            input.write_all(&line)?;
            total += 1;
            if total % 5_000_000 == 0 {
                println!("Shard {index}: {total} input read pairs");
            }
        }
        input.flush()?;
        drop(input);
        let mut data = NamedTempFile::new_in(root)?;
        let mut histogram = NamedTempFile::new_in(root)?;
        let mut writer = Bbgz::new(data.as_file_mut());
        let mut histout = BufWriter::new(histogram.as_file_mut());
        let mut sorted = BufReader::with_capacity(1024 * 1024, child.stdout.take().unwrap());
        let mut cell = Vec::new();
        let mut count = 0u64;
        let mut cells = 0u64;
        let mut out_total = 0u64;
        loop {
            line.clear();
            if sorted.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            ensure!(
                line.ends_with(b"\n") && line.iter().filter(|b| **b == b'\t').count() == 7,
                "Malformed sorted TIRP"
            );
            let next = &line[..line.iter().position(|b| *b == b'\t').unwrap()];
            ensure!(next >= cell.as_slice(), "Sort is not monotonic");
            if next != cell {
                if count > 0 {
                    histout.write_all(&cell)?;
                    writeln!(histout, "\t{count}")?;
                }
                cell.clear();
                cell.extend_from_slice(next);
                count = 0;
                cells += 1;
            }
            writer.record(next, &line)?;
            count += 1;
            out_total += 1;
        }
        if count > 0 {
            histout.write_all(&cell)?;
            writeln!(histout, "\t{count}")?;
        }
        histout.flush()?;
        drop(histout);
        writer.finish()?;
        ensure!(
            child.wait()?.success() && total == out_total && total > 0,
            "Sort/import failed"
        );
        data.as_file().sync_all()?;
        histogram.as_file().sync_all()?;
        id["records"] = json!(total);
        id["cells"] = json!(cells);
        id["tirp_bytes"] = json!(data.as_file().metadata()?.len());
        id["hist_bytes"] = json!(histogram.as_file().metadata()?.len());
        id["tirp_sha256"] = json!(digest(data.path())?);
        id["hist_sha256"] = json!(digest(histogram.path())?);
        fs::hard_link(data.path(), &output)?;
        fs::hard_link(histogram.path(), &hist)?;
        let mut marker = NamedTempFile::new_in(root)?;
        serde_json::to_writer_pretty(&mut marker, &id)?;
        marker.flush()?;
        fs::hard_link(marker.path(), receipt)?;
        println!(
            "Ready {}: {total} read pairs, {cells} observed cells",
            output.display()
        );
        Ok(id)
    })();
    if result.is_err() {
        let _ = child.kill();
        let _ = child.wait();
    }
    result
}
fn run(args: Args) -> Result<()> {
    ensure!(args.workers > 0, "workers must be positive");
    fs::create_dir_all(&args.root)?;
    let lock = File::options()
        .create(true)
        .truncate(false)
        .write(true)
        .open(args.root.join(".import.lock"))?;
    fs2::FileExt::try_lock_exclusive(&lock).context("Another Ling import is running")?;
    let mut groups: BTreeMap<(String, String), BTreeMap<String, Row>> = BTreeMap::new();
    for row in csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.manifest)?
        .deserialize::<Row>()
    {
        let row = row?;
        ensure!(
            Path::new(&row.filename)
                .file_name()
                .and_then(|s| s.to_str())
                == Some(&row.filename),
            "Unsafe filename"
        );
        ensure!(
            row.library
                .bytes()
                .all(|b| b.is_ascii_lowercase() || b == b'_')
                && row.barcode_a.len() == 8
                && row.barcode_a.bytes().all(|b| b"ACGT".contains(&b)),
            "Unsafe library/barcode"
        );
        let path = source(&args.root, &row);
        let saved: Value = serde_json::from_reader(File::open(suffix(&path, ".verified.json"))?)?;
        ensure!(
            path.metadata()?.len() == row.bytes
                && saved["bytes"].as_u64() == Some(row.bytes)
                && saved["md5"].as_str() == Some(&row.md5),
            "Source/receipt mismatch"
        );
        let group = groups
            .entry((row.library.clone(), row.barcode_a.clone()))
            .or_default();
        ensure!(!group.contains_key(&row.role), "Duplicate source role");
        group.insert(row.role.clone(), row);
    }
    let groups: Vec<Group> = groups.into_iter().collect();
    ensure!(!groups.is_empty(), "No source groups");
    for (_, rows) in &groups {
        ensure!(
            rows.len() == 3
                && ["r1", "barcode", "r2_fragment"]
                    .iter()
                    .all(|r| rows.contains_key(*r)),
            "Incomplete triplet"
        );
    }
    let cursor = AtomicUsize::new(0);
    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..args.workers.min(groups.len()))
            .map(|_| {
                scope.spawn(|| -> Result<Vec<(usize, Value)>> {
                    let mut done = Vec::new();
                    loop {
                        let i = cursor.fetch_add(1, Ordering::Relaxed);
                        if i >= groups.len() {
                            break;
                        }
                        done.push((
                            i,
                            import(&args.root, i + 1, &groups[i], &args)
                                .with_context(|| format!("Group {:?}", groups[i].0))?,
                        ));
                    }
                    Ok(done)
                })
            })
            .collect();
        let mut done = Vec::new();
        for handle in handles {
            done.extend(
                handle
                    .join()
                    .map_err(|_| anyhow::anyhow!("Import worker panicked"))??,
            );
        }
        Ok::<_, anyhow::Error>(done)
    })?;
    let results: BTreeMap<_, _> = results.into_iter().collect();
    let mut metadata = NamedTempFile::new_in(&args.root)?;
    {
        let mut csv = csv::Writer::from_writer(metadata.as_file_mut());
        csv.write_record([
            "prefix",
            "r1",
            "r2",
            "dir",
            "barcode_a",
            "barcode_source",
            "shard",
            "records",
            "cells",
        ])?;
        for (i, (group, rows)) in groups.iter().enumerate() {
            let result = &results[&i];
            csv.write_record([
                group.0.clone(),
                rows["r1"].filename.clone(),
                rows["r2_fragment"].filename.clone(),
                source(&args.root, &rows["r1"])
                    .parent()
                    .unwrap()
                    .display()
                    .to_string(),
                group.1.clone(),
                rows["barcode"].filename.clone(),
                (i + 1).to_string(),
                result["records"].to_string(),
                result["cells"].to_string(),
            ])?;
        }
        csv.flush()?;
    }
    let output = args.root.join("debarcoded.meta");
    if output.exists() {
        ensure!(
            fs::read(&output)? == fs::read(metadata.path())?,
            "Metadata differs from existing file"
        );
    } else {
        fs::hard_link(metadata.path(), output)?;
    }
    println!("Ready {} sorted TIRPs and debarcoded.meta", groups.len());
    Ok(())
}
fn main() -> Result<()> {
    run(Args::parse())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn collapsed_mate_headers() {
        assert_eq!(
            read_key(b"@NB:1:10491:N:0:0rbc:ACTG"),
            read_key(b"@NB:1:10492:N:0:0rbc:ACTG")
        );
        assert_eq!(read_key(b"@read1/1"), b"read1");
    }
    #[test]
    fn bbgz_roundtrip() -> Result<()> {
        let mut f = NamedTempFile::new()?;
        let line = b"cell\t1\t1\tACGT\tTGCA\tIIII\tJJJJ\t\n";
        {
            let mut w = Bbgz::new(f.as_file_mut());
            for _ in 0..4000 {
                w.record(b"cell", line)?;
            }
            w.finish()?;
        }
        let mut raw = Vec::new();
        MultiGzDecoder::new(File::open(f.path())?).read_to_end(&mut raw)?;
        assert_eq!(raw, line.repeat(4000));
        Ok(())
    }
}
