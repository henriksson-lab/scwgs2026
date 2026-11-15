use anyhow::{Context, Result, ensure};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
struct Args {
    #[arg(long)]
    r1: PathBuf,
    #[arg(long)]
    r2: PathBuf,
    #[arg(long)]
    tirp: PathBuf,
    #[arg(long)]
    hist: PathBuf,
    #[arg(long)]
    barcodes: PathBuf,
    #[arg(long)]
    report: PathBuf,
    #[arg(long)]
    library: String,
}

type Table = HashMap<[u8; 7], (u8, Option<String>)>;

fn tables(path: &PathBuf) -> Result<[Table; 3]> {
    let mut tables: [Table; 3] = Default::default();
    let mut counts = [0; 3];
    for (i, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if i == 0 {
            ensure!(line == "pos\twell\tseq", "Unexpected whitelist header");
            continue;
        }
        let columns: Vec<_> = line.split('\t').collect();
        ensure!(columns.len() == 3, "Malformed whitelist row");
        let round: usize = columns[0].parse()?;
        ensure!((1..=3).contains(&round), "Invalid barcode round");
        let sequence: [u8; 7] = columns[2]
            .as_bytes()
            .try_into()
            .context("Barcode must have seven bases")?;
        counts[round - 1] += 1;
        for position in 0..7 {
            for base in b"ACGTN" {
                let mut variant = sequence;
                variant[position] = *base;
                let distance = u8::from(*base != sequence[position]);
                let candidate = Some(columns[1].to_owned());
                tables[round - 1]
                    .entry(variant)
                    .and_modify(|previous| {
                        if distance < previous.0 {
                            *previous = (distance, candidate.clone());
                        } else if distance == previous.0 && candidate != previous.1 {
                            previous.1 = None;
                        }
                    })
                    .or_insert((distance, candidate));
            }
        }
    }
    ensure!(
        counts == [96; 3],
        "Expected 96 whitelist entries in each round"
    );
    Ok(tables)
}

fn decode(sequence: &[u8], tables: &[Table; 3], prefix: &str) -> Option<Vec<u8>> {
    if sequence.len() < 57 {
        return None;
    }
    let sequence: Vec<_> = sequence.iter().map(u8::to_ascii_uppercase).collect();
    for (start, expected) in [
        (14, b"GGTCCTTGGCTTCGC".as_slice()),
        (36, b"CCTCCTACGCCAGA".as_slice()),
    ] {
        let differences = sequence[start..start + expected.len()]
            .iter()
            .zip(expected)
            .filter(|(a, b)| a != b)
            .count();
        if differences > 2 {
            return None;
        }
    }
    let mut names = Vec::new();
    for (round, start) in [50, 29, 7].into_iter().enumerate() {
        let key: [u8; 7] = sequence[start..start + 7].try_into().ok()?;
        names.push(tables[round].get(&key)?.1.as_ref()?.as_str());
    }
    Some(format!("{prefix}_{}_{}_{}", names[0], names[1], names[2]).into_bytes())
}

fn trim_end(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\r' | b'\n')) {
        line.pop();
    }
}

#[derive(Default)]
struct Record {
    name: Vec<u8>,
    sequence: Vec<u8>,
    plus: Vec<u8>,
    quality: Vec<u8>,
}
impl Record {
    fn read(&mut self, reader: &mut impl BufRead) -> Result<bool> {
        for line in [
            &mut self.name,
            &mut self.sequence,
            &mut self.plus,
            &mut self.quality,
        ] {
            line.clear();
        }
        if reader.read_until(b'\n', &mut self.name)? == 0 {
            return Ok(false);
        }
        for line in [&mut self.sequence, &mut self.plus, &mut self.quality] {
            ensure!(reader.read_until(b'\n', line)? > 0, "Truncated FASTQ");
            trim_end(line);
        }
        ensure!(
            self.name.starts_with(b"@") && self.plus.starts_with(b"+"),
            "Malformed FASTQ"
        );
        ensure!(
            !self.sequence.is_empty() && self.sequence.len() == self.quality.len(),
            "Invalid FASTQ sequence/quality length"
        );
        if let Some(end) = self.name.iter().position(u8::is_ascii_whitespace) {
            self.name.truncate(end);
        }
        Ok(true)
    }
}

fn gz(path: &PathBuf) -> Result<BufReader<MultiGzDecoder<File>>> {
    Ok(BufReader::with_capacity(
        1 << 20,
        MultiGzDecoder::new(File::open(path)?),
    ))
}

fn add_digest(sum: &mut [u8; 32], line: &[u8]) {
    let digest = Sha256::digest(line);
    let mut carry = 0u16;
    for i in (0..32).rev() {
        let value = u16::from(sum[i]) + u16::from(digest[i]) + carry;
        sum[i] = value as u8;
        carry = value >> 8;
    }
}

fn verify(args: Args) -> Result<()> {
    let started = Instant::now();
    let tables = tables(&args.barcodes)?;
    let (mut left, mut right) = (gz(&args.r1)?, gz(&args.r2)?);
    let (mut r1, mut r2) = (Record::default(), Record::default());
    let mut expected: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut expected_digest = [0u8; 32];
    let (mut total, mut accepted) = (0u64, 0u64);
    let mut line = Vec::new();
    loop {
        let (a, b) = (r1.read(&mut left)?, r2.read(&mut right)?);
        if !a && !b {
            break;
        }
        ensure!(
            a && b && r1.name == r2.name,
            "Mates disagree after {total} pairs"
        );
        total += 1;
        if let Some(cell) = decode(&r1.sequence, &tables, &args.library) {
            accepted += 1;
            line.clear();
            line.extend_from_slice(&cell);
            line.extend_from_slice(b"\t1\t1\t\t");
            line.extend_from_slice(&r2.sequence);
            line.extend_from_slice(b"\t\t");
            line.extend_from_slice(&r2.quality);
            line.push(b'\t');
            line.extend_from_slice(&r1.sequence[..7]);
            add_digest(&mut expected_digest, &line);
            *expected.entry(cell).or_default() += 1;
        }
        if total.is_multiple_of(1_000_000) {
            eprintln!("Oracle: {total} pairs, {accepted} accepted");
        }
    }
    ensure!(accepted > 0, "No accepted input pairs");
    let mut observed: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut observed_digest = [0u8; 32];
    let mut input = gz(&args.tirp)?;
    let mut previous = Vec::new();
    let mut records = 0u64;
    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        trim_end(&mut line);
        let fields: Vec<_> = line.split(|b| *b == b'\t').collect();
        ensure!(fields.len() == 8, "TIRP record must have eight columns");
        ensure!(
            fields[0] >= previous.as_slice(),
            "TIRP is not cell-sorted at record {}",
            records + 1
        );
        previous.clear();
        previous.extend_from_slice(fields[0]);
        *observed.entry(fields[0].to_vec()).or_default() += 1;
        add_digest(&mut observed_digest, &line);
        records += 1;
        if records.is_multiple_of(1_000_000) {
            eprintln!("TIRP: {records} records checked");
        }
    }
    ensure!(
        records == accepted && observed_digest == expected_digest,
        "Output read multiset differs from independently decoded input"
    );
    ensure!(
        observed == expected,
        "Output cell counts differ from independently decoded input"
    );
    let mut histogram = HashMap::new();
    for row in BufReader::new(File::open(&args.hist)?).lines() {
        let row = row?;
        let columns: Vec<_> = row.split('\t').collect();
        ensure!(columns.len() == 2, "Malformed histogram row");
        ensure!(
            histogram
                .insert(columns[0].as_bytes().to_vec(), columns[1].parse::<u64>()?)
                .is_none(),
            "Duplicate histogram cell"
        );
    }
    ensure!(
        histogram == observed,
        "Histogram differs from actual TIRP records"
    );
    let digest: String = expected_digest.iter().map(|b| format!("{b:02x}")).collect();
    let report = serde_json::json!({"status":"passed","library":args.library,
      "input_pairs":total,"accepted_pairs":accepted,"rejected_pairs":total-accepted,
      "output_records":records,"barcode_combinations":observed.len(),"sorted_by_cell":true,
      "all_records_checked":true,"histogram_matches":true,"multiset_sha256_sum":digest,
      "multiset_contents":"All eight TIRP columns; original R2 sequence/quality and seven-base R1 UMI",
      "seconds":started.elapsed().as_secs_f64(),"r1":args.r1,"r2":args.r2,"tirp":args.tirp,
      "implementation":"independent Rust streaming oracle"});
    let text = serde_json::to_string_pretty(&report)?;
    fs::write(&args.report, format!("{text}\n"))?;
    println!("{text}");
    Ok(())
}
fn main() -> Result<()> {
    verify(Args::parse())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn multiset_sum_wraps_as_big_endian_256_bit_integer() {
        let mut sum = [255u8; 32];
        let digest = Sha256::digest(b"record");
        add_digest(&mut sum, b"record");
        let mut expected: [u8; 32] = digest.into();
        for b in expected.iter_mut().rev() {
            let (v, borrow) = b.overflowing_sub(1);
            *b = v;
            if !borrow {
                break;
            }
        }
        assert_eq!(sum, expected);
        let mut a = [0; 32];
        let mut b = [0; 32];
        for line in [b"first".as_slice(), b"second"] {
            add_digest(&mut a, line);
        }
        for line in [b"second".as_slice(), b"first"] {
            add_digest(&mut b, line);
        }
        assert_eq!(a, b);
        add_digest(&mut b, b"first");
        assert_ne!(a, b);
    }
    #[test]
    fn fastq_rejects_truncated_or_wrong_quality_length() {
        let mut record = Record::default();
        assert!(record.read(&mut b"@id\nACGT\n+\nIII\n".as_slice()).is_err());
        assert!(record.read(&mut b"@id\nACGT\n+\n".as_slice()).is_err());
        assert!(
            record
                .read(&mut b"@id annotation\nACGT\n+\nIIII\n".as_slice())
                .unwrap()
        );
        assert_eq!(record.name, b"@id");
    }
}
