use std::{
    fs,
    io::{Read, Write},
    num::NonZeroUsize,
};

use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use noodles::csi::BinningIndex;
use sag2tirp::{Options, convert, sidecar};
use tempfile::tempdir;
use zip::{ZipWriter, write::SimpleFileOptions};

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

fn decode(path: &std::path::Path) -> String {
    let mut result = String::new();
    MultiGzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut result)
        .unwrap();
    result
}

#[test]
fn sorted_zip_preserves_contigs_and_skips_appledouble() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("input.zip");
    let output = dir.path().join("output.tirp.gz");
    let mut zip = ZipWriter::new(fs::File::create(&input).unwrap());
    for (name, bytes) in [
        (
            "GORG-Dark SAGs/SCGC_Z_contigs.fasta.gz",
            gzip(b">z\r\nAC\r\nGT\r\n"),
        ),
        (
            "__MACOSX/GORG-Dark SAGs/._SCGC_A_contigs.fasta.gz",
            b"not gzip".to_vec(),
        ),
        (
            "GORG-Dark SAGs/SCGC_A_contigs.fasta.gz",
            gzip(b">c1 with\ttab\\slash\nacN\n>c2\nRY"),
        ),
    ] {
        zip.start_file(name, SimpleFileOptions::default()).unwrap();
        zip.write_all(&bytes).unwrap();
    }
    zip.finish().unwrap();
    let stats = convert(&input, &output, &Options::default()).unwrap();
    assert_eq!(
        (stats.cells, stats.contigs, stats.records, stats.bases),
        (2, 3, 3, 9)
    );
    assert_eq!(
        decode(&output),
        "SCGC_A\t1\t1\tacN\t\t!!!\t\t\nSCGC_A\t1\t1\tRY\t\t!!\t\t\nSCGC_Z\t1\t1\tACGT\t\t!!!!\t\t\n"
    );
    let index = noodles::tabix::fs::read(sidecar(&output, ".tbi")).unwrap();
    let names: Vec<_> = index
        .header()
        .unwrap()
        .reference_sequence_names()
        .iter()
        .collect();
    assert_eq!(names.len(), 2);
    assert_eq!(&names[0][..], b"SCGC_A");
    let map = fs::read_to_string(sidecar(&output, ".records.tsv")).unwrap();
    assert!(map.contains("SCGC_A\t1\t1\t0\t3\tc1 with\\ttab\\\\slash\n"));
    assert!(convert(&input, &output, &Options::default()).is_err());
}

#[test]
fn fragments_retain_tail_and_multimember_gzip() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("cell.fa.gz");
    let output = dir.path().join("output.tirp.gz");
    let mut bytes = gzip(b">contig\nACGT");
    bytes.extend(gzip(b"N\n"));
    fs::write(&input, bytes).unwrap();
    let opts = Options {
        fragment_length: NonZeroUsize::new(3).unwrap(),
        kmer_length: NonZeroUsize::new(2).unwrap(),
        phred: 40,
    };
    let stats = convert(&input, &output, &opts).unwrap();
    assert_eq!(
        (
            stats.contigs,
            stats.records,
            stats.source_bases,
            stats.bases
        ),
        (1, 2, 5, 6)
    );
    assert_eq!(
        decode(&output),
        "cell\t1\t1\tACG\t\tIII\t\t\ncell\t1\t1\tGTN\t\tIII\t\t\n"
    );
    let map = fs::read_to_string(sidecar(&output, ".records.tsv")).unwrap();
    assert!(map.contains("cell\t2\t1\t2\t5\tcontig\n"));
}

#[test]
fn long_contig_is_split_and_reconstructs_exactly() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("long.fna");
    let output = dir.path().join("output.tirp.gz");
    let sequence = "ACGTN".repeat(100_000);
    fs::write(&input, format!(">large\n{sequence}\n")).unwrap();
    let stats = convert(&input, &output, &Options::default()).unwrap();
    assert_eq!(stats.records, 1 + (500_000_u64 - 150).div_ceil(120));
    let text = decode(&output);
    let mut reconstructed = String::new();
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 8);
        assert!(fields[3].len() <= 150);
        assert_eq!(fields[5].len(), fields[3].len());
        assert_eq!(fields[7], "");
        if reconstructed.is_empty() {
            reconstructed.push_str(fields[3]);
        } else {
            reconstructed.push_str(&fields[3][30..]);
        }
    }
    assert_eq!(reconstructed, sequence);

    // Independently inflate each standard BGZF member. No record may span
    // members, even when a long contig produces thousands of fragments.
    let bytes = fs::read(&output).unwrap();
    let mut offset = 0;
    while offset < bytes.len() {
        assert_eq!(&bytes[offset..offset + 4], &[31, 139, 8, 4]);
        assert_eq!(&bytes[offset + 12..offset + 16], b"BC\x02\x00");
        let size = u16::from_le_bytes([bytes[offset + 16], bytes[offset + 17]]) as usize + 1;
        let mut raw = String::new();
        flate2::read::GzDecoder::new(&bytes[offset..offset + size])
            .read_to_string(&mut raw)
            .unwrap();
        assert!(raw.len() <= 60 * 1024);
        assert!(raw.is_empty() || raw.ends_with('\n'));
        for line in raw.lines() {
            assert_eq!(line.split('\t').count(), 8);
            assert!(line.starts_with("long\t"));
        }
        offset += size;
    }
}

#[test]
fn overlapping_fragments_preserve_exact_kmer_multisets_at_boundaries() {
    use std::collections::BTreeMap;
    fn kmers(sequence: &[u8], k: usize, counts: &mut BTreeMap<Vec<u8>, usize>) {
        for kmer in sequence.windows(k) {
            *counts.entry(kmer.to_vec()).or_default() += 1;
        }
    }
    for k in [1, 31] {
        for length in [k, k + 1, 150] {
            for n in [
                1, 30, 31, 32, 119, 120, 149, 150, 151, 239, 240, 241, 269, 270, 271, 1000,
            ] {
                let dir = tempdir().unwrap();
                let input = dir.path().join("cell.fa");
                let output = dir.path().join("out.tirp.gz");
                let mut state = 17u32;
                let sequence: String = (0..n)
                    .map(|_| {
                        state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                        b"ACGT"[(state >> 24) as usize % 4] as char
                    })
                    .collect();
                fs::write(&input, format!(">a\n{sequence}\n")).unwrap();
                let opts = Options {
                    fragment_length: NonZeroUsize::new(length).unwrap(),
                    kmer_length: NonZeroUsize::new(k).unwrap(),
                    phred: 0,
                };
                let stats = convert(&input, &output, &opts).unwrap();
                let mut expected = BTreeMap::new();
                kmers(sequence.as_bytes(), k, &mut expected);
                let mut actual = BTreeMap::new();
                for line in decode(&output).lines() {
                    let fragment = line.split('\t').nth(3).unwrap();
                    assert!(fragment.len() <= length);
                    kmers(fragment.as_bytes(), k, &mut actual);
                }
                assert_eq!(actual, expected, "n={n}, k={k}, fragment_length={length}");
                assert_eq!(stats.source_bases, n as u64);
            }
        }
    }
}

#[test]
fn malformed_inputs_leave_no_outputs_and_duplicates_fail() {
    for fasta in [b"ACGT\n".as_slice(), b">empty\n", b">bad\nAC X\n", b""] {
        let dir = tempdir().unwrap();
        let input = dir.path().join("cell.fa");
        let output = dir.path().join("output.tirp.gz");
        fs::write(&input, fasta).unwrap();
        assert!(convert(&input, &output, &Options::default()).is_err());
        assert!(!output.exists());
        assert!(!sidecar(&output, ".tbi").exists());
        assert!(!sidecar(&output, ".records.tsv").exists());
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
    }
    let dir = tempdir().unwrap();
    fs::write(dir.path().join("cell.fa"), ">a\nAC\n").unwrap();
    fs::write(dir.path().join("cell_contigs.fasta"), ">b\nGT\n").unwrap();
    let error = convert(
        dir.path(),
        &dir.path().join("out.tirp.gz"),
        &Options::default(),
    )
    .unwrap_err();
    assert!(error.to_string().contains("duplicate cell ID"));
}

#[test]
fn truncated_gzip_fails_without_publishing() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("cell.fa.gz");
    let output = dir.path().join("output.tirp.gz");
    let mut bytes = gzip(b">a\nACGT\n");
    bytes.truncate(bytes.len() - 5);
    fs::write(&input, bytes).unwrap();
    assert!(convert(&input, &output, &Options::default()).is_err());
    assert!(!output.exists());
}
