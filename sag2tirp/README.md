# sag2tirp

Convert the [GORG-Dark SAG assemblies](https://figshare.com/articles/dataset/GORG-Dark_SAGs/29087090)
to a sorted, indexed Bascet TIRP for **31-mer counting**.

For the project dataset, run the repository's executable script:

    ./run_gorg_sag2tirp.sh

It downloads and checksum-verifies the Figshare archive, then writes
/husky/henriksson/atrandi/v9_wgs_gorg/filtered_nohost.1.tirp.gz and its
sidecars. GORG_DIR overrides the destination directory; GORG_INPUT selects
an existing ZIP/directory/FASTA instead of downloading. The filename follows
the project convention; the converter does not perform host filtering.

After conversion, prepare count/sketch inputs for downstream UMAP analysis:

    Rscript /husky/henriksson/atrandi/revision_scripts/prep_gorg.R

This writes `countsketch_mat.feather` (31-mer CountSketch),
`kraken_counts.N.h5` (Kraken taxon counts), and `kmer_counts.N.h5`
(informative 31-mer counts). The Zorn readers are
`BascetLoadCountSketchMatrix(root)`,
`ReadBascetCountMatrix(root, "kraken_counts")`, and
`ReadBascetCountMatrix(root, "kmer_counts")` respectively.
GORG cells are the deposited SAG assemblies; they are not reassembled.
FASTQC, SKESA, GECCO, QUAST, reference alignment, and host filtering are
outside this prep workflow. `GORG_RUN_KRAKEN=0` skips Kraken.

UMAP computation remains a separate analysis step. `post_unified.R` contains
the existing CountSketchUMAP and count-matrix TF-IDF/SVD/RunUMAP examples,
but its mock-data paths and filters require adaptation before using GORG.

For other paths:

    cargo run --release --manifest-path sag2tirp/Cargo.toml -- \
      --input "/path/to/GORG-Dark SAGs.zip" \
      --output gorg-dark.tirp.gz

Input can also be a directory tree of per-cell FASTAs, or one FASTA.
Supported suffixes: .fa, .fna, .fasta, optionally followed by .gz.
ZIP entries are read directly without extracting the archive. macOS
__MACOSX and AppleDouble entries are ignored.

## Fragments and k-mer counts

Defaults: **150-base fragments, k=31, overlap=30, step=120**.
For a contig of length 300, the zero-based, half-open ranges are:

    [0,150)  [120,270)  [240,300)

The first fragment contributes 31-mer starts 0–119, the second 120–239,
and the third 240–269. Every original 31-mer occurrence appears **exactly
once**, including those across fragment boundaries. Repeated genomic
31-mers retain their original multiplicity. This also holds for canonical
31-mers.

The last fragment ends at the contig end. No redundant overlap-only tail
is emitted. A contig shorter than 31 bases is retained as one short record
and contributes no 31-mers. Fragments never cross contig or cell boundaries.
IUPAC ambiguity symbols and sequence case are preserved; interpretation
of ambiguous k-mers remains the downstream tool's responsibility.

Change the fragment length with --fragment-length 250. Change the target k
with --kmer-length K; overlap is always K-1. Count preservation is for that
exact K: using a smaller downstream k can double-count k-mers inside overlaps.
Fragment length must be at least K and no more than 30,000.

These are assembly fragments, not original sequencing reads: they do not
reproduce sequencing depth, pairing, or measured qualities.

## Output layout

Each FASTA file becomes one cell. For example,
SCGC_AH-573-E06_contigs.fasta.gz becomes SCGC_AH-573-E06.
Only filenames are sorted, in ascending byte order, before reading their
contents. All fragments for one cell are written together. No external
sorting, sharding, or merge step is performed. Duplicate cell IDs are errors.

The output has **exactly eight tab-separated columns**, no header, and a
newline after every record:

| Column | Value |
|---|---|
| 1 | Cell ID |
| 2 | Literal 1 |
| 3 | Literal 1 |
| 4 | Fragment sequence (R1) |
| 5 | Empty R2 |
| 6 | Synthetic Phred+33 qualities, same length as R1 |
| 7 | Empty Q2 |
| 8 | Empty UMI |

The default synthetic score is Q0 (!). --phred 40 instead writes I.
These are placeholders; configure downstream quality filtering accordingly.
The empty final UMI column is retained with a trailing tab.
FASTA headers are not put into UMI or cell ID fields.

Compression is **standard BGZF**, using noodles, matching Bascet's indexed
TIRP writer. Every record fits wholly inside one block; blocks are flushed
at cell boundaries. Raw blocks are limited to 60 KiB. This avoids long
records even for very long contigs. The file ends with the standard BGZF
EOF marker. No BBGZ cell-ID extension is written.

Three files are produced:

- gorg-dark.tirp.gz: fragment records.
- gorg-dark.tirp.gz.tbi: tabix BED index, created in Rust. The coordinates
  1,1 are dummy Bascet coordinates, not contig positions. Query by cell:
  tabix gorg-dark.tirp.gz SCGC_AH-573-E06.
- gorg-dark.tirp.gz.records.tsv: source contig header and zero-based,
  half-open fragment coordinates. Record and contig indices are one-based
  within each cell. Headers escape backslash as \\, tab as \t, and carriage
  return as \r. This permits reconstruction without a ninth TIRP field.

Existing output files and sidecars are never overwritten. All three are
staged in temporary files; invalid or truncated inputs leave no published
output. Memory usage scales with one contig, the file list, compression
buffers, and the index. The full dataset is not held in memory.
Source-base counts exclude overlap; emitted-base counts include it.

## Library API

    use std::path::Path;
    use sag2tirp::{convert, Options};

    let stats = convert(
        Path::new("GORG-Dark SAGs.zip"),
        Path::new("gorg-dark.tirp.gz"),
        &Options::default(),
    )?;

## Format references and validation

Checked against the [published TIRP layout](https://zorn.henlab.org/articles/fileformats.html)
and the local Bascet checkout at commit
c378038b1a956c3e7312f680329340c9362fd693:
fileformat/tirp.rs, utils/tabix_bed.rs, and bascet-io/src/parse/tirp/.
The index uses Bascet's dummy-coordinate conversion.

Tests check exact k-mer multisets at fragment/tail boundaries, long-contig
reconstruction, ZIP/gzip handling, sorting, empty fields, malformed inputs,
and truncated gzip cleanup. Local interoperability checks passed with
bgzip -t, tabix cell queries, and Bascet 1.0.0-beta.3 transform and
countsketch -k 31. A 500,000-base contig retained its exact 31-mer multiset
after conversion through Bascet to FASTQ.

    cargo test --manifest-path sag2tirp/Cargo.toml
    cargo clippy --manifest-path sag2tirp/Cargo.toml --all-targets -- -D warnings
