# Yanai / Pountain bacterial PETRI-seq

[Transcription–replication interactions reveal bacterial genome regulation](https://doi.org/10.1038/s41586-023-06974-w),
Pountain et al., Nature 626, 661–669 (2024).
Raw data: GEO GSE217715 / ENA BioProject PRJNA900215.

The pinned manifest maps 88 sequencing runs to 88 GEO library samples,
including dataset, library, organism, raw-mate URLs, byte sizes and MD5s.
There are 176 FASTQs, approximately 146 GB compressed. These are multiplexed
libraries, not per-cell files or files already separated by one barcode.

## Download and prepare Bascet input

From the repository:

    ./run_yanai_download.sh

Default destination: /husky/henriksson/atrandi/v9_rnaseq_yanai.
Override with YANAI_ROOT. YANAI_DOWNLOAD_WORKERS defaults to 8;
YANAI_BGZIP_THREADS defaults to 2 per worker.

The script downloads missing mates resumably with curl, verifies each ENA
MD5 and size, then streams gzip decompression into bgzip. Original read
names, sequences, quality scores and ordering are unchanged. Recompression
is needed because Bascet debarcode requires BGZF rather than ordinary gzip.
Completed outputs have a JSON receipt; rerunning validates and reuses them.
Existing outputs without matching receipts are never overwritten.

Output layout:

    raw_gzip/SRR22254412_1.fastq.gz
    raw_gzip/SRR22254412_2.fastq.gz
    raw_bgzf/D1_lib1_S1_L001_R1_001.fastq.gz
    raw_bgzf/D1_lib1_S1_L001_R2_001.fastq.gz
    raw_bgzf/D1_lib1_S1_L001_R1_001.fastq.gz.json
    rawdata.txt
    fastq_manifest.tsv

The Illumina-style names make Zorn DetectRawFileMeta identify each GEO
library as its own prefix. All files/lanes of one independent library must
use the same prefix; separate libraries must not share one.

For an existing source directory, without downloading:

    Rscript yanai_petriseq/prepare_fastqs.R \
      --root /path/to/dataset --raw-dir /path/to/SRR_fastqs

Use repeated --run SRR... arguments to prepare a subset. The default full
prep_yanai.R workflow expects all libraries listed in its dataset manifest.

## Debarcode and shard with Zorn

To verify a complete library against your local Bascet checkout first:

    Rscript /husky/henriksson/atrandi/revision_scripts/prep_yanai_dev.R

This builds `~/github/bascet` locally with `cargo build --locked --release
--bin bascet`, then uses Zorn's `getBascetDevDir()` and a local runner.
It does not fetch or publish a Bascet release. The default library is
`D1_lib1`; set `YANAI_LIBRARY` to select another manifest library.
`YANAI_BASCET_DIR`, `YANAI_BASCET_TARGET=release|debug`, `YANAI_THREADS`,
and `YANAI_MEM` configure the local run; `YANAI_BUILD=0` skips rebuilding.

Outputs go to `v9_rnaseq_yanai/dev_debarcode/D1_lib1/` by default.
For another validation run, select a fresh `YANAI_OUTPUT_DIR`: existing
outputs are deliberately refused so old results cannot validate a new binary.
`local_bascet.txt` records the exact binary path and checksum.

After debarcoding, an independent Rust decoder checks every input pair
against every output record using an order-independent SHA-256 multiset sum.
This checks multiplicity, corrected cell IDs, UMIs, transcript sequences and
quality scores, as well as lexical sorting and per-cell histogram counts.
The result is saved as `debarcoding_validation.json`. This development script
performs debarcoding and verification only; the full workflow below adds
alignment after local verification succeeds.

The Rust verifier has also checked the complete D1_lib1 library: 44,872,739
input pairs, 36,175,938 accepted records and 359,666 barcode combinations.
Its complete eight-column SHA-256 multiset sum matches the earlier independent
whole-library oracle exactly. The report is retained in
`v9_rnaseq_yanai/validation_full_D1_lib1/debarcoding_validation.rust.json`.

Once all FASTQs and the reference index are ready:

    ./run_yanai_reference_download.sh
    ./run_yanai_reference_index.sh
    Rscript /husky/henriksson/atrandi/revision_scripts/prep_yanai.R

This sources the current Zorn checkout and uses the current Bascet release
binary. Environment overrides are listed at the top of prep_yanai.R.
Build the modified Bascet release before running on another machine.

The script uses chemistry=petriseq, produces debarcoded TIRPs with GEO
library prefixes, and uses Zorn's normal metadata/statistics grouping.
It keeps ALL valid barcode combinations, including single-read combinations,
and writes one filtered.N.tirp.gz shard per independent library.
yanai_libraries.csv and yanai_shard_libraries.csv preserve the mapping.
It then aligns the RNA reads with Bascet BWA-MEM2 against
/husky/henriksson/atrandi/bwamem_ref/ling/all.fa, producing aligned_cell and
aligned_pos BAMs and BAM indexes. It counts reads per cell and gene from the
position-sorted BAMs using the matching `all.gff` annotation (`YANAI_GFF`
overrides this). `gene` features use `locus_tag` as the stable gene ID and
`Name` as the gene name attribute. No RNA expression thresholds, UMI
deduplication, host filtering, or k-mer analysis are performed.

The counting stage explicitly selects Bascet `countfeature --count-mode reads`;
the historical default counts UMI molecules and would collapse repeated UMIs.
Read mode counts primary mapped alignment records, on either strand, including
MAPQ=0 and duplicate-flagged records; it excludes unmapped, secondary and
supplementary records. For PETRI-seq the empty R1 mate is unmapped, so the
mapped transcript R2 contributes one read count. An alignment is counted in
every gene whose interval overlaps its reference span, with no fractional
assignment: overlapping genes can share a read's count. This is unstranded
interval counting, not stranded or exon-union quantification.

Outputs are `genecount.N.h5`: AnnData-compatible CSR matrices with cells in
`obs/_index`, gene locus tags in `var/_index`, and read counts in `X`.
Each matrix retains **all** of that library's valid barcode combinations,
including rows of zeros for cells with no gene-overlapping reads. The
`gene_cells.N.txt` files preserve the row selection, and
`yanai_gene_count_libraries.csv` records the independent library, BAM,
matrix, annotation and counting rules. Libraries are never merged together.
The counting stage is reusable on existing alignments by sourcing
`yanai_petriseq/count_genes.R` and calling `RunYanaiGeneCounts()`.

The Yanai workflow sets --bwamem2-min-score 20 because its transcript reads
are 26 bases long; the stock score threshold 30 would discard even perfect
26-base matches. Other Bascet workflows retain the default threshold 30.
YANAI_BWAMEM2_MIN_SCORE can override the RNA setting.

The new Bascet PETRI-seq decoder keeps the original R2 transcript and
quality scores and the 7-base R1 UMI. Technical R1 sequence is removed.
See Bascet's README for barcode positions and correction rules.

## Validation

The complete `D1_lib1` library (`SRR22254412`) passed `prep_yanai_dev.R`
using the local Bascet build: 44,872,739 input pairs, 36,175,938 accepted
output records, and 359,666 barcode combinations. The independent full-stream
oracle matched every output record's content and multiplicity, lexical cell
sorting, and all histogram counts. The report and exact local binary checksum
are saved in `v9_rnaseq_yanai/dev_debarcode/D1_lib1/`.

The real SRR22254357 prefix was decoded with an independent nearest-barcode
oracle: 9,222 of 10,000 pairs matched, across 5,815 barcode combinations.
Sequences, qualities, UMIs and lexical cell sorting matched exactly.
Five Bascet decoder tests cover error correction, ambiguity, short reads,
invalid layouts and ordering.

A two-library Zorn end-to-end test retained all 9,222 reads and all 5,815
combinations per library after sharding, including 4,061 singleton
combinations in each library. The libraries remained distinct.
The same workflow also completed BWA-MEM2 alignment and BAM indexing.
After extending the counting stage, this two-library workflow also produced
two separate 5,815-cell × 11,800-gene matrices. Each contained 7,674 assigned
read counts and 7,373 nonzero entries, including retained zero-count cells.

The gene-count integration fixture checks repeated UMIs remain separate read
counts, first/last-base overlaps, a single-base gene, opposite-strand reads,
multiple overlapping genes, MAPQ=0, duplicate flags, rejection of secondary,
supplementary and unmapped records, and zero-count cells. It also verifies
88-library shard ordering and that the default UMI mode still deduplicates.
Rust tests cover GFF attribute selection, coordinate conversion and alignment
flag handling. The earlier fixture results above were established during
development. Current download, preparation, verification and analysis entry
points use Rust, R and shell.

    Rscript dataset_download_test.R
    cargo test --locked --offline --manifest-path yanai_petriseq/verifier/Cargo.toml \
      --target-dir /tmp/yanai-verifier-target

The downloader verifies original ENA sizes and MD5s, prepared mate receipts
and complete BGZF CRCs, including when reusing existing mates. Use repeated
`--run` options with `run_yanai_download.sh` to audit selected source libraries.
The wrapper takes the dataset download lock. Use it for concurrent invocations.

After gene counting and UMAP generation:

    Rscript yanai_petriseq/verify_gene_counts.R /path/to/dataset
    Rscript verify_umaps.R yanai

The gene-input audit checks every valid barcode, zero-count cell, reference
gene, sparse count and UMAP input depth, together with BAM integrity and exact
reference contig identities and lengths. It saves `gene_input_validation.tsv`.
The complete audit additionally verifies all 88 libraries and saved UMAP
coordinates, Seurat caches, plots and selection reports.

For a complete run with temporary storage on Husky:

    YANAI_ARRAY_WORKERS=4 bash run_dataset_all.sh yanai

This calls the existing full workflow with `TMPDIR` under the dataset's
`temp/R_workflow` directory, so a full local `/tmp` does not prevent R from
starting later stages. `DATASET_TMPDIR` overrides this storage location.
The same wrapper supports `ling` and `gorg`; it preserves the normal analysis
and completion-audit stages.

To check counting while the remaining full-run libraries are still aligning:

    Rscript validate_yanai_counting.R

This defaults to the complete, already aligned S. aureus `D3_lib1` library.
It calls the same `RunYanaiGeneCounts` function as `prep_yanai.R` in a fresh
`validation_counts_D3_lib1` directory, then independently scans the whole BAM
with Rsamtools and compares every cell/gene count using R interval overlaps.
The check covers primary mapped read counting, repeated UMIs, either strand,
overlapping genes and all seeded zero-count barcodes. It saves
`gene_input_validation.tsv` and `bam_gene_counts.validation.json`.
These validation outputs are kept separate from the main workflow.
Set `YANAI_COUNT_VALIDATION_LIBRARY` and a fresh `YANAI_COUNT_VALIDATION_ROOT`
to test another completed library or repeat the check.
