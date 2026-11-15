# Ling / Aarestrup single-cell metagenomes

Download all four datasets from Ling et al. (2025), DOI
[10.3389/fmicb.2024.1516656](https://doi.org/10.3389/fmicb.2024.1516656):

    ./run_ling_download.sh

Default destination: /husky/henriksson/atrandi/v9_ling/
The script is resumable, preserves existing files, and checks the exact
Zenodo byte count and MD5 for every file. Six downloads run concurrently.
Override concurrency with --workers N. Select one or more datasets with
--library feces_deep, --library feces_shallow, --library sewage_deep,
or --library sewage_shallow.

| Dataset | Zenodo record | Files | Compressed GB |
|---|---|---:|---:|
| Sewage, deep | [13323680](https://zenodo.org/records/13323680) | 27 | 36.396 |
| Feces, deep | [13348502](https://zenodo.org/records/13348502) | 21 | 45.240 |
| Feces, shallow | [13353442](https://zenodo.org/records/13353442) | 21 | 28.454 |
| Sewage, shallow | [13348988](https://zenodo.org/records/13348988) | 27 | 44.974 |

The manifest pins all **96 files, 155.063 GB** with original filenames,
URLs, sizes, checksums, library labels, barcode-A values and file roles.
The four original Zenodo JSON metadata records are retained here and copied
to the destination. Files are stored under source/DATASET/ with per-file
verification receipts.

## The deposited data are split/trimmed intermediates

These are not untouched sequencing FASTQs and are not per-cell files.
Each of the 32 barcode-A groups has three matched files:

- R1.trim.fastq.gz: biological forward read and its measured qualities.
- R2_bc.fastq.gz: three concatenated 8-base cell-barcode components
  (24 bases total), with measured barcode qualities.
- R2_frag.trim.fastq.gz: biological reverse read and its measured qualities.

Barcode A is encoded by the eight-base suffix in the original filename.
The other three components vary between reads inside each group. The
paper used these three components without a whitelist. Retain barcode-A
and independent dataset/library identities when forming cell IDs.

The first 12,000 sewage-deep ATGTCTGC triplets were checked: all three files
match by read identity; R2_bc is exactly 24 bases. The rbc field in the
headers retains a 36-base original barcode/spacer sequence; the 24-base
barcode read is its positions [0,8), [12,20), and [24,32).
R1 and R2 read numbers are concatenated to the final Illumina coordinate
before ":N", so literal complete-header equality is inappropriate.

Current Bascet atrandi-wgs expects a four-barcode cassette embedded in R2
and is not a drop-in decoder for these three-file intermediates. Merely
renaming R2_frag to R2 would lose cell identity. Reconstructing a nominal
raw R2 would require inventing qualities for removed technical sequence,
so the downloader deliberately preserves the deposited source files.
The importer below uses barcode A plus the observed 24-base barcode as
cell ID and keeps both biological mates with their original qualities.

## Import to Bascet

Once all 96 source files have verification receipts:

    ./run_ling_import.sh --root /husky/henriksson/atrandi/v9_ling

The launcher builds and runs the local `ling2tirp` Rust crate. Its target
directory defaults to `/tmp/ling2tirp-target` (override
`LING_CARGO_TARGET_DIR`); no Python converter is invoked by this workflow.

The command refuses to start if any source is missing or its verification
receipt disagrees with the pinned size/MD5. It does not invoke the current
atrandi-wgs decoder: these sources already have separate barcodes.

Each barcode-A group produces one debarcoded.N.tirp.gz, its .hist cell-count
histogram, and a JSON receipt. The 32 outputs are sorted by cell ID using
GNU stable sort with LC_ALL=C. Cell IDs have the form:

    sewage_deep_ATGTCTGC_CAAGCAGTAGACCGATTGTGGTCT
    library      A        observed 24-base barcode

All barcode strings are retained exactly as observed, including N calls.
No barcode whitelist, correction, abundance threshold, host filtering,
UMI, or synthetic quality scores are added. R1/Q1 and R2/Q2 contain the
two biological mates and their measured qualities; the final UMI field
is empty.

The output compression is Bascet BBGZ with cell IDs in every block header:
ID extra field, BC extra field last, FLG=5, at most 60 KiB uncompressed
payload, CRC/ISIZE trailers and a standard EOF block. Raw DEFLATE uses
Z_SYNC_FLUSH followed by final 03 00, allowing Bascet's compressed-block
shard merger to operate without decompression. Complete records stay
inside blocks, and blocks never mix cells. No tabix index is needed.

The debarcoded.meta CSV records source triplets, barcode A, shard indices,
record/cell counts, and the four biological library prefixes. It is
compatible with Zorn PrepareSharding/BascetShardify; grouping by prefix
merges the A groups within a library while full cell IDs prevent collisions.

External sorting uses a 4 GiB memory buffer by default. Adjust it with
--sort-memory (GNU sort syntax, for example 8G); choose fast temporary
storage with --sort-temp. The default is ROOT/temp/ling_import_rust. Eight
groups run concurrently by default; `--workers N` changes concurrency, with
the memory buffer allocated separately to each worker. Sorting spills plain, uncompressed TIRP
text: temporary I/O and disk use can be much greater than compressed input
size. The complete unsorted TSV is not separately retained.

Validation used 12,000 actual sewage-deep triplets and a 1 MiB sort buffer
to force external-sort spills. All 12,000 biological sequence/quality pairs
and their exact observed barcode identities survived; 497 observed cell
IDs, histogram counts, ordering and per-block CRCs were checked. Bascet
shardify then preserved the exact record multiset.

## Prepare count/sketch inputs for UMAP analysis

    Rscript /husky/henriksson/atrandi/revision_scripts/prep_ling.R

The prep script imports the source triplets, selects cells, groups barcode-A
files by their biological library, and prepares three representations:

| Method | Output | Downstream Zorn reader |
|---|---|---|
| CountSketch of 31-mers | `countsketch_mat.feather` | `BascetLoadCountSketchMatrix()` |
| Kraken taxon read counts | `kraken_counts.N.h5` | `ReadBascetCountMatrix(root, "kraken_counts")` |
| Informative 31-mer counts | `kmer_counts.N.h5` | `ReadBascetCountMatrix(root, "kmer_counts")` |

Cell selection follows `prep_mock.R`'s 95th read-count percentile by default.
Set `LING_KEEP_ALL=1` to retain every imported barcode combination, or
`LING_MIN_QUANTILE` to change the percentile. `LING_RUN_KRAKEN=0` skips
the Kraken representation. `LING_ROOT`, `LING_RUNNER`, `LING_THREADS`,
`LING_MEM`, `LING_BASCET_DIR`, and `LING_KRAKEN_DB` configure execution.

The prep script ends with count/sketch inputs. It does not invoke FASTQC,
SKESA, GECCO, QUAST, reference alignment, or host filtering.

## Run through all UMAPs

    ./run_ling_all.sh

This imports with eight workers, runs `prep_ling.R` with a local Bascet
checkout, then runs `post_ling.R`. The postprocessing uses the normal
R workflow in `post_unified.R`: Zorn's readers, Seurat/Signac
`RunTFIDF`, `FindTopFeatures(min.cutoff='q0')`, `RunSVD`, and
`RunUMAP(reduction='lsi', dims=1:30)` for Kraken and informative 31-mer
counts; `BascetLoadCountSketchMatrix` and `CountSketchUMAP` for sketches.
There is no separate Python UMAP implementation.

The normal depth thresholds are retained: Kraken counts >1,000 and
CountSketch depth >40,000, configurable with `LING_MIN_KRAKEN` and
`LING_MIN_SKETCH_KMERS`. All cell selections are recorded. Taxonomy comes
from `KrakenFindConsensusTaxonomy`; library, feces/sewage and deep/shallow
labels are added without batch correction. Mock-only alignment/host plots
and taxon exclusions are not applied.

Seurat caches are `cache_adata_kraken.RDS`, `cache_adata_infokmer.RDS`,
and `cache_adata_cs.RDS` in the dataset directory. Coordinates, cell
selections, run parameters and R session information go in `ROOT/umap/`.
PDF/PNG plots go in `plots/v9_ling/` in the scripts repository, or
`LING_PLOT_DIR`. `all_umaps_library.pdf` shows all three reductions.
For postprocessing alone, run `Rscript post_ling.R`.

The R integration fixture is `Rscript ling_data/test_umaps.R`; it exercises
all three standard reductions on four small synthetic libraries and checks
the saved objects and finite coordinates. The Rust converter also reproduced
all 12,000 previously verified real triplets byte-for-byte after decompression;
Bascet shard merging preserved those same records and their histogram counts.
