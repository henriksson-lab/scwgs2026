# Yanai PETRI-seq reference for BWA-MEM2

The destination requested for this Yanai reference is
/husky/henriksson/atrandi/bwamem_ref/ling/.

    ./run_yanai_reference_download.sh
    ./run_yanai_reference_index.sh

The download script retrieves pinned NCBI FASTA, GFF3 and assembly reports,
checks exact byte counts and published MD5s, validates FASTA/GFF contig
agreement and gene-ID uniqueness, and combines the assemblies into all.fa
and all.gff. The index script uses Bascet's own BWA-MEM2 wrapper and records
the index artifact SHA256s. Both scripts use R and shell, preserve existing
matching files and refuse to replace different references or index artifacts.

## Assemblies

These are the reference genomes named in the methods of Pountain et al.,
[Nature 626, 661–669 (2024)](https://doi.org/10.1038/s41586-023-06974-w).

| Species | Paper reference | Downloaded RefSeq reference |
|---|---|---|
| E. coli K-12 MG1655 | GCA_000005845.2 | [GCF_000005845.2, ASM584v2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/) |
| S. aureus USA300_FPR3757 | GCF_000013465.1 | [GCF_000013465.1, ASM1346v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013465.1/) |
| B. subtilis 168 | GCF_000009045.1 | [GCF_000009045.1, ASM904v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009045.1/) |

NCBI's assembly report identifies the E. coli GenBank and RefSeq assemblies
as identical and maps U00096.3 to NC_000913.3 with "=". This means identical
underlying genome sequence, not identical FASTA headers or compressed
file bytes. The RefSeq FASTA and matching RefSeq GFF are used together.

The combined reference retains the six unique RefSeq sequence accessions:

- NC_000913.3: E. coli chromosome, 4,641,652 bases.
- NC_007793.1: S. aureus chromosome, 2,872,769 bases.
- NC_007790.1, NC_007791.1, NC_007792.1: S. aureus plasmids, 3,125,
  4,439 and 37,136 bases.
- NC_000964.3: B. subtilis chromosome, 4,215,606 bases.

Total: 11,774,727 bases. Matching assembly reports, source files, checksums,
and combined-reference SHA256s are retained beside the index.

## Alignment and counting

Use /husky/henriksson/atrandi/bwamem_ref/ling/all.fa as the BWA-MEM2 index
prefix and /husky/henriksson/atrandi/bwamem_ref/ling/all.gff for counting.
The GFF contains 11,800 gene records with globally unique locus_tag and ID
attributes, and Name attributes for all gene records. Zorn settings:

    useFeature = "gene"
    attrGeneId = "locus_tag"
    attrGeneName = "Name"

The NCBI GFF annotations are pinned by checksum for the downloaded assembly
versions; these do not claim to recreate the exact annotation snapshot
used in the original paper's processed matrices.

## Short-read compatibility validation

The index was loaded by both native bwa-mem2 and Bascet's wrapper. Perfect
reference-derived 100-base controls map for all three species. The same
26-base controls are unmapped at the default minimum alignment score 30,
but map as 26M with -T 20, for both native bwa-mem2 and Bascet's wrapper.

Bascet's align workflow was also tested on a canonical BGZF TIRP with empty
R1/Q1 and the genomic read in R2/Q2. It completes without a panic and retains
CB/UB tags; its paired driver additionally emits empty first-mate records.
Using score 20 for the Yanai RNA workflow is necessary to retain these
short transcripts. Validation files and results are recorded at the
reference destination under validation/.

The completed driver test with --bwamem2-min-score 20 mapped all three
perfect 26-base R2 reads as 26M and preserved CB/UB tags. Empty R1/Q1
records did not cause a panic. The minimum-score option retains its
standard default of 30 for other workflows.
