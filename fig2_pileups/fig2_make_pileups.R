library(duckdb)
library(DBI)

# ── Constants ─────────────────────────────────────────────────────────────────

# BIN_SIZE <- 1000L
BIN_SIZE <- 2500L

CSV_OPTS <- paste(
  "delim='\\t'",
  "header=true",
  "max_line_size=100000000",
  "strict_mode=false",
  "quote=''",
  "auto_detect=true",
  sep = ", "
)

# ── DuckDB setup ──────────────────────────────────────────────────────────────

con <- dbConnect(duckdb())

dbExecute(con, sprintf("
  CREATE VIEW qc_raw AS
  SELECT *,
    CASE
      WHEN reference = 'phix' THEN 'phix'
      WHEN reference LIKE 'NZ_%%' THEN 'NZ_CDNB'
      ELSE regexp_replace(reference, '_contig_\\d+$', '')
    END AS organism_id
  FROM (
    SELECT *, 'MDA' AS \"group\" FROM read_csv('./data/mda/qc_coverage.tsv', %s)
    UNION ALL
    SELECT *, 'PTA' AS \"group\" FROM read_csv('./data/pta/qc_coverage.tsv', %s)
  )
", CSV_OPTS, CSV_OPTS))

# ── Contig offsets ────────────────────────────────────────────────────────────

dbExecute(con, "
  CREATE TABLE contig_offsets AS
  WITH distinct_contigs AS (
    SELECT DISTINCT organism_id, reference, sizeof_reference
    FROM qc_raw
  ),
  ordered AS (
    SELECT *,
      ROW_NUMBER() OVER (PARTITION BY organism_id ORDER BY reference) AS contig_rank
    FROM distinct_contigs
  )
  SELECT *,
    COALESCE(
      SUM(sizeof_reference) OVER (
        PARTITION BY organism_id
        ORDER BY contig_rank
        ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
      ),
      0
    ) AS contig_offset
  FROM ordered
")

# ── RLE decode -> genome bins ─────────────────────────────────────────────────

dbExecute(con, sprintf("
  CREATE TABLE pileup_binned AS
  WITH genome_sizes AS (
    SELECT organism_id, SUM(sizeof_reference) AS genome_size
    FROM contig_offsets
    GROUP BY organism_id
  ),
  contig_bins AS (
    SELECT c.organism_id, c.reference, c.sizeof_reference, c.contig_offset,
      g.genome_size,
      GREATEST(CEIL(c.sizeof_reference::DOUBLE / %d), 1)::BIGINT AS n_contig_bins
    FROM contig_offsets c
    JOIN genome_sizes g USING (organism_id)
  ),
  runs AS (
    SELECT
      q.\"group\", q.id, q.organism_id, q.reference,
      cb.contig_offset, cb.sizeof_reference,
      cb.genome_size, cb.n_contig_bins,
      UNNEST(string_split(q.pileup, ',')) AS rle_token
    FROM qc_raw q
    JOIN contig_bins cb USING (organism_id, reference)
    WHERE q.pileup IS NOT NULL
      AND q.pileup != ''
  ),
  parsed AS (
    SELECT
      \"group\", id, organism_id, reference, genome_size,
      contig_offset, sizeof_reference, n_contig_bins,
      CAST(string_split(rle_token, ':')[1] AS BIGINT) AS local_start,
      CAST(string_split(rle_token, ':')[2] AS DOUBLE)  AS run_depth,
      CAST(string_split(rle_token, ':')[3] AS BIGINT)  AS run_length
    FROM runs
    WHERE rle_token != ''
  ),
  with_ends AS (
    SELECT *,
      local_start + run_length AS local_end,
      LEAST(CAST(FLOOR(local_start::DOUBLE / sizeof_reference * n_contig_bins) AS BIGINT), n_contig_bins - 1) AS first_bin,
      LEAST(CAST(FLOOR((local_start + run_length - 1)::DOUBLE / sizeof_reference * n_contig_bins) AS BIGINT), n_contig_bins - 1) AS last_bin
    FROM parsed
    WHERE run_length > 0
  ),
  expanded AS (
    SELECT
      \"group\", id, organism_id, reference,
      genome_size, contig_offset, sizeof_reference, n_contig_bins, run_depth,
      local_start, local_end,
      UNNEST(generate_series(first_bin, last_bin)) AS local_bin
    FROM with_ends
  )
  SELECT
    \"group\", id, organism_id, reference,
    contig_offset + (local_bin + 0.5) * sizeof_reference / n_contig_bins AS position,
    sizeof_reference::DOUBLE / n_contig_bins AS bin_width,
    SUM(
      run_depth * (
        LEAST(local_end, (local_bin + 1)::DOUBLE * sizeof_reference / n_contig_bins)
        - GREATEST(local_start, local_bin::DOUBLE * sizeof_reference / n_contig_bins)
      ) / (sizeof_reference::DOUBLE / n_contig_bins)
    ) AS weighted_depth
  FROM expanded
  GROUP BY \"group\", id, organism_id, reference, contig_offset, sizeof_reference, n_contig_bins, local_bin, genome_size
", BIN_SIZE))

# ── Read into R ───────────────────────────────────────────────────────────────

pileup.raw <- dbGetQuery(con, "SELECT * FROM pileup_binned")
contigs    <- dbGetQuery(con, "SELECT * FROM contig_offsets")
genome.sizes <- dbGetQuery(con, "
  SELECT organism_id, SUM(sizeof_reference) AS genome_size
  FROM contig_offsets
  GROUP BY organism_id
")

dbDisconnect(con)
