#!/usr/bin/env Rscript
# Run Bascet through Zorn's asynchronous LocalRunner and sample its process RSS.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1 || (length(args) == 1 && args != "--dry-run")) {
    stop("Usage: Rscript run_performance_shard_benchmark.R [--dry-run]")
  }
  root <- path.expand(Sys.getenv("PERFORMANCE_SHARD_ROOT", "/husky/henriksson/atrandi/v9_performance_shard"))
  cells <- as.integer(strsplit(Sys.getenv("BENCHMARK_CELLS", "1 10 100 1000 5000 10000 50000 276243"), "[[:space:]]+")[[1]])
  threads <- as.integer(Sys.getenv("BENCHMARK_THREADS", "40"))
  memory <- Sys.getenv("BENCHMARK_MEM", "50 GB")
  rss_interval <- as.numeric(Sys.getenv("BENCHMARK_RSS_INTERVAL", "0.1"))
  zorn_mode <- Sys.getenv("BENCHMARK_ZORN_MODE", "source")
  keep_outputs <- Sys.getenv("BENCHMARK_KEEP_OUTPUTS", "0") == "1"
  stopifnot(length(cells) > 0, !anyNA(cells), !anyDuplicated(cells),
            all(cells >= 1 & cells <= 331776), !is.na(threads), threads >= 6,
            is.finite(rss_interval), rss_interval > 0,
            zorn_mode %in% c("source", "installed"))
  read_pairs_for_cells <- function(n) ifelse(n == 276243L, 3250000000, n * 11765)
  plan <- data.frame(cells = cells, read_pairs = read_pairs_for_cells(cells),
                     repetitions = ifelse(cells == 276243L, 1L, ifelse(cells <= 5000, 3L, 2L)))
  if (identical(args, "--dry-run")) {
    cat("Threads:", threads, "; memory:", memory, "; shards: 10\n")
    write.table(plan, stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
    return(invisible(NULL))
  }

  required <- c("stringr", "fs", "naturalsort", "jsonlite", "ps")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  zorn_dir <- normalizePath(path.expand(Sys.getenv("ZORN_DIR", "~/github/zorn")), mustWork = TRUE)
  zorn_files <- file.path(zorn_dir, "R", c("job_general.R", "job_local.R", "jobscript.R",
                                          "bascet_file.R", "zorn.R", "shell.R", "container.R", "debarcoding.R"))
  if (zorn_mode == "installed") {
    library(Zorn)
    # Internal file-name helper, also used by Zorn's public wrappers.
    detectShardsForFile <- getFromNamespace("detectShardsForFile", "Zorn")
  } else {
    for (file in zorn_files) source(file, local = .GlobalEnv)
  }
  stopifnot(file.exists("/usr/bin/time"), dir.exists("/proc"))
  bascet_dir <- normalizePath(path.expand(Sys.getenv("BASCET_DIR", "~/github/bascet")), mustWork = TRUE)
  binary <- file.path(bascet_dir, "target", "release", "bascet")
  root <- normalizePath(root, mustWork = TRUE)
  sha256 <- function(paths) {
    result <- system2("sha256sum", shQuote(paths), stdout = TRUE)
    if (!is.null(attr(result, "status"))) stop("SHA256 calculation failed")
    setNames(substr(result, 1, 64), paths)
  }
  binary_hash <- unname(sha256(binary))
  version <- system2(binary, "--version", stdout = TRUE)
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  script_hash <- if (length(script) == 1) unname(sha256(script)) else NA_character_

  # Verify the generated inputs before starting any timed step.
  inputs <- lapply(cells, function(n) {
    dataset <- file.path(root, paste0("cells_", n))
    parameters <- read.delim(file.path(dataset, "parameters.tsv"))
    values <- setNames(parameters$value, parameters$parameter)
    stopifnot(as.integer(values[["cells"]]) == n,
              as.numeric(values[["read_pairs"]]) == read_pairs_for_cells(n),
              values[["compression"]] == "BGZF")
    expected <- readLines(file.path(dataset, "checksums.sha256"))
    names <- sub("^[0-9a-f]{64} [ *]", "", expected)
    stopifnot(setequal(names, c("reads_R1.fastq.gz", "reads_R2.fastq.gz")), length(names) == 2)
    message("Verifying inputs for ", n, " cells")
    actual <- sha256(file.path(dataset, names))
    if (!identical(unname(actual), substr(expected, 1, 64))) stop("Input checksum mismatch: ", dataset)
    list(directory = dataset, parameters = values, sha256 = actual)
  })
  names(inputs) <- as.character(cells)

  started <- format(Sys.time(), "%Y%m%dT%H%M%OS6Z", tz = "UTC")
  output <- path.expand(Sys.getenv("BENCHMARK_OUT", file.path(root, "benchmarks", paste0("run_", started))))
  if (file.exists(output)) stop("Results directory already exists: ", output)
  if (!dir.create(output, recursive = TRUE)) stop("Cannot create results directory: ", output)
  output <- normalizePath(output)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  runner <- LocalRunner(direct = FALSE, showScript = FALSE, ncpu = as.character(threads), mem = memory)
  metadata <- list(started_utc = started, status = "running", bascet = binary,
                   bascet_version = version, bascet_sha256 = binary_hash,
                   script_sha256 = script_hash,
                   zorn_dir = zorn_dir, zorn_source_sha256 = sha256(zorn_files),
                   zorn_mode = zorn_mode,
                   zorn_installed_path = if(zorn_mode == "installed") find.package("Zorn") else NULL,
                   zorn_installed_version = if(zorn_mode == "installed") as.character(packageVersion("Zorn")) else NULL,
                   threads = threads, memory = memory, shards = 10L, compression_level = 6L,
                   plan = plan, keep_outputs = keep_outputs, inputs = inputs,
                   host = as.list(Sys.info()), R_version = R.version.string,
                   cpuinfo = if (file.exists("/proc/cpuinfo")) readLines("/proc/cpuinfo") else NULL,
                   environment = as.list(Sys.getenv(c("OMP_NUM_THREADS", "RAYON_NUM_THREADS", "TMPDIR"))),
                   cache_policy = "No cache dropping or warmup; input hashes verified first. Sharding follows debarcoding immediately.",
                   rss_interval_seconds = rss_interval,
                   rss_measurement = "Sampled Bascet PID RSS from ps/procfs; R launcher and worker excluded. GNU time records peak RSS separately, including peaks between samples.",
                   run_note = Sys.getenv("BENCHMARK_NOTE", ""),
                   timing = "GNU time Bascet subprocess wall seconds, including startup and I/O; wrapper_wall_seconds separately records the full async Zorn call and wait.")
  save_metadata <- function() jsonlite::write_json(metadata, file.path(output, "metadata.json"),
                                                  auto_unbox = TRUE, pretty = TRUE)
  save_metadata()
  capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
  on.exit({
    if (metadata$status == "running") metadata$status <- "interrupted"
    metadata$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    save_metadata()
  }, add = TRUE)
  timing_file <- file.path(output, "timings.tsv")
  timings <- data.frame(cells = integer(), read_pairs = numeric(), `repeat` = integer(),
                        step = character(), wall_seconds = numeric(), wrapper_wall_seconds = numeric(),
                        peak_rss_bytes = numeric(), exit_code = integer(), check.names = FALSE)
  write.table(timings, timing_file, sep = "\t", row.names = FALSE, quote = FALSE)
  rss_file <- file.path(output, "rss.tsv")
  empty_rss <- data.frame(cells = integer(), `repeat` = integer(), step = character(),
                          elapsed_seconds = numeric(), pid = integer(), rss_bytes = numeric(), check.names = FALSE)
  write.table(empty_rss, rss_file, sep = "\t", row.names = FALSE, quote = FALSE)
  wrappers <- character()
  on.exit(unlink(wrappers), add = TRUE)
  profile_instance <- function(instance, step) {
    resource_file <- file.path(getwd(), paste0(step, ".resources.tsv"))
    wrapper <- tempfile(pattern = "zorn-benchmark-time-", tmpdir = "/tmp", fileext = ".sh")
    writeLines(c("#!/usr/bin/env bash", paste("exec /usr/bin/time -f", shQuote("%e\\t%M"),
                                                        "-o", shQuote(resource_file), '"$@"')), wrapper)
    Sys.chmod(wrapper, "0755")
    wrappers <<- c(wrappers, wrapper)
    file.copy(wrapper, paste0(step, ".time-wrapper.sh"))
    instance@prependCmd <- wrapper
    instance
  }
  timed_step <- function(n, repetition, step, expr) {
    message("cells=", n, " repeat=", repetition, " step=", step)
    error <- NULL
    job <- NULL
    samples <- list()
    on.exit(if(!is.null(job) && is(job, "LocalJob") && job@proc$is_alive()) CancelJob(job), add = TRUE)
    start <- proc.time()[["elapsed"]]
    tryCatch({
      job <- force(expr)
      if(!is(job, "LocalJob")) stop("Expected a new asynchronous LocalJob")
      worker <- job@proc$as_ps_handle()
      while(job@proc$is_alive()) {
        descendants <- tryCatch(ps::ps_children(worker, recursive = TRUE), error = function(e) list())
        for(handle in descendants) {
          value <- tryCatch({
            if(ps::ps_exe(handle) != binary) NULL else {
              data.frame(cells = n, `repeat` = repetition, step = step,
                         elapsed_seconds = proc.time()[["elapsed"]] - start,
                         pid = ps::ps_pid(handle), rss_bytes = unname(ps::ps_memory_info(handle)[["rss"]]),
                         check.names = FALSE)
            }
          }, error = function(e) NULL)
          if(!is.null(value)) samples[[length(samples) + 1L]] <- value
        }
        Sys.sleep(rss_interval)
      }
      WaitForJob(job)
    }, error = function(e) error <<- e)
    elapsed <- proc.time()[["elapsed"]] - start
    trace <- if(length(samples)) do.call(rbind, samples) else empty_rss
    write.table(trace, paste0(step, ".rss.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    if(nrow(trace)) write.table(trace, rss_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    resource_file <- paste0(step, ".resources.tsv")
    resources <- if(file.exists(resource_file)) readLines(resource_file, warn = FALSE) else character()
    resources <- resources[grepl("^[0-9.]+\\t[0-9]+$", resources)]
    metrics <- if(length(resources)) as.numeric(strsplit(tail(resources, 1), "\t")[[1]]) else c(NA_real_, NA_real_)
    if(is.null(error) && anyNA(metrics)) error <- simpleError("Missing GNU time metrics")
    row <- data.frame(cells = n, read_pairs = read_pairs_for_cells(n), `repeat` = repetition,
                       step = step, wall_seconds = metrics[1], wrapper_wall_seconds = elapsed,
                       peak_rss_bytes = metrics[2] * 1024, exit_code = as.integer(!is.null(error)), check.names = FALSE)
    timings <<- rbind(timings, row)
    write.table(row, timing_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    message("  ", metrics[1], " Bascet seconds; peak RSS ", round(metrics[2] / 1024, 1), " MiB; ", nrow(trace), " samples")
    if (!is.null(error)) stop(error)
  }
  summarize <- function() {
    successful <- timings[timings$exit_code == 0, , drop = FALSE]
    if (!nrow(successful)) return(invisible(NULL))
    summary <- do.call(rbind, lapply(split(successful, list(successful$cells, successful$step), drop = TRUE), function(d) {
      data.frame(cells = d$cells[1], read_pairs = d$read_pairs[1], step = d$step[1],
                 repetitions = nrow(d), median_seconds = median(d$wall_seconds),
                 min_seconds = min(d$wall_seconds), max_seconds = max(d$wall_seconds))
    }))
    write.table(summary, file.path(output, "summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  on.exit(summarize(), add = TRUE)
  message("Results: ", output)

  tryCatch({
    for (i in seq_len(nrow(plan))) {
      n <- plan$cells[i]
      dataset <- inputs[[as.character(n)]]$directory
      for (repetition in seq_len(plan$repetitions[i])) {
        repeat_dir <- file.path(output, paste0("cells_", n), paste0("rep_", repetition))
        work <- file.path(repeat_dir, "work")
        temp <- file.path(work, "temp")
        dir.create(temp, recursive = TRUE)
        # Keep Zorn logs outside work, so they survive intermediate cleanup.
        setwd(repeat_dir)
        instance <- getBascetDevDir(bascet_dir, tempdir = temp)
        rawmeta <- data.frame(dir = dataset, r1 = "reads_R1.fastq.gz", r2 = "reads_R2.fastq.gz", prefix = "sim")
        saveRDS(list(rawmeta = rawmeta, runner = runner, instance = instance,
                     threads = threads, totalMem = memory, compressionLevel = 6L,
                     maxShardSize = "1t", numOutputShards = 10L), "zorn_parameters.rds")
        debarcode_instance <- profile_instance(instance, "debarcode")
        timed_step(n, repetition, "debarcode", BascetDebarcode(
          work, rawmeta, maxShardSize = "1t", outputName = "debarcoded", chemistry = "atrandi-wgs",
          numThreads = threads, totalMem = memory, compressionLevel = 6L,
          runner = runner, bascetInstance = debarcode_instance
        ))
        parsed <- detectShardsForFile(work, "debarcoded")
        if (length(parsed) != 1) stop("Expected one debarcoded output")
        histogram <- read.delim(paste0(file.path(work, parsed), ".hist"), header = FALSE,
                                 col.names = c("cell", "pairs"))
        if (!nrow(histogram)) stop("No cells retained after debarcoding")
        # Construct the sharding plan from every recovered cell. PrepareSharding's
        # knee filter uses a strict count cutoff, which can exclude the sole cell.
        debstat <- structure(list(bascetRoot = work, numgroup = 1L,
                                  meta = data.frame(group = 1L), inputFiles = parsed,
                                  list_picked_cells = list(as.character(histogram$cell))), class = "debstat")
        saveRDS(debstat, "sharding_plan.rds")
        shardify_instance <- profile_instance(instance, "shardify")
        timed_step(n, repetition, "shardify", BascetShardify(
          debstat, numOutputShards = 10L, outputName = "sharded", numThreads = threads,
          totalMem = memory, runner = runner, bascetInstance = shardify_instance
        ))
        if (length(detectShardsForFile(work, "sharded")) != 10) stop("Expected exactly 10 shard outputs")
        summarize()
        if (!keep_outputs) unlink(work, recursive = TRUE)
      }
    }
    if (!identical(unname(sha256(binary)), binary_hash) || !identical(sha256(zorn_files), metadata$zorn_source_sha256)) {
      stop("Bascet binary or Zorn source changed during the benchmark")
    }
    metadata$status <- "complete"
  }, error = function(e) {
    metadata$status <<- "failed"
    metadata$error <<- conditionMessage(e)
    stop(e)
  })
}

main()
