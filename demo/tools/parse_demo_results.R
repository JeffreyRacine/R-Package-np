#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else "."
out_dir <- if (length(args) >= 2L && nzchar(args[[2L]])) args[[2L]] else "."

root <- normalizePath(root, mustWork = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

parse_key_values <- function(line) {
  line <- sub("^.*DEMO_RESULT[[:space:]]+", "", line)
  pieces <- strsplit(line, "[[:space:]]+", perl = TRUE)[[1L]]
  out <- list()
  for (piece in pieces) {
    pos <- regexpr("=", piece, fixed = TRUE)[[1L]]
    if (pos <= 1L) next
    key <- substr(piece, 1L, pos - 1L)
    value <- substr(piece, pos + 1L, nchar(piece))
    out[[key]] <- value
  }
  out
}

path_info <- function(path) {
  rel <- substring(normalizePath(path, mustWork = TRUE), nchar(root) + 2L)
  parts <- strsplit(rel, .Platform$file.sep, fixed = TRUE)[[1L]]
  info <- list(engine = NA_character_, mode = NA_character_,
               ranks = NA_integer_, slaves = NA_integer_)

  if ("serial" %in% parts) {
    info$engine <- "serial"
    info$mode <- "serial"
    return(info)
  }

  session_idx <- match("session", parts)
  if (!is.na(session_idx)) {
    info$engine <- "session"
    info$mode <- "session"
    slave_part <- parts[grepl("^slaves_[0-9]+$", parts)]
    if (length(slave_part)) info$slaves <- as.integer(sub("^slaves_", "", slave_part[[1L]]))
    return(info)
  }

  rank_part <- parts[grepl("^ranks_[0-9]+$", parts)]
  if (length(rank_part)) {
    info$engine <- "mpi_launch"
    info$ranks <- as.integer(sub("^ranks_", "", rank_part[[1L]]))
    info$slaves <- info$ranks - 1L
    if ("attach" %in% parts) info$mode <- "attach"
    if ("profile" %in% parts) info$mode <- "profile"
    return(info)
  }

  legacy <- parts[grepl("^(.*_)?n_[0-9]+_(attach|profile)$", parts)]
  if (length(legacy)) {
    info$engine <- "mpi_launch_legacy"
    info$ranks <- as.integer(sub("^(.*_)?n_([0-9]+)_(attach|profile)$", "\\2", legacy[[1L]], perl = TRUE))
    info$slaves <- info$ranks - 1L
    info$mode <- sub("^(.*_)?n_[0-9]+_", "", legacy[[1L]], perl = TRUE)
    return(info)
  }

  info
}

num_or_na <- function(x) {
  if (is.null(x) || !length(x) || !nzchar(x[[1L]]) || identical(x[[1L]], "NA")) return(NA_real_)
  suppressWarnings(as.numeric(x[[1L]]))
}

int_or_na <- function(x) {
  out <- num_or_na(x)
  if (is.na(out)) NA_integer_ else as.integer(out)
}

first_match <- function(pattern, text, group = 1L) {
  m <- regexec(pattern, text, perl = TRUE)
  r <- regmatches(text, m)[[1L]]
  if (length(r) <= group) NA_character_ else r[[group + 1L]]
}

fallback_record <- function(path, text) {
  elapsed <- first_match("Elapsed time[[:space:]]*=[[:space:]]*([0-9.]+)", text)
  if (is.na(elapsed)) return(NULL)

  file <- basename(path)
  demo <- sub("_(serial|npRmpi_attach|npRmpi_profile|npRmpi_session)\\.Rout$", "", file, perl = TRUE)
  mode <- if (grepl("_serial\\.Rout$", file)) "serial" else
    if (grepl("_npRmpi_attach\\.Rout$", file)) "attach" else
      if (grepl("_npRmpi_profile\\.Rout$", file)) "profile" else
        if (grepl("_npRmpi_session\\.Rout$", file)) "session" else NA_character_

  n <- first_match("\\(([0-9]+)[[:space:]]+observations", text)
  if (is.na(n)) n <- first_match("Data:[[:space:]]+([0-9]+)[[:space:]]+training", text)
  if (is.na(n)) n <- first_match("Data\\):[[:space:]]+([0-9]+)[[:space:]]+training", text)
  if (is.na(n)) n <- first_match("n[[:space:]]*<-[[:space:]]*([0-9]+)", text)

  info <- path_info(path)
  if (is.na(info$mode) && !is.na(mode)) info$mode <- mode

  data.frame(
    run_id = basename(root),
    demo = demo,
    mode = info$mode,
    engine = info$engine,
    ranks = info$ranks,
    slaves = info$slaves,
    n = int_or_na(n),
    default_n = NA_integer_,
    elapsed = num_or_na(elapsed),
    bwmethod = NA_character_,
    method = NA_character_,
    family = NA_character_,
    case = NA_character_,
    tier = NA_character_,
    regtype = NA_character_,
    nomad = NA_character_,
    degree = NA_character_,
    degree.max = NA_character_,
    selected.degree = NA_real_,
    bwtype = NA_character_,
    objective = NA_real_,
    bandwidth = NA_real_,
    bandwidth.second = NA_real_,
    beta.second = NA_real_,
    coef.first = NA_real_,
    coef.second = NA_real_,
    density.first = NA_real_,
    distribution.first = NA_real_,
    statistic.first = NA_real_,
    statistic.second = NA_real_,
    p.value = NA_real_,
    estimate.first = NA_real_,
    estimate.second = NA_real_,
    num.feval = NA_real_,
    num.feval.fast = NA_real_,
    nomad.time = NA_real_,
    powell.time = NA_real_,
    fitted.length = NA_integer_,
    source_file = normalizePath(path, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

records <- list()
files <- list.files(root, pattern = "\\.Rout$", recursive = TRUE, full.names = TRUE)

for (path in files) {
  text <- read_text(path)
  lines <- grep("DEMO_RESULT", strsplit(text, "\n", fixed = TRUE)[[1L]], value = TRUE)
  info <- path_info(path)

  if (length(lines)) {
    for (line in lines) {
      kv <- parse_key_values(line)
      mode <- if (!is.null(kv$mode)) kv$mode else info$mode
      engine <- info$engine
      if (is.na(engine) && identical(mode, "session")) engine <- "session_legacy"
      records[[length(records) + 1L]] <- data.frame(
        run_id = basename(root),
        demo = if (!is.null(kv$demo)) kv$demo else sub("\\.Rout$", "", basename(path)),
        mode = mode,
        engine = engine,
        ranks = if (!is.na(info$ranks)) info$ranks else int_or_na(kv$ranks),
        slaves = if (!is.na(info$slaves)) info$slaves else int_or_na(kv$slaves),
        n = int_or_na(kv$n),
        default_n = int_or_na(kv$default_n),
        elapsed = num_or_na(kv$elapsed),
        bwmethod = if (!is.null(kv$bwmethod)) kv$bwmethod else NA_character_,
        method = if (!is.null(kv$method)) kv$method else NA_character_,
        family = if (!is.null(kv$family)) kv$family else NA_character_,
        case = if (!is.null(kv$case)) kv$case else NA_character_,
        tier = if (!is.null(kv$tier)) kv$tier else NA_character_,
        regtype = if (!is.null(kv$regtype)) kv$regtype else NA_character_,
        nomad = if (!is.null(kv$nomad)) kv$nomad else NA_character_,
        degree = if (!is.null(kv$degree)) kv$degree else NA_character_,
        degree.max = if (!is.null(kv$degree.max)) kv$degree.max else NA_character_,
        selected.degree = num_or_na(kv$selected.degree),
        bwtype = if (!is.null(kv$bwtype)) kv$bwtype else NA_character_,
        objective = num_or_na(kv$objective),
        bandwidth = num_or_na(kv$bandwidth),
        bandwidth.second = num_or_na(kv$bandwidth.second),
        beta.second = num_or_na(kv$beta.second),
        coef.first = num_or_na(kv$coef.first),
        coef.second = num_or_na(kv$coef.second),
        density.first = num_or_na(kv$density.first),
        distribution.first = num_or_na(kv$distribution.first),
        statistic.first = num_or_na(kv$statistic.first),
        statistic.second = num_or_na(kv$statistic.second),
        p.value = num_or_na(kv$p.value),
        estimate.first = num_or_na(kv$estimate.first),
        estimate.second = num_or_na(kv$estimate.second),
        num.feval = num_or_na(kv$num.feval),
        num.feval.fast = num_or_na(kv$num.feval.fast),
        nomad.time = num_or_na(kv$nomad.time),
        powell.time = num_or_na(kv$powell.time),
        fitted.length = int_or_na(kv$fitted.length),
        source_file = normalizePath(path, mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    }
  } else {
    rec <- fallback_record(path, text)
    if (!is.null(rec)) records[[length(records) + 1L]] <- rec
  }
}

results <- if (length(records)) do.call(rbind, records) else {
  data.frame(run_id = character(), demo = character(), mode = character(),
             engine = character(), ranks = integer(), slaves = integer(),
             n = integer(), default_n = integer(), elapsed = numeric(),
             bwmethod = character(), method = character(),
             family = character(), case = character(),
             tier = character(), regtype = character(), nomad = character(),
             degree = character(), degree.max = character(),
             selected.degree = numeric(),
             bwtype = character(), objective = numeric(), bandwidth = numeric(),
             bandwidth.second = numeric(), beta.second = numeric(),
             coef.first = numeric(), coef.second = numeric(),
             density.first = numeric(), distribution.first = numeric(),
             statistic.first = numeric(), statistic.second = numeric(),
             p.value = numeric(),
             estimate.first = numeric(), estimate.second = numeric(),
             num.feval = numeric(),
             num.feval.fast = numeric(), nomad.time = numeric(),
             powell.time = numeric(), fitted.length = integer(),
             source_file = character(),
             stringsAsFactors = FALSE)
}

for (name in c("method", "family", "case", "tier", "regtype", "nomad", "degree",
               "degree.max", "selected.degree", "bwtype", "objective", "num.feval",
               "bandwidth", "bandwidth.second", "beta.second", "coef.first", "coef.second",
               "density.first", "distribution.first",
               "statistic.first", "statistic.second", "p.value",
               "estimate.first", "estimate.second",
               "num.feval.fast", "nomad.time", "powell.time",
               "fitted.length")) {
  if (!name %in% names(results)) results[[name]] <- NA
}

results <- results[order(results$family, results$demo, results$case,
                         results$mode, results$ranks, results$slaves,
                         results$source_file), ]

write.csv(results, file.path(out_dir, "demo_results.csv"), row.names = FALSE)

compute_label <- function(row) {
  mode <- row[["mode"]]
  if (is.na(mode) || !nzchar(mode)) mode <- "unknown"
  if (identical(mode, "serial")) return("serial")
  if (identical(mode, "session")) {
    slaves <- row[["slaves"]]
    return(sprintf("session_s%02d", ifelse(is.na(slaves), 0L, slaves)))
  }
  ranks <- row[["ranks"]]
  if (is.na(ranks)) ranks <- row[["slaves"]] + 1L
  if (is.na(ranks)) return(mode)
  sprintf("%s_r%02d", mode, ranks)
}

write_wide <- function(results, path) {
  if (!nrow(results)) {
    write.csv(results, path, row.names = FALSE)
    return(invisible(TRUE))
  }

  wide <- results
  wide$compute <- vapply(seq_len(nrow(wide)), function(i) compute_label(wide[i, ]),
                         character(1L))
  id_cols <- c("run_id", "family", "case", "tier", "regtype", "bwmethod",
               "nomad", "degree", "degree.max", "selected.degree", "bwtype",
               "n", "default_n")
  id_cols <- id_cols[id_cols %in% names(wide)]
  wide <- wide[, c(id_cols, "compute", "elapsed"), drop = FALSE]
  names(wide)[names(wide) == "elapsed"] <- "seconds"
  out <- reshape(wide, idvar = id_cols, timevar = "compute",
                 direction = "wide")
  names(out) <- sub("^seconds\\.", "seconds_", names(out))
  write.csv(out, path, row.names = FALSE)
  invisible(TRUE)
}

write_wide(results, file.path(out_dir, "demo_results_wide.csv"))

line_for <- function(row) {
  compute <- if (!is.na(row[["slaves"]])) paste0("slaves=", row[["slaves"]]) else
    if (!is.na(row[["ranks"]])) paste0("ranks=", row[["ranks"]]) else
      if (identical(row[["mode"]], "serial")) "serial" else "compute=NA"
  sprintf("%s mode=%s %s n=%s elapsed=%.3f",
          row[["demo"]], row[["mode"]], compute,
          ifelse(is.na(row[["n"]]), "NA", as.character(row[["n"]])),
          row[["elapsed"]])
}

dat <- if (nrow(results)) vapply(seq_len(nrow(results)), function(i) line_for(results[i, ]), character(1L)) else character()
writeLines(dat, file.path(out_dir, "timing_all.dat"))

qmd <- c(
  "---",
  "title: \"npRmpi Demo Timing Summary\"",
  "format:",
  "  html:",
  "    toc: true",
  "    code-fold: true",
  "execute:",
  "  echo: false",
  "  warning: false",
  "  message: false",
  "---",
  "",
  "```{r}",
  "results <- read.csv('demo_results.csv', stringsAsFactors = FALSE)",
  "wide <- read.csv('demo_results_wide.csv', stringsAsFactors = FALSE)",
  "metadata_path <- file.path('..', 'RUN_METADATA.txt')",
  "metadata <- if (file.exists(metadata_path)) readLines(metadata_path, warn = FALSE) else character()",
  "num <- function(x) suppressWarnings(as.numeric(x))",
  "results$elapsed <- num(results$elapsed)",
  "results$slaves <- num(results$slaves)",
  "results$ranks <- num(results$ranks)",
  "results$n <- num(results$n)",
  "results$family[!nzchar(results$family)] <- results$demo[!nzchar(results$family)]",
  "```",
  "",
  "## Run Metadata",
  "",
  "```{r}",
  "if (length(metadata)) {",
  "  cat(paste(metadata, collapse = '\\n'))",
  "} else {",
  "  cat('No RUN_METADATA.txt found.')",
  "}",
  "```",
  "",
  "## Coverage",
  "",
  "```{r}",
  "coverage <- aggregate(elapsed ~ family + mode, results, length)",
  "names(coverage)[names(coverage) == 'elapsed'] <- 'rows'",
  "coverage <- coverage[order(coverage$family, coverage$mode), ]",
  "knitr::kable(coverage, row.names = FALSE)",
  "```",
  "",
  "## Elapsed Time",
  "",
  "```{r}",
  "elapsed <- results[, c('family', 'case', 'mode', 'ranks', 'slaves', 'n', 'elapsed')]",
  "elapsed <- elapsed[order(elapsed$family, elapsed$case, elapsed$mode, elapsed$ranks, elapsed$slaves), ]",
  "knitr::kable(elapsed, row.names = FALSE, digits = 3)",
  "```",
  "",
  "## Wide Elapsed Table",
  "",
  "```{r}",
  "knitr::kable(wide, row.names = FALSE, digits = 3)",
  "```",
  "",
  "## Option Drift Sentinels",
  "",
  "```{r}",
  "sentinel_cols <- c('family', 'case', 'mode', 'regtype', 'bwmethod', 'nomad', 'degree',",
  "                   'degree.max', 'selected.degree', 'bwtype', 'objective', 'bandwidth',",
  "                   'bandwidth.second', 'beta.second', 'coef.first', 'coef.second',",
  "                   'density.first', 'distribution.first', 'statistic.first',",
  "                   'statistic.second', 'p.value', 'estimate.first', 'estimate.second',",
  "                   'num.feval', 'num.feval.fast', 'fitted.length')",
  "sentinel_cols <- sentinel_cols[sentinel_cols %in% names(results)]",
  "sentinels <- results[, sentinel_cols, drop = FALSE]",
  "sentinels <- sentinels[order(sentinels$family, sentinels$case, sentinels$mode), ]",
  "knitr::kable(sentinels, row.names = FALSE, digits = 6)",
  "```"
)
writeLines(qmd, file.path(out_dir, "demo_summary.qmd"))

cat("Parsed", nrow(results), "demo result rows from", length(files), "transcripts\n")
