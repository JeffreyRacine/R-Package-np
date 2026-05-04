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
    family = NA_character_,
    case = NA_character_,
    tier = NA_character_,
    regtype = NA_character_,
    nomad = NA_character_,
    degree = NA_character_,
    selected.degree = NA_real_,
    bwtype = NA_character_,
    objective = NA_real_,
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
        family = if (!is.null(kv$family)) kv$family else NA_character_,
        case = if (!is.null(kv$case)) kv$case else NA_character_,
        tier = if (!is.null(kv$tier)) kv$tier else NA_character_,
        regtype = if (!is.null(kv$regtype)) kv$regtype else NA_character_,
        nomad = if (!is.null(kv$nomad)) kv$nomad else NA_character_,
        degree = if (!is.null(kv$degree)) kv$degree else NA_character_,
        selected.degree = num_or_na(kv$selected.degree),
        bwtype = if (!is.null(kv$bwtype)) kv$bwtype else NA_character_,
        objective = num_or_na(kv$objective),
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
             bwmethod = character(), family = character(), case = character(),
             tier = character(), regtype = character(), nomad = character(),
             degree = character(), selected.degree = numeric(),
             bwtype = character(), objective = numeric(), num.feval = numeric(),
             num.feval.fast = numeric(), nomad.time = numeric(),
             powell.time = numeric(), fitted.length = integer(),
             source_file = character(),
             stringsAsFactors = FALSE)
}

for (name in c("family", "case", "tier", "regtype", "nomad", "degree",
               "selected.degree", "bwtype", "objective", "num.feval",
               "num.feval.fast", "nomad.time", "powell.time",
               "fitted.length")) {
  if (!name %in% names(results)) results[[name]] <- NA
}

results <- results[order(results$family, results$demo, results$case,
                         results$mode, results$ranks, results$slaves,
                         results$source_file), ]

write.csv(results, file.path(out_dir, "demo_results.csv"), row.names = FALSE)

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

tex_escape <- function(x) gsub("_", "\\\\_", x, fixed = TRUE)
tex <- c(
  "\\begin{tabular}{llrrrr}",
  "Demo & Mode & Ranks & Slaves & $n$ & Seconds\\\\",
  "\\hline"
)
if (nrow(results)) {
  for (i in seq_len(nrow(results))) {
    row <- results[i, ]
    tex <- c(tex, sprintf("%s & %s & %s & %s & %s & %.3f\\\\",
                          tex_escape(row$demo), tex_escape(row$mode),
                          ifelse(is.na(row$ranks), "--", row$ranks),
                          ifelse(is.na(row$slaves), "--", row$slaves),
                          ifelse(is.na(row$n), "--", row$n),
                          row$elapsed))
  }
}
tex <- c(tex, "\\hline", "\\end{tabular}")
writeLines(tex, file.path(out_dir, "timing_all.tex"))

cat("Parsed", nrow(results), "demo result rows from", length(files), "transcripts\n")
