#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("usage: validate_demo_matrix.R <demo_results.csv> <observed.tsv> <status.tsv> <require_full_matrix>",
       call. = FALSE)
}

results_path <- args[[1L]]
observed_path <- args[[2L]]
status_path <- args[[3L]]
require_full_matrix <- tolower(args[[4L]]) %in% c("1", "true", "yes")

write_status <- function(effective_status, proof, details = "") {
  dir.create(dirname(status_path), showWarnings = FALSE, recursive = TRUE)
  lines <- c(
    paste("effective_status", effective_status, sep = "\t"),
    paste("require_full_matrix", require_full_matrix, sep = "\t"),
    paste("proof", proof, sep = "\t"),
    paste("details", details, sep = "\t")
  )
  writeLines(lines, status_path)
}

if (!file.exists(results_path)) {
  write_status("FAIL", "missing_demo_results_csv", results_path)
  quit(status = 1L)
}

results <- read.csv(results_path, stringsAsFactors = FALSE)
required <- c("mode", "ranks", "slaves")
missing <- setdiff(required, names(results))
if (length(missing)) {
  write_status("FAIL", "demo_results_missing_required_columns",
               paste(missing, collapse = ","))
  quit(status = 1L)
}

if (!nrow(results)) {
  write_status("FAIL", "demo_results_has_no_rows", results_path)
  quit(status = 1L)
}

norm_int <- function(x) {
  out <- suppressWarnings(as.integer(x))
  out
}

results$ranks <- norm_int(results$ranks)
results$slaves <- norm_int(results$slaves)

key <- unique(results[, c("mode", "ranks", "slaves"), drop = FALSE])
key$rows <- vapply(seq_len(nrow(key)), function(i) {
  sum(results$mode == key$mode[[i]] &
        ((is.na(results$ranks) & is.na(key$ranks[[i]])) |
           (!is.na(results$ranks) & !is.na(key$ranks[[i]]) &
              results$ranks == key$ranks[[i]])) &
        ((is.na(results$slaves) & is.na(key$slaves[[i]])) |
           (!is.na(results$slaves) & !is.na(key$slaves[[i]]) &
              results$slaves == key$slaves[[i]])))
}, integer(1L))
key <- key[order(key$mode, key$ranks, key$slaves), , drop = FALSE]

dir.create(dirname(observed_path), showWarnings = FALSE, recursive = TRUE)
write.table(key, observed_path, sep = "\t", quote = FALSE, row.names = FALSE)

problems <- character()
if (!"serial" %in% results$mode) {
  problems <- c(problems, "missing serial rows")
}

if (require_full_matrix) {
  session_slaves <- sort(unique(results$slaves[results$mode == "session"]))
  missing_session <- setdiff(c(1L, 2L, 3L), session_slaves)
  if (length(missing_session)) {
    problems <- c(problems,
                  paste0("missing session nslaves=", paste(missing_session, collapse = ",")))
  }

  attach_ranks <- sort(unique(results$ranks[results$mode == "attach"]))
  attach_ranks <- attach_ranks[!is.na(attach_ranks)]
  if (length(attach_ranks) < 2L) {
    problems <- c(problems, "attach has fewer than two rank counts")
  }

  profile_ranks <- sort(unique(results$ranks[results$mode == "profile"]))
  profile_ranks <- profile_ranks[!is.na(profile_ranks)]
  if (length(profile_ranks) < 2L) {
    problems <- c(problems, "profile has fewer than two rank counts")
  }
}

if (length(problems)) {
  write_status("FAIL", "matrix_coverage_failed", paste(problems, collapse = "; "))
  quit(status = 1L)
}

write_status("PASS", "matrix_coverage_parsed", sprintf("rows=%d", nrow(results)))
