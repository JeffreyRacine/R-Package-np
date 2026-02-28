#!/usr/bin/env Rscript

parse_args <- function(args) {
  cfg <- list(
    repo = getwd(),
    baseline = NULL,
    tests = character(),
    out_dir = tempfile('testthat_delta_gate_')
  )

  if (length(args) == 0L) return(cfg)

  for (a in args) {
    if (!startsWith(a, '--')) stop('Bad arg: ', a)
    key <- sub('^--([^=]+)=?.*$', '\\1', a)
    has_eq <- grepl('=', a, fixed = TRUE)
    val <- if (has_eq) sub('^--[^=]+=','', a) else 'TRUE'

    if (key == 'repo') cfg$repo <- val
    else if (key == 'baseline') cfg$baseline <- val
    else if (key == 'tests') cfg$tests <- strsplit(val, ',', fixed = TRUE)[[1L]]
    else if (key == 'out_dir') cfg$out_dir <- val
    else stop('Unknown arg: ', key)
  }

  cfg$repo <- normalizePath(cfg$repo, winslash = '/', mustWork = TRUE)
  cfg$out_dir <- normalizePath(cfg$out_dir, winslash = '/', mustWork = FALSE)
  cfg$tests <- trimws(cfg$tests)
  cfg$tests <- cfg$tests[nzchar(cfg$tests)]
  cfg
}

resolve_tests <- function(repo, tests_arg) {
  test_dir <- file.path(repo, 'tests', 'testthat')
  if (!dir.exists(test_dir)) stop('Missing testthat dir: ', test_dir)

  all_files <- list.files(test_dir, pattern = '^test-.*\\.R$', full.names = TRUE)
  if (length(tests_arg) == 0L) return(all_files)

  out <- character()
  for (t in tests_arg) {
    if (file.exists(t)) {
      out <- c(out, normalizePath(t, winslash = '/', mustWork = TRUE))
      next
    }
    cand <- file.path(test_dir, t)
    if (file.exists(cand)) {
      out <- c(out, normalizePath(cand, winslash = '/', mustWork = TRUE))
      next
    }
    stop('Unable to resolve test file: ', t)
  }
  unique(out)
}

is_mpi_sensitive_test <- function(path) {
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  if (!length(txt)) return(FALSE)
  any(grepl("spawn_mpi_slaves\\(|npRmpi\\.init\\(|mpi\\.", txt))
}

read_baseline <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(data.frame(
      file = character(),
      allow_failures = integer(),
      allow_errors = integer(),
      stringsAsFactors = FALSE
    ))
  }

  b <- read.csv(path, stringsAsFactors = FALSE)
  req <- c('file', 'allow_failures', 'allow_errors')
  miss <- setdiff(req, names(b))
  if (length(miss)) stop('Baseline file missing columns: ', paste(miss, collapse = ','))
  b$file <- basename(trimws(b$file))
  b$allow_failures <- as.integer(b$allow_failures)
  b$allow_errors <- as.integer(b$allow_errors)
  b[req]
}

count_expectations <- function(results) {
  failures <- 0L
  errors <- 0L
  warnings <- 0L
  skips <- 0L

  for (tr in results) {
    rr <- tr$results
    if (is.null(rr) || length(rr) == 0L) next
    failures <- failures + sum(vapply(rr, function(x) inherits(x, 'expectation_failure'), logical(1)))
    errors <- errors + sum(vapply(rr, function(x) inherits(x, 'expectation_error'), logical(1)))
    warnings <- warnings + sum(vapply(rr, function(x) inherits(x, 'expectation_warning'), logical(1)))
    skips <- skips + sum(vapply(rr, function(x) inherits(x, 'expectation_skip'), logical(1)))
  }

  list(failures = failures, errors = errors, warnings = warnings, skips = skips)
}

run_one <- function(path) {
  r <- testthat::ListReporter$new()
  run_error <- ''

  tryCatch(
    testthat::test_file(path, reporter = r),
    error = function(e) {
      run_error <<- conditionMessage(e)
    }
  )

  counts <- count_expectations(r$get_results())
  data.frame(
    file = basename(path),
    path = path,
    failures = counts$failures,
    errors = counts$errors,
    warnings = counts$warnings,
    skips = counts$skips,
    run_error = run_error,
    stringsAsFactors = FALSE
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- parse_args(args)
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

  baseline <- read_baseline(cfg$baseline)
  tests <- resolve_tests(cfg$repo, cfg$tests)

  if (basename(cfg$repo) == "np-npRmpi" &&
      !identical(toupper(Sys.getenv("NP_RMPI_ALLOW_INPROC_TESTS", "FALSE")), "TRUE")) {
    sens <- tests[vapply(tests, is_mpi_sensitive_test, logical(1))]
    if (length(sens) > 0L) {
      cat("Refusing in-process testthat run for MPI-sensitive files in npRmpi:\n")
      cat(paste0(" - ", basename(sens), collapse = "\n"), "\n")
      stop("Use session/attach/manual route validators (or set NP_RMPI_ALLOW_INPROC_TESTS=TRUE explicitly).")
    }
  }

  suppressPackageStartupMessages(library(pkgload))
  suppressPackageStartupMessages(library(testthat))
  old_skip <- Sys.getenv("NP_RMPI_SKIP_INIT", unset = NA_character_)
  Sys.setenv(NP_RMPI_SKIP_INIT = "1")
  pkgload::load_all(cfg$repo, quiet = TRUE, export_all = FALSE)
  if (is.na(old_skip)) {
    Sys.unsetenv("NP_RMPI_SKIP_INIT")
  } else {
    Sys.setenv(NP_RMPI_SKIP_INIT = old_skip)
  }

  rows <- lapply(tests, run_one)
  res <- do.call(rbind, rows)

  if (nrow(baseline) > 0L) {
    res <- merge(res, baseline, by = 'file', all.x = TRUE)
  } else {
    res$allow_failures <- NA_integer_
    res$allow_errors <- NA_integer_
  }
  res$allow_failures[is.na(res$allow_failures)] <- 0L
  res$allow_errors[is.na(res$allow_errors)] <- 0L

  res$violate_failures <- res$failures > res$allow_failures
  res$violate_errors <- res$errors > res$allow_errors
  res$violate_run_error <- nzchar(res$run_error)
  res$status <- ifelse(
    res$violate_failures | res$violate_errors | res$violate_run_error,
    'FAIL',
    'PASS'
  )

  summary_path <- file.path(cfg$out_dir, 'testthat_delta_summary.csv')
  write.csv(res[order(res$file), ], summary_path, row.names = FALSE)

  cat('repo:', cfg$repo, '\n')
  cat('baseline:', if (is.null(cfg$baseline)) '<none>' else cfg$baseline, '\n')
  cat('summary:', summary_path, '\n')
  cat('status table:\n')
  print(table(res$status, useNA = 'ifany'))

  bad <- res[res$status == 'FAIL', c('file', 'failures', 'allow_failures', 'errors', 'allow_errors', 'run_error')]
  if (nrow(bad) > 0L) {
    cat('\nviolations:\n')
    print(bad, row.names = FALSE)
    stop('testthat delta gate failed')
  }

  cat('TESTTHAT_DELTA_GATE_OK\n')
}

main()
