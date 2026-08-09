parse_positive_integer <- function(value, where) {
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    stop(where, " must be a positive integer")
  }
  parsed
}

regex_escape <- function(value) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", value)
}

npRmpi_shard_pool_active <- function(nslaves) {
  isTRUE(getOption("npRmpi.mpi.initialized", FALSE)) &&
    isTRUE(try(mpi.comm.size(1L) == nslaves + 1L, silent = TRUE))
}

npRmpi_shard_has_failures <- function(result_summary) {
  if (!is.data.frame(result_summary) || !nrow(result_summary))
    stop("npRmpi shard produced no test results")
  required <- c("failed", "error")
  missing <- setdiff(required, names(result_summary))
  if (length(missing)) {
    stop(
      "npRmpi shard result schema is missing: ",
      paste(missing, collapse = ", ")
    )
  }

  failed <- suppressWarnings(as.integer(result_summary[["failed"]]))
  errors <- as.logical(result_summary[["error"]])
  if (anyNA(failed) || any(failed < 0L) || anyNA(errors))
    stop("npRmpi shard result summary contains indeterminate status")
  any(failed > 0L) || any(errors)
}

npRmpi_run_full_test_shard <- function(args) {
  if (length(args) != 6L) {
    stop(
      paste(
        "usage: run_full_mpi_test_shard.R",
        "SHARD SHARD_COUNT SHARD_SIZE NSLAVES WITNESS TEST_DIR"
      )
    )
  }

  shard <- parse_positive_integer(args[[1L]], "SHARD")
  shard_count <- parse_positive_integer(args[[2L]], "SHARD_COUNT")
  shard_size <- parse_positive_integer(args[[3L]], "SHARD_SIZE")
  nslaves <- parse_positive_integer(args[[4L]], "NSLAVES")
  witness_dir <- normalizePath(dirname(args[[5L]]), mustWork = TRUE)
  witness <- file.path(witness_dir, basename(args[[5L]]))
  witness_tmp <- paste0(witness, ".tmp-", Sys.getpid())
  test_dir_path <- normalizePath(args[[6L]], mustWork = TRUE)
  pool_owned <- FALSE
  success <- FALSE

  unlink(c(witness, witness_tmp), force = TRUE)
  on.exit({
    if (!success)
      unlink(c(witness, witness_tmp), force = TRUE)
    if (pool_owned)
      try(npRmpi.quit(force = TRUE), silent = TRUE)
  }, add = TRUE)

  if (shard > shard_count)
    stop("SHARD exceeds SHARD_COUNT")

  suppressPackageStartupMessages(library(testthat))
  suppressPackageStartupMessages(library(npRmpi))

  test_files <- sort(list.files(
    test_dir_path, pattern = "^test-.*\\.[rR]$", full.names = FALSE
  ))
  if (!length(test_files) || anyDuplicated(test_files))
    stop("npRmpi shard found no unique installed test files")

  expected_count <- as.integer(ceiling(length(test_files) / shard_size))
  if (!identical(shard_count, expected_count))
    stop("npRmpi shard count disagrees with the installed test inventory")
  first <- (shard - 1L) * shard_size + 1L
  last <- min(shard * shard_size, length(test_files))
  assigned <- test_files[seq.int(first, last)]

  stems <- sub("\\.[rR]$", "", sub("^test-", "", assigned))
  filter <- paste0("^(", paste(regex_escape(stems), collapse = "|"), ")$")

  cat(
    sprintf(
      "NP_RMPI_FULL_SHARD_START %d/%d files=%d first=%s last=%s\n",
      shard, shard_count, length(assigned), assigned[[1L]],
      assigned[[length(assigned)]]
    )
  )
  flush.console()

  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)
  Sys.setenv(NP_RMPI_TEST_SUITE_POOL = "1")
  pool_owned <- TRUE
  npRmpi.init(nslaves = nslaves, quiet = TRUE)
  if (!npRmpi_shard_pool_active(nslaves))
    stop("npRmpi shard failed to establish its owned pool")

  options(cli.hyperlink = FALSE)
  withr::local_envvar(TESTTHAT_IS_CHECKING = "true")
  results <- test_dir(
    test_dir_path, filter = filter, package = "npRmpi",
    reporter = check_reporter(), load_package = "installed",
    stop_on_failure = FALSE, stop_on_warning = FALSE
  )

  result_summary <- as.data.frame(results)
  if (npRmpi_shard_has_failures(result_summary))
    stop("npRmpi shard contains failed expectations")
  if (!npRmpi_shard_pool_active(nslaves))
    stop("npRmpi shard lost its owned pool")

  cleanup_error <- tryCatch({
    npRmpi.quit(force = TRUE)
    NULL
  }, error = identity)
  if (inherits(cleanup_error, "error"))
    stop("npRmpi shard pool cleanup failed: ", conditionMessage(cleanup_error))
  if (npRmpi_shard_pool_active(nslaves))
    stop("npRmpi shard pool remained active after cleanup")
  pool_owned <- FALSE

  marker <- sprintf("NP_RMPI_FULL_SHARD_OK %d/%d", shard, shard_count)
  writeLines(marker, witness_tmp, useBytes = TRUE)
  if (!identical(readLines(witness_tmp, warn = FALSE), marker))
    stop("npRmpi shard temporary witness verification failed")
  if (!file.rename(witness_tmp, witness))
    stop("npRmpi shard could not publish its success witness")
  if (!identical(readLines(witness, warn = FALSE), marker))
    stop("npRmpi shard published an invalid success witness")

  success <- TRUE
  cat(marker, "\n", sep = "")
  flush.console()
  0L
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 5L && nzchar(args[[5L]]))
  unlink(c(args[[5L]], paste0(args[[5L]], ".tmp-", Sys.getpid())), force = TRUE)

exit_status <- tryCatch(
  npRmpi_run_full_test_shard(args),
  error = function(error) {
    message("Error: ", conditionMessage(error))
    1L
  }
)
quit(save = "no", status = exit_status, runLast = FALSE)
