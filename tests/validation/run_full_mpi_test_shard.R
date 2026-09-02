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

npRmpi_shard_open_pool <- function(nslaves) {
  npRmpi.init(nslaves = nslaves, quiet = TRUE)
  if (!npRmpi_shard_pool_active(nslaves))
    stop("npRmpi shard failed to establish its owned pool")
  invisible(TRUE)
}

npRmpi_shard_close_pool <- function(nslaves) {
  cleanup_error <- tryCatch({
    npRmpi.quit(force = TRUE)
    NULL
  }, error = identity)
  if (inherits(cleanup_error, "error"))
    stop("npRmpi shard pool cleanup failed: ", conditionMessage(cleanup_error))
  if (npRmpi_shard_pool_active(nslaves))
    stop("npRmpi shard pool remained active after cleanup")
  invisible(TRUE)
}

npRmpi_shard_local_only_files <- function() {
  c(
    "test-adaptive-conditional-plot-session-contract.R",
    "test-attach-gennn-exdat-contract.R",
    "test-mpi-finalization-lifecycle-contract.R"
  )
}

npRmpi_shard_local_mode_files <- function() {
  c(
    "test-adaptive-nn-exact-conditional-density-cvls-contract.R",
    "test-adaptive-nn-exact-conditional-density-cvml-contract.R",
    "test-adaptive-nn-exact-conditional-distribution-cv-contract.R",
    "test-adaptive-nn-exact-density-cv-contract.R",
    "test-adaptive-nn-exact-distribution-cv-contract.R",
    "test-adaptive-nn-exact-regression-cvls-contract.R",
    "test-adaptive-nn-exact-regression-deleteone-objectives-contract.R",
    "test-beta-bandwidth-objectives-contract.R",
    "test-beta-conditional-contract.R",
    "test-beta-conditional-count-canonical-contract.R",
    "test-beta-gradient-estimator-contract.R",
    "test-beta-higher-order-regression-contract.R",
    "test-beta-regression-objective-canonical-contract.R"
  )
}

npRmpi_shard_run_files <- function(test_dir_path, test_files) {
  if (!length(test_files))
    return(NULL)
  stems <- sub("\\.[rR]$", "", sub("^test-", "", test_files))
  filter <- paste0(
    "^(", paste(regex_escape(stems), collapse = "|"), ")$"
  )
  reporter <- check_reporter()
  test_dir(
    test_dir_path,
    filter = filter,
    package = "npRmpi",
    reporter = reporter,
    load_package = "installed",
    stop_on_failure = FALSE,
    stop_on_warning = FALSE
  )
}

npRmpi_shard_run_isolated_files <- function(test_dir_path, test_files) {
  if (!length(test_files))
    return(NULL)
  result.frames <- vector("list", length(test_files))
  for (index in seq_along(test_files)) {
    test.file <- test_files[[index]]
    started <- proc.time()[["elapsed"]]
    cat(sprintf(
      "NP_RMPI_LOCAL_FILE_START %d/%d file=%s\n",
      index, length(test_files), test.file
    ))
    flush.console()
    result <- test_file(
      file.path(test_dir_path, test.file),
      package = "npRmpi",
      reporter = check_reporter(),
      load_package = "installed",
      stop_on_failure = FALSE,
      stop_on_warning = FALSE
    )
    frame <- as.data.frame(result)
    if (npRmpi_shard_has_failures(frame))
      stop("npRmpi local-only file contains failed expectations: ", test.file)
    result.frames[[index]] <- frame
    elapsed <- proc.time()[["elapsed"]] - started
    cat(sprintf(
      "NP_RMPI_LOCAL_FILE_OK %d/%d file=%s elapsed=%.3f\n",
      index, length(test_files), test.file, elapsed
    ))
    flush.console()
  }
  do.call(rbind, result.frames)
}

npRmpi_shard_run_local_mode_files <- function(test_dir_path, test_files) {
  withr::local_envvar(NP_RMPI_TEST_SUITE_LOCAL_MODE = "1")
  old.local <- getOption("npRmpi.local.regression.mode", FALSE)
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.context <- getOption("npRmpi.autodispatch.context", FALSE)
  old.native <- FALSE
  native.active <- FALSE
  on.exit({
    if (native.active) {
      try(.Call(
        "C_np_set_local_regression_mode", old.native,
        PACKAGE = "npRmpi"
      ), silent = TRUE)
    }
    options(
      npRmpi.local.regression.mode = old.local,
      npRmpi.autodispatch.disable = old.disable,
      npRmpi.autodispatch.context = old.context
    )
  }, add = TRUE)

  options(
    npRmpi.local.regression.mode = TRUE,
    npRmpi.autodispatch.disable = TRUE,
    npRmpi.autodispatch.context = TRUE
  )
  old.native <- .Call(
    "C_np_set_local_regression_mode", TRUE, PACKAGE = "npRmpi"
  )
  native.active <- TRUE
  npRmpi_shard_run_files(test_dir_path, test_files)
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

  plan_file <- normalizePath(file.path(
    dirname(test_dir_path), "validation", "full_mpi_test_plan.R"
  ), mustWork = TRUE)
  source(plan_file, local = TRUE)
  execution_plan <- npRmpi_full_test_plan(test_files, shard_size)
  expected_count <- length(execution_plan)
  if (!identical(shard_count, expected_count))
    stop("npRmpi shard count disagrees with the installed test inventory")
  assigned <- execution_plan[[shard]]
  local_registry <- npRmpi_shard_local_only_files()
  local_mode_registry <- npRmpi_shard_local_mode_files()
  if (anyDuplicated(c(local_registry, local_mode_registry)) ||
      !all(c(local_registry, local_mode_registry) %in% test_files))
    stop("npRmpi shard execution-owner registry is invalid")
  local_assigned <- intersect(assigned, local_registry)
  local_mode_assigned <- intersect(assigned, local_mode_registry)
  pooled_assigned <- setdiff(
    assigned, c(local_assigned, local_mode_assigned)
  )
  if (!setequal(
    assigned, c(local_assigned, local_mode_assigned, pooled_assigned)
  ))
    stop("npRmpi shard lane partition lost or duplicated a test file")

  cat(
    sprintf(
      "NP_RMPI_FULL_SHARD_START %d/%d files=%d first=%s last=%s\n",
      shard, shard_count, length(assigned), assigned[[1L]],
      assigned[[length(assigned)]]
    )
  )
  flush.console()

  options(cli.hyperlink = FALSE)
  withr::local_envvar(TESTTHAT_IS_CHECKING = "true")
  if (length(pooled_assigned)) {
    options(npRmpi.autodispatch = TRUE, np.messages = FALSE)
    Sys.setenv(NP_RMPI_TEST_SUITE_POOL = "1")
    pool_owned <- TRUE
    npRmpi_shard_open_pool(nslaves)
    pooled_results <- npRmpi_shard_run_files(
      test_dir_path, pooled_assigned
    )
    if (npRmpi_shard_has_failures(as.data.frame(pooled_results)))
      stop("npRmpi shard pooled lane contains failed expectations")
    if (!npRmpi_shard_pool_active(nslaves))
      stop("npRmpi shard lost its owned pool")
    npRmpi_shard_close_pool(nslaves)
    pool_owned <- FALSE
  }

  if (length(local_mode_assigned)) {
    if (length(pooled_assigned)) {
      cat(sprintf(
        "NP_RMPI_FULL_SHARD_POOL_RESTART %d/%d\n", shard, shard_count
      ))
      flush.console()
    }
    options(npRmpi.autodispatch = TRUE, np.messages = FALSE)
    Sys.setenv(NP_RMPI_TEST_SUITE_POOL = "1")
    pool_owned <- TRUE
    npRmpi_shard_open_pool(nslaves)
    cat(sprintf(
      "NP_RMPI_FULL_SHARD_LOCAL_MODE %d/%d files=%s\n",
      shard, shard_count, paste(local_mode_assigned, collapse = ",")
    ))
    flush.console()
    local_mode_results <- npRmpi_shard_run_local_mode_files(
      test_dir_path, local_mode_assigned
    )
    if (npRmpi_shard_has_failures(as.data.frame(local_mode_results)))
      stop("npRmpi shard local-mode lane contains failed expectations")
    if (!npRmpi_shard_pool_active(nslaves))
      stop("npRmpi shard local-mode lane damaged its owned pool")
    npRmpi_shard_close_pool(nslaves)
    pool_owned <- FALSE
  }

  if (length(local_assigned)) {
    cat(sprintf(
      "NP_RMPI_FULL_SHARD_LOCAL_ONLY %d/%d files=%s\n",
      shard, shard_count, paste(local_assigned, collapse = ",")
    ))
    flush.console()
    local_results <- withr::with_envvar(
      c(NP_RMPI_TEST_SUITE_POOL = ""),
      npRmpi_shard_run_isolated_files(test_dir_path, local_assigned)
    )
    if (npRmpi_shard_has_failures(as.data.frame(local_results)))
      stop("npRmpi shard local-only lane contains failed expectations")
    if (npRmpi_shard_pool_active(nslaves))
      stop("npRmpi shard local-only lane established a slave pool")
  }

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
