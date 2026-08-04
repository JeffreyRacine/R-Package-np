args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    paste(
      "usage: run_full_mpi_test_shard.R",
      "SHARD SHARD_COUNT SHARD_SIZE NSLAVES WITNESS TEST_DIR"
    )
  )
}

parse_positive_integer <- function(value, where) {
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    stop(where, " must be a positive integer")
  }
  parsed
}

shard <- parse_positive_integer(args[[1L]], "SHARD")
shard_count <- parse_positive_integer(args[[2L]], "SHARD_COUNT")
shard_size <- parse_positive_integer(args[[3L]], "SHARD_SIZE")
nslaves <- parse_positive_integer(args[[4L]], "NSLAVES")
witness <- args[[5L]]
test_dir_path <- normalizePath(args[[6L]], mustWork = TRUE)

if (file.exists(witness)) unlink(witness)

if (shard > shard_count) stop("SHARD exceeds SHARD_COUNT")

suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(npRmpi))

test_files <- sort(list.files(
  test_dir_path, pattern = "^test-.*\\.[rR]$", full.names = FALSE
))
if (!length(test_files) || anyDuplicated(test_files)) {
  stop("npRmpi shard found no unique installed test files")
}

expected_count <- as.integer(ceiling(length(test_files) / shard_size))
if (!identical(shard_count, expected_count)) {
  stop("npRmpi shard count disagrees with the installed test inventory")
}
first <- (shard - 1L) * shard_size + 1L
last <- min(shard * shard_size, length(test_files))
assigned <- test_files[seq.int(first, last)]

regex_escape <- function(value) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", value)
}
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
on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)

npRmpi.init(nslaves = nslaves, quiet = TRUE)
pool_active <- function() {
  isTRUE(getOption("npRmpi.mpi.initialized", FALSE)) &&
    isTRUE(try(mpi.comm.size(1L) == nslaves + 1L, silent = TRUE))
}
if (!pool_active()) stop("npRmpi shard failed to establish its owned pool")

options(cli.hyperlink = FALSE)
withr::local_envvar(TESTTHAT_IS_CHECKING = "true")
results <- test_dir(
  test_dir_path, filter = filter, package = "npRmpi",
  reporter = check_reporter(), load_package = "installed",
  stop_on_failure = FALSE, stop_on_warning = FALSE
)

result_summary <- as.data.frame(results)
if (sum(result_summary$failed) > 0L || any(result_summary$error)) {
  stop("npRmpi shard contains failed expectations")
}

if (!pool_active()) stop("npRmpi shard lost its owned pool")
marker <- sprintf("NP_RMPI_FULL_SHARD_OK %d/%d", shard, shard_count)
writeLines(marker, witness)
cat(marker, "\n", sep = "")
flush.console()
