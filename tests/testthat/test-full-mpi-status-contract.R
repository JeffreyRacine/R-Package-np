test_that("full-suite raw receipts preserve exits independently of witnesses", {
  helper <- testthat::test_path("..", "validation", "full_mpi_test_status.R")
  skip_if_not(file.exists(helper), "source status helper unavailable")
  scope <- new.env(parent = globalenv())
  sys.source(helper, envir = scope)
  directory <- withr::local_tempdir()
  for (status in c(0L, 7L, 124L)) {
    witness <- file.path(directory, paste0("status-", status, ".ok"))
    output <- capture.output(scope$npRmpi_full_test_record_status(
      witness, 2L, 4L, status
    ))
    expected <- sprintf("NP_RMPI_FULL_SHARD_RAW_STATUS 2/4 status=%d", status)
    expect_identical(output, expected)
    expect_identical(readLines(paste0(witness, ".raw-status")), expected)
    expect_false(file.exists(witness))
    expect_false(file.exists(paste0(witness, ".raw-status.tmp")))
  }
  before <- list.files(directory, all.files = TRUE)
  for (bad in list(NA_integer_, -1L, 0.5, Inf, 2^32, integer(), c(0L, 1L))) {
    expect_error(scope$npRmpi_full_test_record_status(
      file.path(directory, "invalid.ok"), 1L, 4L, bad
    ), "invalid npRmpi full-suite raw-status identity", fixed = TRUE)
  }
  expect_error(scope$npRmpi_full_test_record_status(
    file.path(directory, "invalid.ok"), 5L, 4L, 0L
  ), "invalid npRmpi full-suite raw-status identity", fixed = TRUE)
  expect_identical(list.files(directory, all.files = TRUE), before)
})

test_that("raw receipts survive a later failing parent exit", {
  helper <- testthat::test_path("..", "validation", "full_mpi_test_status.R")
  skip_if_not(file.exists(helper), "source status helper unavailable")
  directory <- withr::local_tempdir()
  witness <- file.path(directory, "survives.ok")
  result <- npRmpi_run_rscript_subprocess(c(
    sprintf("source(%s)", deparse(normalizePath(helper))),
    sprintf("npRmpi_full_test_record_status(%s, 1L, 2L, 0L)", deparse(witness)),
    "quit(save='no', status=9L, runLast=FALSE)"
  ), cleanup = FALSE)
  expect_identical(result$status, 9L)
  expect_identical(readLines(paste0(witness, ".raw-status")),
                   "NP_RMPI_FULL_SHARD_RAW_STATUS 1/2 status=0")
  expect_false(file.exists(witness))
})

test_that("the full-suite parent records direct exits before aggregate coercion", {
  parent <- testthat::test_path("..", "testthat.R")
  skip_if_not(file.exists(parent), "source parent unavailable")
  lines <- readLines(parent, warn = FALSE)
  launch <- grep("statuses[[shard]] <- system2(", lines, fixed = TRUE)
  record <- grep("npRmpi_full_test_record_status(", lines, fixed = TRUE)
  aggregate <- grep("if (!witnessed || statuses[[shard]] != 0L)", lines, fixed = TRUE)
  clear <- grep("if (file.exists(raw_receipt)) unlink(raw_receipt)", lines, fixed = TRUE)
  expect_length(launch, 1L)
  expect_length(record, 1L)
  expect_length(aggregate, 1L)
  expect_length(clear, 1L)
  expect_lt(clear, launch)
  expect_lt(launch, record)
  expect_lt(record, aggregate)
  expect_match(paste(lines[seq.int(record, record+2L)], collapse = "\n"),
               "witness, shard, shard_count, statuses[[shard]]", fixed = TRUE)
})
