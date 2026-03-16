test_that("progress core avoids legacy emitters", {
  src_path <- testthat::test_path("..", "..", "R", "progress.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_render_single_line <- function\\(", src))
  expect_true(grepl("base::cat\\(", src))
  expect_false(grepl("\\bprint\\(", src))
  expect_true(grepl("\\.np_progress_select_bandwidth <- function\\(", src))
  expect_true(grepl("\\.np_progress_select_bandwidth_enhanced <- function\\(", src))
  expect_true(grepl("\\.np_progress_select_iv <- function\\(", src))
  expect_true(grepl("\\.np_warning <- function\\(", src))
})

test_that("R layer routes warnings through unified helper", {
  r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "source R files unavailable in installed test context")
  files <- setdiff(list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE), file.path(r_dir, "progress.R"))

  for (src_path in files) {
    src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
    expect_false(grepl("\\bwarning\\(", src), info = basename(src_path))
    expect_false(grepl("\\bmessage\\(", src), info = basename(src_path))
  }
})

test_that("core estimator wrappers emit top-level bandwidth-selection notes", {
  cases <- list(
    list(file = "np.regression.R", label = "Selecting regression bandwidth", helper = ".np_progress_select_bandwidth_enhanced"),
    list(file = "np.density.R", label = "Selecting density bandwidth", helper = ".np_progress_select_bandwidth_enhanced"),
    list(file = "np.distribution.R", label = "Selecting distribution bandwidth", helper = ".np_progress_select_bandwidth_enhanced"),
    list(file = "np.condensity.R", label = "Selecting conditional density bandwidth", helper = ".np_progress_select_bandwidth_enhanced"),
    list(file = "np.condistribution.R", label = "Selecting conditional distribution bandwidth", helper = ".np_progress_select_bandwidth_enhanced"),
    list(file = "np.plregression.R", label = "Selecting partially linear regression bandwidth", helper = ".np_progress_select_bandwidth_enhanced"),
    list(file = "np.singleindex.R", label = "Selecting single-index bandwidth", helper = ".np_progress_select_bandwidth_enhanced"),
    list(file = "np.smoothcoef.R", label = "Selecting smooth coefficient bandwidth", helper = ".np_progress_select_bandwidth_enhanced"),
    list(file = "np.conmode.R", label = "Selecting conditional density bandwidth"),
    list(file = "np.qregression.R", label = "Selecting conditional distribution bandwidth")
  )

  for (case in cases) {
    src_path <- testthat::test_path("..", "..", "R", case$file)
    skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
    src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

    helper <- if (is.null(case$helper)) ".np_progress_select_bandwidth" else case$helper
    expect_true(
      grepl(sprintf("%s\\(\\s*\"%s\"", helper, case$label), src),
      info = case$file
    )
  }
})

test_that("npregivderiv no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "npregivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Iterating Landweber-Fridman derivative solve\"", src))
})

test_that("npregiv no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "npregiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_select_iv\\(", src))
  expect_true(grepl("\\.np_progress_iv_set_object\\(", src))
  expect_true(grepl("\\.np_progress_iv_activity_step\\(", src))
  expect_true(grepl("\\.np_progress_with_legacy_suppressed\\(", src))
})

test_that("npunitest no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.unitest.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Computing bandwidths\"\\)", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Bootstrap replications\"", src))
})

test_that("npdeneqtest no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.deneqtest.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Computing bandwidths\"\\)", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Bootstrap replications\"", src))
  expect_true(grepl("\\.np_progress_with_legacy_suppressed\\(", src))
})

test_that("npdeptest no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.deptest.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Computing bandwidths\"\\)", src))
  expect_true(grepl("\\.np_progress_note\\(\"Constructing metric entropy\"\\)", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Bootstrap replications\"", src))
  expect_true(grepl("\\.np_progress_with_legacy_suppressed\\(", src))
})

test_that("npcmstest no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.cmstest.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Computing bandwidths\"\\)", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Bootstrap replications\"", src))
  expect_true(grepl("\\.np_progress_with_legacy_suppressed\\(", src))
})

test_that("npqcmstest no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.qcmstest.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_false(grepl("\\bcat\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Computing bandwidths\"\\)", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Bootstrap replications\"", src))
  expect_true(grepl("\\.np_progress_with_legacy_suppressed\\(", src))
})

test_that("npsymtest no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.symtest.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Computing bandwidths\"\\)", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Bootstrap replications\"", src))
})

test_that("npsdeptest no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.sdeptest.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Computing bandwidths\"\\)", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Constructing metric entropy by lag\"", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Bootstrap replications\"", src))
  expect_true(grepl("\\.np_progress_with_legacy_suppressed\\(", src))
})

test_that("npsigtest no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.sigtest.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Testing joint significance\"\\)", src))
  expect_true(grepl("\\.np_progress_note\\(sprintf\\(\"Testing variable %s of \\(%s\\)\",", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Bootstrap replications\"", src))
  expect_true(grepl("\\.np_progress_with_legacy_suppressed\\(", src))
})

test_that("npcopula no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.copula.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_note\\(\"Computing the copula for the sample realizations\"\\)", src))
  expect_true(grepl("\\.np_progress_note\\(\"Expanding the u matrix\"\\)", src))
  expect_true(grepl("\\.np_progress_note\\(\\s*sprintf\\(\\s*\"Computing the marginal of %s for the sample realizations\"", src))
  expect_true(grepl("\\.np_progress_note\\(\\s*sprintf\\(\\s*\"Computing the quasi-inverse for the marginal of %s\"", src))
  expect_true(grepl("\\.np_progress_note\\(\\s*sprintf\\(\\s*\"Computing the marginal of %s for the expanded grid\"", src))
})

test_that("npindexbw no longer uses legacy console helpers", {
  src_path <- testthat::test_path("..", "..", "R", "np.singleindex.bw.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_select_bandwidth_enhanced\\(\\s*\"Selecting single-index bandwidth\"", src))
  expect_true(grepl("\\.np_progress_bandwidth_activity_step\\(done = bandwidth_eval_count\\)", src))
  expect_true(grepl("\\.np_progress_bandwidth_multistart_step\\(done = i, total = nmulti\\)", src))
})

test_that("npscoefbw multistart path uses central progress core", {
  src_path <- testthat::test_path("..", "..", "R", "np.smoothcoef.bw.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("Optimizing smooth coefficient bandwidth", src, fixed = TRUE))
  expect_true(grepl("\\.np_progress_select_bandwidth_enhanced\\(\\s*\"Selecting smooth coefficient bandwidth\"", src))
  expect_true(grepl("\\.np_progress_bandwidth_activity_step\\(done = cv_state\\$optim_eval\\)", src))
  expect_true(grepl("\\.np_progress_bandwidth_multistart_step\\(done = i, total = nmulti\\)", src))
  expect_true(grepl("Backfitting smooth coefficient bandwidth", src, fixed = TRUE))
  expect_true(grepl("Optimizing partial residual bandwidth", src, fixed = TRUE))
  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
})

test_that("npplregbw direct path uses coordinated bandwidth selector", {
  src_path <- testthat::test_path("..", "..", "R", "np.plregression.bw.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("\\.np_progress_bandwidth_set_coordinator\\(", src))
  expect_true(grepl("\\.np_progress_bandwidth_set_coordinator_group\\(1L, \"y~z\"\\)", src))
  expect_true(grepl("\\.np_progress_select_bandwidth_enhanced\\(\\s*\"Selecting partially linear regression bandwidth\"", src))
})

test_that("compiled multistart progress uses shared bandwidth helper", {
  src_path <- testthat::test_path("..", "..", "src", "np.c")
  skip_if_not(file.exists(src_path), "source C files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("np_progress_bandwidth_multistart_step\\(", src))
  expect_true(grepl("np_progress_bandwidth_activity_step\\(", src))
  expect_true(grepl("bwm_maybe_signal_activity\\(", src))
  expect_true(grepl("np_progress_signal\\(\"bandwidth_activity_step\",\\s*\"bandwidth\"", src))
  expect_false(grepl("Rprintf\\(\"\\\\rMultistart", src))
  expect_false(grepl("Rprintf\\(\"\\\\r +\\\\r\"", src))
})

test_that("plot helpers use append-only progress core", {
  src_path <- testthat::test_path("..", "..", "R", "np.plot.helpers.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_begin\\(", src))
  expect_true(grepl("\\.np_progress_step\\(", src))
  expect_true(grepl("\\.np_progress_end\\(", src))
})
