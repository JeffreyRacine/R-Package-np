test_that("progress core avoids legacy emitters", {
  src_path <- testthat::test_path("..", "..", "R", "progress.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_false(grepl("\\bcat\\(", src))
  expect_false(grepl("\\bprint\\(", src))
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
  expect_true(grepl("\\.np_progress_begin\\(\"Iterating Landweber-Fridman solve\"", src))
  expect_true(grepl("\\.np_progress_begin\\(\"Iterating Tikhonov solve\"", src))
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
  expect_true(grepl("\\.np_progress_begin\\(\"Multistart optimization\"", src))
  expect_true(grepl("detail = sprintf\\(\"multistart %d\", i\\)", src))
})

test_that("npscoefbw multistart path uses append-only progress core", {
  src_path <- testthat::test_path("..", "..", "R", "np.smoothcoef.bw.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_true(grepl("Optimizing smooth coefficient bandwidth", src, fixed = TRUE))
  expect_true(grepl("Multistart optimization", src, fixed = TRUE))
})

test_that("plot helpers use append-only progress core", {
  src_path <- testthat::test_path("..", "..", "R", "np.plot.helpers.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_false(grepl("printPush\\(|printPop\\(|printClear\\(|newLineConsole\\(", src))
  expect_true(grepl("\\.np_progress_begin\\(", src))
  expect_true(grepl("\\.np_progress_step\\(", src))
  expect_true(grepl("\\.np_progress_end\\(", src))
  expect_true(grepl("\\.np_progress_note\\(", src))
})
