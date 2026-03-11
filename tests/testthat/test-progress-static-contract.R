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
