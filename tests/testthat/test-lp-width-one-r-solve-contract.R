test_that("shared evaluation solver dispatches width one before solve", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  solve_eval <- getFromNamespace(".npreghat_solve_eval", package)
  W <- matrix(c(1.0, 0.5, 2.0, 1.5), ncol = 1L)
  k <- c(0.8, 0.2, 0.7, 0.4)
  w.eval <- 1.25

  result <- solve_eval(W = W, w.eval = w.eval, k = k, ridge.base = 0.0)
  denominator <- drop(crossprod(W, W * k))
  expect_identical(result$v, w.eval / denominator)
  expect_identical(result$ridge, 0.0)

  source_file <- testthat::test_path("..", "..", "R", "np.reghat.R")
  skip_if_not(file.exists(source_file), "source R/np.reghat.R unavailable")
  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  start <- regexpr(".npreghat_solve_eval <- function(", source, fixed = TRUE)
  stop <- regexpr(
    ".npreghat_exact_matrix_from_core <- function(",
    source,
    fixed = TRUE
  )
  expect_gt(start, 0L)
  expect_gt(stop, start)
  region <- substr(source, start, stop - 1L)
  expect_true(grepl("if (p == 1L)", region, fixed = TRUE))
  expect_true(grepl("v <- w.eval[1L] / denominator", region, fixed = TRUE))
  expect_lt(
    regexpr("if (p == 1L)", region, fixed = TRUE),
    regexpr("solve(t(A)", region, fixed = TRUE)
  )
})
