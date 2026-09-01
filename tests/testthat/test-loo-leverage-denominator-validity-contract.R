loo_denominator_run_local <- function(package, expression) {
  if (identical(package, "npRmpi"))
    getFromNamespace(".npRmpi_with_local_regression", package)(expression)
  else
    expression
}

loo_denominator_bandwidth <- function(package, x, y) {
  owner <- getFromNamespace("npregbw", package)
  loo_denominator_run_local(package, owner(
    xdat = x,
    ydat = y,
    regtype = "lp",
    degree = 1L,
    degree.select = "manual",
    basis = "glp",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckertype = "epanechnikov",
    bws = 1e8,
    bandwidth.compute = FALSE
  ))
}

loo_denominator_objective <- function(package, x, y, bandwidth, largeh) {
  options(np.largeh = largeh)
  evaluate <- getFromNamespace(".npregbw_eval_only", package)
  loo_denominator_run_local(
    package,
    evaluate(x, y, bandwidth, invalid.penalty = "dbmax")
  )
}

test_that("positive-small LOO leverage denominators are never epsilon-floored", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(0, 1e-12, 1))
  y <- as.double(seq_len(nrow(x)))
  bandwidth <- loo_denominator_bandwidth(package, x, y)
  accelerated <- loo_denominator_objective(package, x, y, bandwidth, TRUE)
  canonical <- loo_denominator_objective(package, x, y, bandwidth, FALSE)

  expect_true(is.finite(accelerated$objective))
  expect_gt(accelerated$objective, 1e20)
  expect_equal(accelerated$objective, canonical$objective, tolerance = 2e-14)
  expect_identical(as.integer(accelerated$num.feval.fast), 0L)
})

test_that("ordinary all-large LOO objectives retain fast-path parity", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(0, 0.2, 1))
  y <- as.double(seq_len(nrow(x)))
  bandwidth <- loo_denominator_bandwidth(package, x, y)
  accelerated <- loo_denominator_objective(package, x, y, bandwidth, TRUE)
  canonical <- loo_denominator_objective(package, x, y, bandwidth, FALSE)

  expect_equal(accelerated$objective, canonical$objective, tolerance = 2e-14)
  expect_identical(as.integer(accelerated$num.feval.fast), 1L)
})
