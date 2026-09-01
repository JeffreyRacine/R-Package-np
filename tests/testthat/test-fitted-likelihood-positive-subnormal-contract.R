fitted_subnormal_count_fixed <- function(text, pattern) {
  match <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  if (identical(match, -1L)) 0L else length(match)
}

fitted_subnormal_run_local <- function(package, expression) {
  if (identical(package, "npRmpi"))
    getFromNamespace(".npRmpi_with_local_regression", package)(expression)
  else
    expression
}

test_that("one fitted-likelihood helper owns the complete current floor sweep", {
  header <- paste(readLines(test_path("..", "..", "src", "headers.h"),
                            warn = FALSE), collapse = "\n")
  jksum <- paste(readLines(test_path("..", "..", "src", "jksum.c"),
                            warn = FALSE), collapse = "\n")
  np_c <- paste(readLines(test_path("..", "..", "src", "np.c"),
                          warn = FALSE), collapse = "\n")

  expect_match(
    header,
    "if\\(fit > 0\\.0 \\|\\| ISNAN\\(fit\\)\\)\\n    return log\\(fit\\);",
    perl = TRUE
  )
  expect_identical(fitted_subnormal_count_fixed(
    jksum, "np_fitted_log_likelihood_contribution("
  ), 5L)
  expect_identical(fitted_subnormal_count_fixed(
    np_c, "np_fitted_log_likelihood_contribution("
  ), 1L)
  expect_false(grepl("*log_likelihood += (pdf[i] < DBL_MIN)",
                     jksum, fixed = TRUE))
  expect_false(grepl("const double val = (pdf[j] < DBL_MIN)",
                     np_c, fixed = TRUE))
})

test_that("represented fitted-density subnormals contribute literal logs", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  density <- getFromNamespace("npudens", package)
  conditional_density <- getFromNamespace("npcdens", package)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  fit <- fitted_subnormal_run_local(package, density(
    bws = 1,
    tdat = data.frame(x = seq(-0.01, 0.01, length.out = 41L)),
    edat = data.frame(x = 38),
    bwtype = "fixed",
    bwscaling = FALSE,
    ckertype = "gaussian",
    ckerorder = 2L
  ))

  x <- data.frame(x = factor(rep(c("a", "b"), length.out = 41L)))
  conditional_fit <- fitted_subnormal_run_local(package, conditional_density(
    bws = c(1, 0.2),
    txdat = x,
    tydat = data.frame(y = seq(-0.02, 0.02, length.out = 41L)),
    exdat = data.frame(x = factor("a", levels = levels(x$x))),
    eydat = data.frame(y = 38),
    bwtype = "fixed",
    bwscaling = FALSE,
    regtype = "lc",
    cykertype = "gaussian",
    cykerorder = 2L,
    uxkertype = "aitchisonaitken"
  ))

  expect_true(fit$dens > 0 && fit$dens < .Machine$double.xmin)
  expect_true(conditional_fit$condens > 0 &&
              conditional_fit$condens < .Machine$double.xmin)
  expect_identical(fit$log_likelihood, log(fit$dens))
  expect_identical(conditional_fit$log_likelihood,
                   log(conditional_fit$condens))
})

test_that("fitted densities already underflowed to zero retain the placeholder", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  density <- getFromNamespace("npudens", package)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  fit <- fitted_subnormal_run_local(package, density(
    bws = 1,
    tdat = data.frame(x = seq(-0.01, 0.01, length.out = 41L)),
    edat = data.frame(x = 40),
    bwtype = "fixed",
    bwscaling = FALSE,
    ckertype = "gaussian",
    ckerorder = 2L
  ))

  expect_identical(fit$dens, 0)
  expect_identical(fit$log_likelihood, log(.Machine$double.xmin))
})
