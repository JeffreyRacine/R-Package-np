cvml_subnormal_run_local <- function(package, expression) {
  if (identical(package, "npRmpi"))
    getFromNamespace(".npRmpi_with_local_regression", package)(expression)
  else
    expression
}

cvml_subnormal_eval <- function(package, separation) {
  owner <- getFromNamespace("npudensbw", package)
  data <- data.frame(x = c(0, separation))
  bandwidth <- cvml_subnormal_run_local(package, owner(
    dat = data,
    bws = 1,
    bandwidth.compute = FALSE,
    bwmethod = "cv.ml",
    bwtype = "fixed",
    bwscaling = FALSE,
    ckertype = "gaussian",
    ckerorder = 2L
  ))
  cvml_subnormal_run_local(package, owner(
    dat = data,
    bws = bandwidth,
    bandwidth.compute = TRUE,
    eval.only = TRUE,
    invalid.penalty = "dbmax"
  ))
}

test_that("CVML compositors distinguish positive range from invalid sign", {
  header <- paste(
    readLines(test_path("..", "..", "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(header,
               "if\\(fit > 0\\.0\\)\\n    return -log\\(fit\\);",
               perl = TRUE)
  expect_match(header,
               "if\\(sign > 0 && log_fit > -INFINITY\\)\\n    return -log_fit;",
               perl = TRUE)
  expect_match(header,
               "if\\(fit < -DBL_MIN\\)\\n    return log\\(-fit\\) - 2\\.0\\*log\\(DBL_MIN\\);",
               perl = TRUE)
})

test_that("represented positive subnormal CVML values use their true logs", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  result <- cvml_subnormal_eval(package, 38)
  expected <- 2 * dnorm(38, log = TRUE)

  expect_identical(as.numeric(result$num.feval.guarded), 0)
  expect_lt(abs(result$fval - expected), 5e-11)
  expect_lt(result$fval, 2 * log(.Machine$double.xmin))
})

test_that("underflowed zero CVML values retain the guarded mapping", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  result <- cvml_subnormal_eval(package, 40)

  expect_identical(as.numeric(result$num.feval.guarded), 1)
  expect_identical(as.numeric(result$fval), 2 * log(.Machine$double.xmin))
})
