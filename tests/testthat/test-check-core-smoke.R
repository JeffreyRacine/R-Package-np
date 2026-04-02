check_core_smoke_env <- function() {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess core smoke")
  env
}

run_check_core_smoke <- function(lines, timeout = 45L) {
  npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      lines
    ),
    timeout = timeout,
    env = check_core_smoke_env()
  )
}

expect_check_core_smoke_marker <- function(res, marker) {
  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(marker, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npudens core smoke stays alive", {
  res <- run_check_core_smoke(
    c(
      "set.seed(1)",
      "x <- data.frame(x = seq(0.1, 0.9, length.out = 24L))",
      "bw <- npudensbw(dat = x, bws = 0.2, bandwidth.compute = FALSE)",
      "fit <- npudens(bws = bw)",
      "stopifnot(inherits(fit, 'npdensity'))",
      "pred <- predict(fit)",
      "stopifnot(length(pred) == nrow(x))",
      "stopifnot(all(is.finite(pred)))",
      "cat('CHECK_CORE_SMOKE_NPUDENS_OK\\n')"
    )
  )

  expect_check_core_smoke_marker(res, "CHECK_CORE_SMOKE_NPUDENS_OK")
})

test_that("npudist core smoke stays alive", {
  res <- run_check_core_smoke(
    c(
      "set.seed(2)",
      "x <- data.frame(x = seq(0.1, 0.9, length.out = 24L))",
      "bw <- npudistbw(dat = x, bws = 0.2, bandwidth.compute = FALSE)",
      "fit <- npudist(bws = bw)",
      "stopifnot(inherits(fit, 'npdistribution'))",
      "pred <- predict(fit)",
      "stopifnot(length(pred) == nrow(x))",
      "stopifnot(all(pred >= 0 & pred <= 1))",
      "cat('CHECK_CORE_SMOKE_NPUDIST_OK\\n')"
    )
  )

  expect_check_core_smoke_marker(res, "CHECK_CORE_SMOKE_NPUDIST_OK")
})

test_that("npreg core smoke stays alive", {
  res <- run_check_core_smoke(
    c(
      "set.seed(3)",
      "x <- data.frame(x = seq(0.1, 1.0, length.out = 28L))",
      "y <- sin(2 * pi * x$x)",
      "bw <- npregbw(xdat = x, ydat = y, bws = 0.25, bandwidth.compute = FALSE)",
      "fit <- npreg(bws = bw)",
      "stopifnot(inherits(fit, 'npregression'))",
      "pred <- predict(fit)",
      "stopifnot(length(pred) == nrow(x))",
      "stopifnot(all(is.finite(pred)))",
      "cat('CHECK_CORE_SMOKE_NPREG_OK\\n')"
    )
  )

  expect_check_core_smoke_marker(res, "CHECK_CORE_SMOKE_NPREG_OK")
})

test_that("npcdens core smoke stays alive", {
  res <- run_check_core_smoke(
    c(
      "set.seed(4)",
      "x <- data.frame(x = seq(0.1, 1.0, length.out = 24L))",
      "y <- data.frame(y = x$x^2)",
      "bw <- npcdensbw(xdat = x, ydat = y, bws = c(0.25, 0.25), bandwidth.compute = FALSE)",
      "fit <- npcdens(bws = bw)",
      "stopifnot(inherits(fit, 'condensity'))",
      "pred <- predict(fit)",
      "stopifnot(length(pred) == nrow(x))",
      "stopifnot(all(is.finite(pred)))",
      "cat('CHECK_CORE_SMOKE_NPCDENS_OK\\n')"
    )
  )

  expect_check_core_smoke_marker(res, "CHECK_CORE_SMOKE_NPCDENS_OK")
})

test_that("npcdist core smoke stays alive", {
  res <- run_check_core_smoke(
    c(
      "set.seed(5)",
      "x <- data.frame(x = seq(0.1, 1.0, length.out = 24L))",
      "y <- data.frame(y = x$x^2)",
      "bw <- npcdistbw(xdat = x, ydat = y, bws = c(0.25, 0.25), bandwidth.compute = FALSE)",
      "fit <- npcdist(bws = bw)",
      "stopifnot(inherits(fit, 'condistribution'))",
      "pred <- predict(fit)",
      "stopifnot(length(pred) == nrow(x))",
      "stopifnot(all(pred >= 0 & pred <= 1))",
      "cat('CHECK_CORE_SMOKE_NPCDIST_OK\\n')"
    )
  )

  expect_check_core_smoke_marker(res, "CHECK_CORE_SMOKE_NPCDIST_OK")
})

test_that("npplreg core smoke stays alive", {
  res <- run_check_core_smoke(
    c(
      "set.seed(6)",
      "n <- 24L",
      "z <- data.frame(z = seq(0.1, 1.0, length.out = n))",
      "x <- data.frame(x = seq(0.2, 1.1, length.out = n))",
      "y <- z$z^2 + 2 * x$x",
      "bw <- npplregbw(xdat = x, zdat = z, ydat = y, bws = matrix(c(0.25, 0.25), nrow = 2L), bandwidth.compute = FALSE)",
      "fit <- npplreg(bws = bw)",
      "stopifnot(inherits(fit, 'plregression'))",
      "pred <- predict(fit)",
      "stopifnot(length(pred) == n)",
      "stopifnot(all(is.finite(pred)))",
      "cat('CHECK_CORE_SMOKE_NPPLREG_OK\\n')"
    )
  )

  expect_check_core_smoke_marker(res, "CHECK_CORE_SMOKE_NPPLREG_OK")
})

test_that("npindex core smoke stays alive", {
  res <- run_check_core_smoke(
    c(
      "set.seed(7)",
      "n <- 24L",
      "x <- data.frame(x1 = seq(0.1, 1.0, length.out = n), x2 = seq(1.0, 0.1, length.out = n))",
      "y <- x$x1 - x$x2",
      "bw <- npindexbw(xdat = x, ydat = y, method = 'ichimura', bws = c(1, 0.25, 0.25), bandwidth.compute = FALSE)",
      "fit <- npindex(bws = bw)",
      "stopifnot(inherits(fit, 'singleindex'))",
      "pred <- predict(fit)",
      "stopifnot(length(pred) == n)",
      "stopifnot(all(is.finite(pred)))",
      "cat('CHECK_CORE_SMOKE_NPINDEX_OK\\n')"
    )
  )

  expect_check_core_smoke_marker(res, "CHECK_CORE_SMOKE_NPINDEX_OK")
})

test_that("npscoef core smoke stays alive", {
  res <- run_check_core_smoke(
    c(
      "set.seed(8)",
      "n <- 24L",
      "x <- data.frame(x = seq(0.1, 1.0, length.out = n))",
      "z <- data.frame(z = seq(1.0, 0.1, length.out = n))",
      "y <- x$x * z$z",
      "bw <- npscoefbw(xdat = x, zdat = z, ydat = y, bws = 0.25, bandwidth.compute = FALSE)",
      "fit <- npscoef(bws = bw, iterate = FALSE, errors = FALSE)",
      "stopifnot(inherits(fit, 'smoothcoefficient'))",
      "pred <- predict(fit)",
      "stopifnot(length(pred) == n)",
      "stopifnot(all(is.finite(pred)))",
      "cat('CHECK_CORE_SMOKE_NPSCOEF_OK\\n')"
    )
  )

  expect_check_core_smoke_marker(res, "CHECK_CORE_SMOKE_NPSCOEF_OK")
})
