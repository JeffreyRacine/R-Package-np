nomad_shortcut_env <- function() {
  skip_on_cran()
  skip_if_not(
    tolower(Sys.getenv("NOT_CRAN", "")) %in% c("true", "1", "yes"),
    "extended local subprocess NOMAD shortcut contract"
  )
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")
  env
}

run_nomad_shortcut_subprocess <- function(lines, timeout = 120L) {
  npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      lines
    ),
    timeout = timeout,
    env = nomad_shortcut_env()
  )
}

expect_nomad_shortcut_marker <- function(res, marker) {
  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(marker, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

nomad_shortcut_timing_lines <- function(fit_expr = "fit") {
  c(
    sprintf("stopifnot(is.finite(%s$bws$nomad.time))", fit_expr),
    sprintf("stopifnot(is.finite(%s$bws$powell.time))", fit_expr),
    sprintf("stopifnot(isTRUE(all.equal(as.double(%s$bws$total.time), as.double(%s$bws$nomad.time + %s$bws$powell.time), tolerance = 1e-8)))", fit_expr, fit_expr, fit_expr),
    sprintf("if (!is.null(%s$optim.time) && length(%s$optim.time)) stopifnot(isTRUE(all.equal(as.double(%s$optim.time), as.double(%s$bws$total.time), tolerance = 1e-8)))", fit_expr, fit_expr, fit_expr, fit_expr),
    sprintf("if (!is.null(%s$total.time) && length(%s$total.time) && !is.null(%s$fit.time) && length(%s$fit.time) && !is.null(%s$optim.time) && length(%s$optim.time)) stopifnot(isTRUE(all.equal(as.double(%s$total.time), as.double(%s$optim.time + %s$fit.time), tolerance = 1e-8)))", fit_expr, fit_expr, fit_expr, fit_expr, fit_expr, fit_expr, fit_expr, fit_expr, fit_expr)
  )
}

test_that("nomad shortcut smoke covers npreg in subprocess session mode", {
  res <- run_nomad_shortcut_subprocess(
    lines = c(
      "set.seed(20260322)",
      "n <- 12L",
      "x <- runif(n, -1, 1)",
      "y <- x + 0.4 * x^2 + rnorm(n, sd = 0.18)",
      "dat <- data.frame(y = y, x = x)",
      "fit <- npreg(y ~ x, data = dat, nomad = TRUE, degree.max = 1L, nmulti = 1L)",
      "stopifnot(inherits(fit, 'npregression'))",
      "stopifnot(isTRUE(fit$bws$nomad.shortcut$enabled))",
      "stopifnot(identical(fit$bws$nomad.shortcut$preset, 'lp_nomad'))",
      nomad_shortcut_timing_lines(),
      "cat('NOMAD_SHORTCUT_NPREG_OK\\n')"
    )
  )

  expect_nomad_shortcut_marker(res, "NOMAD_SHORTCUT_NPREG_OK")
})

test_that("nomad shortcut smoke covers npcdens in subprocess session mode", {
  res <- run_nomad_shortcut_subprocess(
    lines = c(
      "set.seed(20260323)",
      "n <- 12L",
      "x <- runif(n, -1, 1)",
      "y <- x + 0.35 * x^2 + rnorm(n, sd = 0.18)",
      "dat <- data.frame(y = y, x = x)",
      "fit <- npcdens(y ~ x, data = dat, nomad = TRUE, degree.max = 1L, nmulti = 1L)",
      "stopifnot(inherits(fit, 'condensity'))",
      "stopifnot(isTRUE(fit$bws$nomad.shortcut$enabled))",
      "stopifnot(isTRUE(fit$bws$bernstein.basis))",
      "stopifnot(identical(fit$bws$nomad.shortcut$preset, 'lp_nomad'))",
      nomad_shortcut_timing_lines(),
      "cat('NOMAD_SHORTCUT_NPCDENS_OK\\n')"
    )
  )

  expect_nomad_shortcut_marker(res, "NOMAD_SHORTCUT_NPCDENS_OK")
})

test_that("nomad shortcut smoke covers npcdist in subprocess session mode", {
  res <- run_nomad_shortcut_subprocess(
    lines = c(
      "set.seed(20260324)",
      "n <- 12L",
      "x <- runif(n, -1, 1)",
      "y <- x + 0.25 * x^2 + rnorm(n, sd = 0.18)",
      "dat <- data.frame(y = y, x = x)",
      "fit <- npcdist(y ~ x, data = dat, nomad = TRUE, degree.max = 1L, nmulti = 1L, ngrid = 7L)",
      "stopifnot(inherits(fit, 'condistribution'))",
      "stopifnot(isTRUE(fit$bws$nomad.shortcut$enabled))",
      "stopifnot(isTRUE(fit$bws$bernstein.basis))",
      "stopifnot(identical(fit$bws$nomad.shortcut$preset, 'lp_nomad'))",
      nomad_shortcut_timing_lines(),
      "cat('NOMAD_SHORTCUT_NPCDIST_OK\\n')"
    )
  )

  expect_nomad_shortcut_marker(res, "NOMAD_SHORTCUT_NPCDIST_OK")
})

test_that("nomad shortcut smoke covers npplreg in subprocess session mode", {
  res <- run_nomad_shortcut_subprocess(
    lines = c(
      "set.seed(20260325)",
      "n <- 12L",
      "x <- runif(n)",
      "z <- runif(n, -1, 1)",
      "y <- 1 + x + sin(pi * z) + rnorm(n, sd = 0.18)",
      "dat <- data.frame(y = y, x = x, z = z)",
      "fit <- npplreg(y ~ x | z, data = dat, nomad = TRUE, degree.max = 1L, nmulti = 1L)",
      "stopifnot(inherits(fit, 'plregression'))",
      "stopifnot(isTRUE(fit$bws$nomad.shortcut$enabled))",
      "stopifnot(identical(fit$bws$nomad.shortcut$preset, 'lp_nomad'))",
      nomad_shortcut_timing_lines(),
      "cat('NOMAD_SHORTCUT_NPPLREG_OK\\n')"
    )
  )

  expect_nomad_shortcut_marker(res, "NOMAD_SHORTCUT_NPPLREG_OK")
})

test_that("nomad shortcut smoke covers npscoef in subprocess session mode", {
  res <- run_nomad_shortcut_subprocess(
    lines = c(
      "set.seed(20260326)",
      "n <- 12L",
      "x <- runif(n)",
      "z <- runif(n, -1, 1)",
      "y <- 1 + x + sin(pi * z) + rnorm(n, sd = 0.18)",
      "dat <- data.frame(y = y, x = x, z = z)",
      "fit <- npscoef(y ~ x | z, data = dat, nomad = TRUE, degree.max = 1L, nmulti = 1L)",
      "stopifnot(inherits(fit, 'smoothcoefficient'))",
      "stopifnot(isTRUE(fit$bws$nomad.shortcut$enabled))",
      "stopifnot(identical(fit$bws$nomad.shortcut$preset, 'lp_nomad'))",
      nomad_shortcut_timing_lines(),
      "cat('NOMAD_SHORTCUT_NPSCOEF_OK\\n')"
    )
  )

  expect_nomad_shortcut_marker(res, "NOMAD_SHORTCUT_NPSCOEF_OK")
})

test_that("nomad shortcut smoke covers npindex in subprocess session mode", {
  res <- run_nomad_shortcut_subprocess(
    lines = c(
      "set.seed(20260327)",
      "n <- 12L",
      "x1 <- runif(n, -1, 1)",
      "x2 <- runif(n, -1, 1)",
      "idx <- x1 + 0.75 * x2",
      "y <- sin(idx) + 0.25 * idx^2 + rnorm(n, sd = 0.08)",
      "dat <- data.frame(y = y, x1 = x1, x2 = x2)",
      "fit <- npindex(y ~ x1 + x2, data = dat, method = 'ichimura', nomad = TRUE, degree.max = 1L, nmulti = 1L)",
      "stopifnot(inherits(fit, 'singleindex'))",
      "stopifnot(isTRUE(fit$bws$nomad.shortcut$enabled))",
      "stopifnot(identical(fit$bws$nomad.shortcut$preset, 'lp_nomad'))",
      "stopifnot(identical(fit$bws$ynames, 'y'))",
      nomad_shortcut_timing_lines(),
      "cat('NOMAD_SHORTCUT_NPINDEX_OK\\n')"
    )
  )

  expect_nomad_shortcut_marker(res, "NOMAD_SHORTCUT_NPINDEX_OK")
})

test_that("nomad.nmulti fails fast outside active NOMAD degree search in subprocess session mode", {
  res <- run_nomad_shortcut_subprocess(
    lines = c(
      "x <- data.frame(x = seq(0, 1, length.out = 8))",
      "y <- x$x^2",
      "msg <- tryCatch({",
      "  npregbw(xdat = x, ydat = y, nomad.nmulti = 1L)",
      "  'NO_ERROR'",
      "}, error = conditionMessage)",
      "stopifnot(grepl('nomad.nmulti is only supported', msg, fixed = TRUE))",
      "cat('NOMAD_SHORTCUT_FAILFAST_OK\\n')"
    ),
    timeout = 60L
  )

  expect_nomad_shortcut_marker(res, "NOMAD_SHORTCUT_FAILFAST_OK")
})
