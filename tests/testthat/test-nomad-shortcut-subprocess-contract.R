test_that("nomad shortcut smoke covers regression and conditional families in subprocess session mode", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(20260322)",
      "n <- 28L",
      "x <- runif(n, -1, 1)",
      "y <- x + 0.4 * x^2 + rnorm(n, sd = 0.18)",
      "dat <- data.frame(y = y, x = x)",
      "fit_reg <- npreg(y ~ x, data = dat, nomad = TRUE, nmulti = 1)",
      "stopifnot(isTRUE(fit_reg$bws$nomad.shortcut$enabled))",
      "stopifnot(identical(fit_reg$bws$nomad.shortcut$preset, 'lp_nomad'))",
      "bw_cd <- npcdensbw(y ~ x, data = dat, nomad = TRUE, nmulti = 1)",
      "stopifnot(isTRUE(bw_cd$nomad.shortcut$enabled))",
      "stopifnot(identical(bw_cd$nomad.shortcut$preset, 'lp_nomad'))",
      "fit_cd <- npcdens(bws = bw_cd)",
      "stopifnot(inherits(fit_cd, 'condensity'))",
      "bw_cdf <- npcdistbw(y ~ x, data = dat, nomad = TRUE, nmulti = 1, ngrid = 25L)",
      "stopifnot(isTRUE(bw_cdf$nomad.shortcut$enabled))",
      "stopifnot(identical(bw_cdf$nomad.shortcut$preset, 'lp_nomad'))",
      "fit_cdf <- npcdist(bws = bw_cdf)",
      "stopifnot(inherits(fit_cdf, 'condistribution'))",
      "cat('NOMAD_SHORTCUT_CORE_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NOMAD_SHORTCUT_CORE_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("nomad shortcut smoke covers semiparametric families in subprocess session mode", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(20260323)",
      "n <- 28L",
      "x <- runif(n)",
      "z <- runif(n, -1, 1)",
      "y <- 1 + x + sin(pi * z) + rnorm(n, sd = 0.18)",
      "dat_pl <- data.frame(y = y, x = x, z = z)",
      "fit_pl <- npplreg(y ~ x | z, data = dat_pl, nomad = TRUE, nmulti = 1)",
      "stopifnot(isTRUE(fit_pl$bws$nomad.shortcut$enabled))",
      "stopifnot(identical(fit_pl$bws$nomad.shortcut$preset, 'lp_nomad'))",
      "fit_sc <- npscoef(y ~ x | z, data = dat_pl, nomad = TRUE, nmulti = 1)",
      "stopifnot(isTRUE(fit_sc$bws$nomad.shortcut$enabled))",
      "stopifnot(identical(fit_sc$bws$nomad.shortcut$preset, 'lp_nomad'))",
      "dat_si <- data.frame(y = y, x = x)",
      "fit_si <- npindex(y ~ x, data = dat_si, nomad = TRUE, nmulti = 1)",
      "stopifnot(isTRUE(fit_si$bws$nomad.shortcut$enabled))",
      "stopifnot(identical(fit_si$bws$nomad.shortcut$preset, 'lp_nomad'))",
      "cat('NOMAD_SHORTCUT_SEMIPARAM_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NOMAD_SHORTCUT_SEMIPARAM_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
