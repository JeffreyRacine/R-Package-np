local_npRmpi_semiparam_formula_env <- function(extra = character()) {
  pkg.root <- tryCatch(
    normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
    error = function(e) ""
  )
  if (!nzchar(pkg.root))
    return(NULL)

  src.copy <- tempfile("npRmpi-semiparam-formula-src-")
  dir.create(src.copy, recursive = TRUE, showWarnings = FALSE)
  pkg.copy <- file.path(src.copy, "np-npRmpi")
  if (!isTRUE(file.copy(pkg.root, src.copy, recursive = TRUE)))
    return(NULL)

  vignette.index <- file.path(pkg.copy, "build", "vignette.rds")
  if (file.exists(vignette.index))
    unlink(vignette.index, force = TRUE)

  build.dir <- tempfile("npRmpi-semiparam-formula-build-")
  dir.create(build.dir, recursive = TRUE, showWarnings = FALSE)

  old.wd <- getwd()
  on.exit(setwd(old.wd), add = TRUE)
  setwd(build.dir)

  build.out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "build", "--no-build-vignettes", "--no-manual", pkg.copy),
    stdout = TRUE,
    stderr = TRUE
  ))
  build.status <- attr(build.out, "status")
  if (is.null(build.status))
    build.status <- 0L
  if (build.status != 0L)
    return(NULL)

  tarballs <- Sys.glob(file.path(build.dir, "npRmpi_*.tar.gz"))
  if (!length(tarballs))
    return(NULL)

  lib.path <- tempfile("npRmpi-semiparam-formula-lib-")
  dir.create(lib.path, recursive = TRUE, showWarnings = FALSE)

  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-test-load", "-l", lib.path, tail(tarballs, 1L)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L

  if (status != 0L) {
    unlink(lib.path, recursive = TRUE, force = TRUE)
    return(NULL)
  }

  c(
    sprintf("R_LIBS=%s", paste(c(lib.path, .libPaths()), collapse = .Platform$path.sep)),
    "NP_RMPI_NO_REUSE_SLAVES=1",
    extra
  )
}

test_that("npRmpi semiparametric formula fits preserve response-name metadata", {
  skip_on_cran()
  env <- local_npRmpi_semiparam_formula_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(123)",
      "n <- 40L",
      "dat_pl <- data.frame(x = rnorm(n), z = sort(runif(n)))",
      "dat_pl$y <- 1 + dat_pl$x + sin(2 * pi * dat_pl$z) + rnorm(n, sd = 0.08)",
      "fit_pl <- npplreg(y ~ x | z, data = dat_pl, nmulti = 1, regtype = 'lp', degree = 1, bwtype = 'adaptive_nn')",
      "stopifnot(identical(fit_pl$bws$ynames, 'y'))",
      "dat_sc <- data.frame(x = runif(n), z = sort(runif(n)))",
      "dat_sc$y <- (1 + dat_sc$z^2) * dat_sc$x + rnorm(n, sd = 0.08)",
      "fit_sc <- npscoef(y ~ x | z, data = dat_sc, nmulti = 1, regtype = 'lp', degree = 1, bwtype = 'adaptive_nn', errors = FALSE, betas = FALSE)",
      "stopifnot(identical(fit_sc$bws$ynames, 'y'))",
      "dat_fn <- data.frame(z1 = factor(rbinom(n, 1L, 0.5)), z2 = factor(rbinom(n, 1L, 0.5)), z3 = factor(rbinom(n, 1L, 0.5)))",
      "dat_fn$x1 <- rnorm(n)",
      "dat_fn$x2 <- rnorm(n)",
      "dat_fn$y <- 0.75 * dat_fn$x1 - 0.5 * dat_fn$x2 + 0.6 * as.integer(dat_fn$z1) - 0.3 * as.integer(dat_fn$z2) + rnorm(n, sd = 0.25)",
      "xdat_fn <- dat_fn[, c('x1', 'x2')]",
      "zdat_fn <- dat_fn[, c('z1', 'z2', 'z3')]",
      "fit_sc_formula <- npscoef(y ~ x1 + x2 | z1 + z2 + z3, data = dat_fn, bws = rep(0.2, 3L), errors = FALSE, betas = FALSE)",
      "bw_sc_object <- npscoefbw(xdat = xdat_fn, ydat = dat_fn$y, zdat = zdat_fn, bws = rep(0.2, 3L), bandwidth.compute = FALSE)",
      "fit_sc_object <- npscoef(bws = bw_sc_object, txdat = xdat_fn, tydat = dat_fn$y, tzdat = zdat_fn, errors = FALSE, betas = FALSE)",
      "stopifnot(max(abs(as.numeric(fitted(fit_sc_formula)) - as.numeric(fitted(fit_sc_object)))) <= 1e-10)",
      "pl_bws_fn <- matrix(0.2, nrow = 3L, ncol = 3L)",
      "fit_pl_formula <- npplreg(y ~ x1 + x2 | z1 + z2 + z3, data = dat_fn, bws = pl_bws_fn)",
      "bw_pl_object <- npplregbw(xdat = xdat_fn, ydat = dat_fn$y, zdat = zdat_fn, bws = pl_bws_fn, bandwidth.compute = FALSE)",
      "fit_pl_object <- npplreg(bws = bw_pl_object, txdat = xdat_fn, tydat = dat_fn$y, tzdat = zdat_fn)",
      "stopifnot(max(abs(as.numeric(fitted(fit_pl_formula)) - as.numeric(fitted(fit_pl_object)))) <= 1e-10)",
      "dat_si <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))",
      "idx <- dat_si$x1 + 0.5 * dat_si$x2",
      "dat_si$y <- sin(idx) + 0.25 * idx^2 + rnorm(n, sd = 0.05)",
      "fit_si <- npindex(y ~ x1 + x2, data = dat_si, method = 'ichimura', nmulti = 1, regtype = 'lp', degree = 1, bwtype = 'adaptive_nn')",
      "stopifnot(identical(fit_si$bws$ynames, 'y'))",
      "tf <- tempfile(fileext = '.png')",
      "grDevices::png(tf, width = 800, height = 600)",
      "plot(fit_pl)",
      "plot(fit_sc)",
      "plot(fit_si)",
      "grDevices::dev.off()",
      "unlink(tf)",
      "cat('NP_RMPI_SEMIPARAM_FORMULA_YNAME_OK\\n')"
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NP_RMPI_SEMIPARAM_FORMULA_YNAME_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
