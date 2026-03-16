npscoef_semiparam_frozen_contract_case <- function() {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess contract")

  ok_tag <- "NPSCOEF_SEMIPARAM_FROZEN_CONTRACT_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(np))",
    "suppressPackageStartupMessages(library(npRmpi))",
    "run_case <- function() {",
    "  npRmpi.init(nslaves = 1, quiet = TRUE)",
    "  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    "  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
    "  xdat <- data.frame(x = c(0.05, 0.05, 0.25, 0.25, 0.60, 0.60, 0.90, 0.90))",
    "  zdat <- data.frame(z = c(0.10, 0.10, 0.35, 0.35, 0.65, 0.65, 0.95, 0.95))",
    "  y <- with(xdat, 1 + x * (1 + zdat$z))",
    "  exdat <- data.frame(x = c(0.10, 0.30, 0.70))",
    "  ezdat <- data.frame(z = c(0.12, 0.40, 0.85))",
    "  counts <- cbind(c(2, 0, 2, 0, 1, 1, 1, 1), c(0, 2, 0, 2, 1, 1, 1, 1), c(1, 1, 3, 0, 0, 2, 1, 0))",
    "  bw_np <- np::npscoefbw(xdat = xdat, zdat = zdat, ydat = y, bws = 3L, bwtype = 'adaptive_nn', bandwidth.compute = FALSE, regtype = 'll')",
    "  bw <- npscoefbw(xdat = xdat, zdat = zdat, ydat = y, bws = 3L, bwtype = 'adaptive_nn', bandwidth.compute = FALSE, regtype = 'll')",
    "  boot.fun <- getFromNamespace('.np_inid_boot_from_scoef', 'npRmpi')",
    "  exact.out <- boot.fun(txdat = xdat, ydat = y, tzdat = zdat, exdat = exdat, ezdat = ezdat, bws = bw, B = ncol(counts), counts = counts, mode = 'exact')",
    "  frozen.out <- boot.fun(txdat = xdat, ydat = y, tzdat = zdat, exdat = exdat, ezdat = ezdat, bws = bw, B = ncol(counts), counts = counts, mode = 'frozen')",
    "  exact.oracle <- vapply(seq_len(ncol(counts)), function(j) {",
    "    idx <- np:::.np_counts_to_indices(counts[, j])",
    "    as.vector(np::npscoef(bws = bw_np, txdat = xdat[idx, , drop = FALSE], tzdat = zdat[idx, , drop = FALSE], tydat = y[idx], exdat = exdat, ezdat = ezdat, iterate = FALSE, errors = FALSE)$mean)",
    "  }, numeric(nrow(exdat)))",
    "  exact.oracle <- t(exact.oracle)",
    "  t0.oracle <- as.vector(np::npscoef(bws = bw_np, txdat = xdat, tzdat = zdat, tydat = y, exdat = exdat, ezdat = ezdat, iterate = FALSE, errors = FALSE)$mean)",
    "  stopifnot(isTRUE(all.equal(exact.out$t, exact.oracle, tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(as.vector(exact.out$t0), t0.oracle, tolerance = 1e-10)))",
    "  if (isTRUE(all.equal(exact.out$t, frozen.out$t, tolerance = 1e-10, check.attributes = FALSE))) stop('frozen matched exact unexpectedly')",
    "  expected.plot <- plot(bw, xdat = xdat, ydat = y, zdat = zdat, neval = 6L, coef = FALSE, plot.behavior = 'data')",
    "  if (!is.null(expected.plot$r1)) expected.plot <- expected.plot$r1",
    "  ns <- asNamespace('npRmpi')",
    "  orig <- getFromNamespace('.np_inid_boot_from_scoef', 'npRmpi')",
    "  orig.compute <- getFromNamespace('compute.bootstrap.errors.scbandwidth', 'npRmpi')",
    "  modes <- character()",
    "  seen.ex <- list()",
    "  seen.ez <- list()",
    "  seen.t0 <- list()",
    "  assignInNamespace('.np_inid_boot_from_scoef', function(..., mode = c('exact', 'frozen')) {",
    "    mode <- match.arg(mode)",
    "    dots <- list(...)",
    "    modes <<- c(modes, mode)",
    "    seen.ex[[length(seen.ex) + 1L]] <<- dots$exdat",
    "    seen.ez[[length(seen.ez) + 1L]] <<- dots$ezdat",
    "    orig(..., mode = mode)",
    "  }, ns = 'npRmpi')",
    "  assignInNamespace('compute.bootstrap.errors.scbandwidth', function(..., t0.override = NULL) {",
    "    seen.t0[[length(seen.t0) + 1L]] <<- t0.override",
    "    orig.compute(..., t0.override = t0.override)",
    "  }, ns = 'npRmpi')",
    "  on.exit(assignInNamespace('.np_inid_boot_from_scoef', orig, ns = 'npRmpi'), add = TRUE)",
    "  on.exit(assignInNamespace('compute.bootstrap.errors.scbandwidth', orig.compute, ns = 'npRmpi'), add = TRUE)",
    "  capture.output(plot(bw, xdat = xdat, ydat = y, zdat = zdat, neval = 6L, coef = FALSE, plot.behavior = 'data', plot.errors.method = 'bootstrap', plot.errors.boot.method = 'inid', plot.errors.boot.nonfixed = 'frozen', plot.errors.boot.num = 41L, plot.errors.type = 'pointwise'))",
    "  stopifnot(length(modes) >= 1L)",
    "  stopifnot(all(modes == 'frozen'))",
    "  stopifnot(isTRUE(all.equal(seen.ex[[1L]], expected.plot$eval$exdat, tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(seen.ez[[1L]], expected.plot$eval$ezdat, tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(as.vector(seen.t0[[1L]]), as.vector(expected.plot$mean), tolerance = 1e-10)))",
    "}",
    "run_case()",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi npscoef exact/frozen split is correct and frozen is forwarded", {
  npscoef_semiparam_frozen_contract_case()
})

npscoef_fixed_semiparam_frozen_contract_case <- function() {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess contract")

  ok_tag <- "NPSCOEF_FIXED_SEMIPARAM_FROZEN_CONTRACT_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "run_case <- function() {",
    "  npRmpi.init(nslaves = 1, quiet = TRUE)",
    "  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
    "  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
    "  set.seed(42)",
    "  n <- 30L",
    "  xdat <- data.frame(x = runif(n, -1, 1))",
    "  zdat <- data.frame(z = rnorm(n))",
    "  y <- with(xdat, x^2 + rnorm(n, sd = 0.1))",
    "  exdat <- data.frame(x = seq(-0.9, 0.9, length.out = 9L))",
    "  ezdat <- data.frame(z = seq(-1.0, 1.0, length.out = 9L))",
    "  counts <- rmultinom(n = 5L, size = n, prob = rep.int(1 / n, n))",
    "  bw <- npscoefbw(xdat = xdat, zdat = zdat, ydat = y, bws = c(0.6), bwtype = 'fixed', bandwidth.compute = FALSE, regtype = 'll')",
    "  boot.fun <- getFromNamespace('.np_inid_boot_from_scoef', 'npRmpi')",
    "  exact.out <- boot.fun(txdat = xdat, ydat = y, tzdat = zdat, exdat = exdat, ezdat = ezdat, bws = bw, B = ncol(counts), counts = counts, mode = 'exact')",
    "  frozen.out <- boot.fun(txdat = xdat, ydat = y, tzdat = zdat, exdat = exdat, ezdat = ezdat, bws = bw, B = ncol(counts), counts = counts, mode = 'frozen')",
    "  stopifnot(isTRUE(all.equal(exact.out$t0, frozen.out$t0, tolerance = 1e-12)))",
    "  stopifnot(isTRUE(all.equal(exact.out$t, frozen.out$t, tolerance = 1e-12)))",
    "}",
    "run_case()",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi npscoef fixed exact matches frozen", {
  npscoef_fixed_semiparam_frozen_contract_case()
})
