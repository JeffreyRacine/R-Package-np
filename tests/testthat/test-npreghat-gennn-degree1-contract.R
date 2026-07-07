test_that("generalized-nn ll and canonical lp degree-1 share the same public/apply regression route", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess npreghat smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(npRmpi.autodispatch = TRUE, np.messages = FALSE)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(20260308)",
      "n <- 90L",
      "x <- sort(runif(n))",
      "y <- sin(2 * pi * x) + 0.25 * x + rnorm(n, sd = 0.04)",
      "tx <- data.frame(x = x)",
      "ex <- data.frame(x = seq(0.05, 0.95, length.out = 25L))",
      "bw.ll <- npregbw(",
      "  xdat = tx,",
      "  ydat = y,",
      "  regtype = 'll',",
      "  bwtype = 'generalized_nn',",
      "  bws = 9L,",
      "  bandwidth.compute = FALSE",
      ")",
      "bw.lp <- npregbw(",
      "  xdat = tx,",
      "  ydat = y,",
      "  regtype = 'lp',",
      "  basis = 'glp',",
      "  degree = 1L,",
      "  bernstein.basis = FALSE,",
      "  bwtype = 'generalized_nn',",
      "  bws = 9L,",
      "  bandwidth.compute = FALSE",
      ")",
      "fit.ll <- npreg(bws = bw.ll, txdat = tx, tydat = y, exdat = ex, warn.glp.gradient = FALSE)",
      "fit.lp <- npreg(bws = bw.lp, txdat = tx, tydat = y, exdat = ex, warn.glp.gradient = FALSE)",
      "hat.apply.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex, y = y, output = 'apply')",
      "hat.apply.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex, y = y, output = 'apply')",
      "hat.matrix.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex)",
      "hat.matrix.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex)",
      "stopifnot(isTRUE(all.equal(as.vector(fit.ll$mean), as.vector(fit.lp$mean), tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(as.vector(hat.apply.ll), as.vector(hat.apply.lp), tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(as.vector(hat.apply.ll), as.vector(hat.matrix.ll %*% y), tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(as.vector(hat.apply.lp), as.vector(hat.matrix.lp %*% y), tolerance = 1e-10)))",
      "cat('NPREGHAT_GENNN_DEGREE1_SUBPROCESS_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_true(res$status %in% c(0L, 137L), info = info)
  expect_true(any(grepl("NPREGHAT_GENNN_DEGREE1_SUBPROCESS_OK",
                        res$output,
                        fixed = TRUE)),
              info = info)
})
