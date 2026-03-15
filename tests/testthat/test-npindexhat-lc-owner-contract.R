npindexhat_lc_owner_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPINDEXHAT_LC_OWNER_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(20260315)",
    "n <- 90L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- rnorm(n)",
    "y <- sin(2 * (0.7 * x1 - 0.3 * x2)) + 0.25 * x1 + rnorm(n, sd = 0.04)",
    "tx <- data.frame(x1 = x1, x2 = x2)",
    "ex <- data.frame(x1 = seq(min(x1) * 0.9, max(x1) * 0.9, length.out = 30L), x2 = seq(quantile(x2, 0.2), quantile(x2, 0.8), length.out = 30L))",
    "for (bwtype in c('fixed', 'generalized_nn', 'adaptive_nn')) {",
    "  bw <- npindexbw(xdat = tx, ydat = y, regtype = 'lc', bwtype = bwtype, bandwidth.compute = FALSE, bws = c(1, 1, if (identical(bwtype, 'fixed')) 0.22 else 9L))",
    "  fit.mean <- npindex(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = FALSE, errors = FALSE)",
    "  fit.grad <- npindex(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = TRUE, errors = FALSE)",
    "  apply.mean <- npindexhat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = 0L)",
    "  matrix.mean <- npindexhat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = 0L)",
    "  apply.grad <- npindexhat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = 1L)",
    "  matrix.grad <- npindexhat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = 1L)",
    "  stopifnot(isTRUE(all.equal(as.vector(apply.mean), as.vector(fit.mean$mean), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(matrix.mean %*% y), as.vector(fit.mean$mean), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(apply.mean), as.vector(matrix.mean %*% y), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(as.vector(apply.grad), as.vector(fit.grad$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(matrix.grad %*% y), as.vector(fit.grad$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(apply.grad), as.vector(matrix.grad %*% y), tolerance = 1e-10)))",
    "}",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi npindexhat lc owner preserves fit and matrix/apply parity", {
  npindexhat_lc_owner_case()
})
