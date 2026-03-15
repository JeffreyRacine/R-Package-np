npindexhat_ll_owner_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPINDEXHAT_LL_OWNER_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(20260314)",
    "n <- 70L",
    "x1 <- runif(n)",
    "x2 <- runif(n)",
    "y <- x1 + 0.5 * x2 + rnorm(n, sd = 0.04)",
    "tx <- data.frame(x1 = x1, x2 = x2)",
    "ex <- tx[seq_len(18L), , drop = FALSE]",
    "for (bwtype in c('fixed', 'generalized_nn', 'adaptive_nn')) {",
    "  bw <- npindexbw(xdat = tx, ydat = y, regtype = 'll', bwtype = bwtype, bandwidth.compute = FALSE, bws = c(1, 1, if (identical(bwtype, 'fixed')) 0.18 else 9))",
    "  fit.mean <- npindex(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = FALSE)",
    "  fit.grad <- npindex(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = TRUE)",
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

test_that("npRmpi npindexhat ll owner preserves fit and matrix/apply parity", {
  npindexhat_ll_owner_case()
})
