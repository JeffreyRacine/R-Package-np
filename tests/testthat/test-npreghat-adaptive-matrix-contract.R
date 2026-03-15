npreghat_adaptive_matrix_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_ADAPTIVE_MATRIX_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(42)",
    "n <- 80L",
    "x <- runif(n)",
    "y <- x + rnorm(n)",
    "tx <- data.frame(x = x)",
    "ex <- data.frame(x = seq(min(x), max(x), length.out = 25L))",
    "bw <- npregbw(y ~ x, regtype = 'll', bwtype = 'adaptive_nn')",
    "mean.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply')",
    "mean.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix')",
    "grad.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = 1L)",
    "grad.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = 1L)",
    "np.mean <- npreg(txdat = tx, tydat = y, bws = bw, exdat = ex)$mean",
    "np.grad <- npreg(txdat = tx, tydat = y, bws = bw, exdat = ex, gradients = TRUE)$grad[, 1L]",
    "stopifnot(isTRUE(all.equal(drop(mean.matrix %*% y), mean.apply, tolerance = 1e-12)))",
    "stopifnot(isTRUE(all.equal(drop(mean.matrix %*% y), np.mean, tolerance = 1e-12)))",
    "stopifnot(isTRUE(all.equal(drop(grad.matrix %*% y), grad.apply, tolerance = 1e-12)))",
    "stopifnot(isTRUE(all.equal(drop(grad.matrix %*% y), np.grad, tolerance = 1e-12)))",
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

test_that("npRmpi adaptive ll matrix hat matches apply and npreg for mean and gradient", {
  npreghat_adaptive_matrix_case()
})
