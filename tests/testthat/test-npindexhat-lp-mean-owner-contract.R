npindexhat_lp_mean_owner_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPINDEXHAT_LP_MEAN_OWNER_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(20260315)",
    "n <- 90L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- rnorm(n)",
    "y <- sin(2 * (0.7 * x1 - 0.3 * x2)) + 0.35 * x1 * x2 + rnorm(n, sd = 0.04)",
    "tx <- data.frame(x1 = x1, x2 = x2)",
    "ex <- data.frame(x1 = seq(min(x1) * 0.9, max(x1) * 0.9, length.out = 30L), x2 = seq(quantile(x2, 0.2), quantile(x2, 0.8), length.out = 30L))",
    "check_case <- function(label, bw) {",
    "  fit.mean <- npindex(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = FALSE, errors = FALSE)",
    "  apply.mean <- npindexhat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = 0L)",
    "  matrix.mean <- npindexhat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = 0L)",
    "  stopifnot(isTRUE(all.equal(as.vector(apply.mean), as.vector(fit.mean$mean), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(matrix.mean %*% y), as.vector(fit.mean$mean), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(apply.mean), as.vector(matrix.mean %*% y), tolerance = 1e-10)))",
    "}",
    "bw <- npindexbw(xdat = tx, ydat = y, regtype = 'lp', degree = 1L, basis = 'glp', bernstein.basis = FALSE, bwtype = 'generalized_nn', bandwidth.compute = FALSE, bws = c(1, 1, 9L))",
    "check_case('lp1 generalized canonical', bw)",
    "bw <- npindexbw(xdat = tx, ydat = y, regtype = 'lp', degree = 1L, basis = 'tensor', bernstein.basis = FALSE, bwtype = 'generalized_nn', bandwidth.compute = FALSE, bws = c(1, 1, 9L))",
    "check_case('lp1 generalized tensor legacy', bw)",
    "bw <- npindexbw(xdat = tx, ydat = y, regtype = 'lp', degree = 1L, basis = 'glp', bernstein.basis = TRUE, bwtype = 'generalized_nn', bandwidth.compute = FALSE, bws = c(1, 1, 9L))",
    "check_case('lp1 generalized bernstein legacy', bw)",
    "bw <- npindexbw(xdat = tx, ydat = y, regtype = 'lp', degree = 1L, basis = 'tensor', bernstein.basis = FALSE, bwtype = 'adaptive_nn', bandwidth.compute = FALSE, bws = c(1, 1, 9L))",
    "check_case('lp1 adaptive tensor', bw)",
    "bw <- npindexbw(xdat = tx, ydat = y, regtype = 'lp', degree = 2L, basis = 'tensor', bernstein.basis = FALSE, bwtype = 'generalized_nn', bandwidth.compute = FALSE, bws = c(1, 1, 11L))",
    "check_case('lp2 generalized tensor', bw)",
    "bw <- npindexbw(xdat = tx, ydat = y, regtype = 'lp', degree = 2L, basis = 'glp', bernstein.basis = TRUE, bwtype = 'generalized_nn', bandwidth.compute = FALSE, bws = c(1, 1, 11L))",
    "check_case('lp2 generalized bernstein', bw)",
    "bw <- npindexbw(xdat = tx, ydat = y, regtype = 'lp', degree = 2L, basis = 'tensor', bernstein.basis = FALSE, bwtype = 'adaptive_nn', bandwidth.compute = FALSE, bws = c(1, 1, 11L))",
    "check_case('lp2 adaptive tensor', bw)",
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

test_that("npRmpi npindexhat lp mean owner preserves fit and matrix/apply parity", {
  npindexhat_lp_mean_owner_case()
})
