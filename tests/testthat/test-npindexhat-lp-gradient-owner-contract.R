npindexhat_lp_gradient_owner_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPINDEXHAT_LP_GRAD_OWNER_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(20260316)",
    "n <- 90L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- rnorm(n)",
    "y <- sin(2 * (0.7 * x1 - 0.3 * x2)) + 0.35 * x1 * x2 + rnorm(n, sd = 0.04)",
    "tx <- data.frame(x1 = x1, x2 = x2)",
    "ex <- data.frame(x1 = seq(min(x1) * 0.9, max(x1) * 0.9, length.out = 30L), x2 = seq(quantile(x2, 0.2), quantile(x2, 0.8), length.out = 30L))",
    "cases <- list(",
    "  list(label = 'lp1 generalized canonical', bwtype = 'generalized_nn', degree = 1L, basis = 'glp', bern = FALSE, bws = c(1, 1, 9L)),",
    "  list(label = 'lp1 generalized tensor legacy', bwtype = 'generalized_nn', degree = 1L, basis = 'tensor', bern = FALSE, bws = c(1, 1, 9L)),",
    "  list(label = 'lp1 generalized bernstein legacy', bwtype = 'generalized_nn', degree = 1L, basis = 'glp', bern = TRUE, bws = c(1, 1, 9L)),",
    "  list(label = 'lp1 adaptive tensor', bwtype = 'adaptive_nn', degree = 1L, basis = 'tensor', bern = FALSE, bws = c(1, 1, 9L)),",
    "  list(label = 'lp2 generalized tensor', bwtype = 'generalized_nn', degree = 2L, basis = 'tensor', bern = FALSE, bws = c(1, 1, 11L)),",
    "  list(label = 'lp2 generalized bernstein', bwtype = 'generalized_nn', degree = 2L, basis = 'glp', bern = TRUE, bws = c(1, 1, 11L)),",
    "  list(label = 'lp2 adaptive tensor', bwtype = 'adaptive_nn', degree = 2L, basis = 'tensor', bern = FALSE, bws = c(1, 1, 11L))",
    ")",
    "for (cfg in cases) {",
    "  bw.args <- list(xdat = tx, ydat = y, regtype = 'lp', degree = cfg$degree, basis = cfg$basis, bernstein.basis = cfg$bern, bwtype = cfg$bwtype, bandwidth.compute = FALSE, bws = cfg$bws)",
    "  bw <- do.call(npindexbw, bw.args)",
    "  fit.grad <- npindex(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = TRUE, errors = FALSE)",
    "  apply.grad <- npindexhat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = 1L)",
    "  matrix.grad <- npindexhat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = 1L)",
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

test_that("npRmpi npindexhat lp gradient owner preserves fit and matrix/apply parity", {
  npindexhat_lp_gradient_owner_case()
})
