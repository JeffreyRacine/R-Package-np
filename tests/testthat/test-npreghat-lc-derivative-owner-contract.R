npreghat_lc_derivative_owner_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_LC_DERIV_OWNER_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "set.seed(410)",
    "n <- 60L",
    "x <- rnorm(n)",
    "y <- x^2 + rnorm(n, sd = 0.2)",
    "train <- data.frame(x = x)",
    "eval <- data.frame(x = seq(min(x), max(x), length.out = 31L))",
    "run_case <- function(bt, bv, tol = 1e-10) {",
    "  bw <- npregbw(y ~ x, data = data.frame(train, y = y), regtype = 'lc', bwtype = bt, bandwidth.compute = FALSE, bws = bv)",
    "  Ht <- npreghat(bws = bw, txdat = train, s = 1L)",
    "  He <- npreghat(bws = bw, txdat = train, exdat = eval, s = 1L)",
    "  at <- npreghat(bws = bw, txdat = train, y = y, s = 1L, output = 'apply')",
    "  ae <- npreghat(bws = bw, txdat = train, exdat = eval, y = y, s = 1L, output = 'apply')",
    "  gt <- npreg(bws = bw, gradients = TRUE)",
    "  ge <- npreg(bws = bw, exdat = eval, gradients = TRUE)",
    "  stopifnot(isTRUE(all.equal(drop(Ht %*% y), at, tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(drop(Ht %*% y), gt$grad[, 1L], tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(drop(He %*% y), ae, tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(drop(He %*% y), ge$grad[, 1L], tolerance = tol)))",
    "}",
    "run_case('fixed', 0.45)",
    "run_case('generalized_nn', 9)",
    "run_case('adaptive_nn', 9)",
    "set.seed(20260630)",
    "n <- 36L",
    "idx <- seq_len(n)",
    "tx <- data.frame(",
    "  x1 = seq(0.08, 0.92, length.out = n),",
    "  x2 = sin(seq(0.2, 2.8, length.out = n)) + idx / 200,",
    "  u = factor(rep(c('a', 'b', 'c'), length.out = n), levels = c('a', 'b', 'c')),",
    "  o = ordered(rep(c('low', 'mid', 'high'), length.out = n), levels = c('low', 'mid', 'high'))",
    ")",
    "ex <- tx[seq_len(14L), , drop = FALSE]",
    "ex$x1 <- pmin(0.98, pmax(0.02, ex$x1 + 0.015))",
    "ex$x2 <- pmin(0.98, pmax(0.02, ex$x2 + 0.015))",
    "y <- 1.0 + sin(2 * pi * tx$x1) + 0.25 * tx$x2^2 + 0.12 * as.integer(tx$u) - 0.08 * as.integer(tx$o)",
    "Y <- cbind(y = y, y2 = y^2)",
    "s <- c(1L, 0L)",
    "run_mixed_case <- function(bt, bv, eval_data = NULL, tol = 1e-8) {",
    "  bw <- npregbw(xdat = tx, ydat = y, bws = bv, bandwidth.compute = FALSE, bwtype = bt, regtype = 'lc')",
    "  base_args <- list(bws = bw, txdat = tx, s = s)",
    "  if (!is.null(eval_data)) base_args$exdat <- eval_data",
    "  H <- do.call(npreghat, c(base_args, list(output = 'matrix')))",
    "  A1 <- do.call(npreghat, c(base_args, list(y = y, output = 'apply')))",
    "  A2 <- do.call(npreghat, c(base_args, list(y = Y, output = 'apply')))",
    "  C <- do.call(npreghat, c(base_args, list(y = y, output = 'constraint')))",
    "  Hmat <- matrix(as.numeric(H), nrow = nrow(H), ncol = ncol(H))",
    "  stopifnot(is.matrix(H), all(is.finite(as.matrix(H))))",
    "  stopifnot(ncol(H) == nrow(tx))",
    "  stopifnot(nrow(H) == if (is.null(eval_data)) nrow(tx) else nrow(eval_data))",
    "  stopifnot(isTRUE(all.equal(as.vector(A1), as.vector(Hmat %*% y), tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(unname(as.matrix(A2)), unname(Hmat %*% Y), tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(unname(as.matrix(C)), unname(t(Hmat) * as.vector(y)), tolerance = tol)))",
    "}",
    "for (eval_data in list(NULL, ex)) {",
    "  run_mixed_case('fixed', c(0.55, 0.55, 0.45, 0.45), eval_data = eval_data)",
    "  run_mixed_case('generalized_nn', c(7, 7, 0.45, 0.45), eval_data = eval_data)",
    "  run_mixed_case('adaptive_nn', c(7, 7, 0.45, 0.45), eval_data = eval_data)",
    "}",
    "bw <- npregbw(xdat = tx[, c('x1', 'x2')], ydat = y, bws = c(0.55, 0.55), bandwidth.compute = FALSE, bwtype = 'fixed', regtype = 'lc')",
    "err <- try(npreghat(bws = bw, txdat = tx[, c('x1', 'x2')], s = c(2L, 0L)), silent = TRUE)",
    "stopifnot(inherits(err, 'try-error'))",
    "stopifnot(grepl(\"requested derivative order in 's' exceeds local polynomial degree\", as.character(err), fixed = TRUE))",
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

test_that("npRmpi lc derivative owner matches npreg and apply on common cells", {
  npreghat_lc_derivative_owner_case()
})
