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

npreghat_adaptive_lc_matrix_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_ADAPTIVE_LC_MATRIX_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(20260315)",
    "n <- 64L",
    "x <- sort(runif(n))",
    "y <- sin(2 * pi * x) + 0.2 * x + rnorm(n, sd = 0.03)",
    "tx <- data.frame(x = x)",
    "ex <- data.frame(x = seq(0.05, 0.95, length.out = 21L))",
    "bw <- npregbw(xdat = tx, ydat = y, regtype = 'lc', bwtype = 'adaptive_nn', bandwidth.compute = FALSE, bws = 9L)",
    "fit.train <- npreg(bws = bw, txdat = tx, tydat = y)",
    "fit.eval <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex)",
    "hat.train.matrix <- npreghat(bws = bw, txdat = tx, output = 'matrix')",
    "hat.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = 'apply')",
    "hat.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix')",
    "hat.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply')",
    "stopifnot(isTRUE(all.equal(drop(hat.train.matrix %*% y), as.vector(fit.train$mean), tolerance = 1e-12)))",
    "stopifnot(isTRUE(all.equal(as.vector(hat.train.apply), as.vector(fit.train$mean), tolerance = 1e-12)))",
    "stopifnot(isTRUE(all.equal(drop(hat.train.matrix %*% y), as.vector(hat.train.apply), tolerance = 1e-12)))",
    "stopifnot(isTRUE(all.equal(drop(hat.eval.matrix %*% y), as.vector(fit.eval$mean), tolerance = 1e-12)))",
    "stopifnot(isTRUE(all.equal(as.vector(hat.eval.apply), as.vector(fit.eval$mean), tolerance = 1e-12)))",
    "stopifnot(isTRUE(all.equal(drop(hat.eval.matrix %*% y), as.vector(hat.eval.apply), tolerance = 1e-12)))",
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

test_that("npRmpi adaptive lc matrix hat matches fitted values on train and eval data", {
  npreghat_adaptive_lc_matrix_case()
})

npreghat_adaptive_lp_matrix_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_ADAPTIVE_LP_MATRIX_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(20260315)",
    "n <- 72L",
    "x <- sort(rnorm(n))",
    "y <- x^2 + 0.3 * x + rnorm(n, sd = 0.05)",
    "tx <- data.frame(x = x)",
    "ex <- data.frame(x = seq(min(x), max(x), length.out = 27L))",
    "for (degree in 1:3) {",
    "  bw <- npregbw(xdat = tx, ydat = y, regtype = 'lp', degree = degree, basis = 'glp', bernstein.basis = FALSE, bwtype = 'adaptive_nn', bandwidth.compute = FALSE, bws = 9L)",
    "  fit.train <- npreg(bws = bw, txdat = tx, tydat = y, gradients = TRUE)",
    "  fit.eval <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = TRUE)",
    "  mean.train.matrix <- npreghat(bws = bw, txdat = tx, output = 'matrix')",
    "  mean.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = 'apply')",
    "  mean.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix')",
    "  mean.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply')",
    "  grad.train.matrix <- npreghat(bws = bw, txdat = tx, output = 'matrix', s = 1L)",
    "  grad.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = 'apply', s = 1L)",
    "  grad.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = 1L)",
    "  grad.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = 1L)",
    "  stopifnot(isTRUE(all.equal(drop(mean.train.matrix %*% y), as.vector(fit.train$mean), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(as.vector(mean.train.apply), as.vector(fit.train$mean), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(drop(mean.train.matrix %*% y), as.vector(mean.train.apply), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(drop(mean.eval.matrix %*% y), as.vector(fit.eval$mean), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(as.vector(mean.eval.apply), as.vector(fit.eval$mean), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(drop(mean.eval.matrix %*% y), as.vector(mean.eval.apply), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(drop(grad.train.matrix %*% y), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(grad.train.apply), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(drop(grad.train.matrix %*% y), as.vector(grad.train.apply), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(drop(grad.eval.matrix %*% y), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(grad.eval.apply), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(drop(grad.eval.matrix %*% y), as.vector(grad.eval.apply), tolerance = 1e-8)))",
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

test_that("npRmpi adaptive lp matrix hat matches apply and npreg across degrees", {
  npreghat_adaptive_lp_matrix_case()
})

npreghat_nonadaptive_lcll_matrix_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_NONADAPTIVE_LCLL_MATRIX_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(20260316)",
    "n <- 96L",
    "x <- sort(rnorm(n))",
    "y <- x^2 + 0.15 * x + rnorm(n, sd = 0.05)",
    "tx <- data.frame(x = x)",
    "ex <- data.frame(x = seq(min(x), max(x), length.out = 31L))",
    "for (regtype in c('lc', 'll')) {",
    "  for (bwtype in c('fixed', 'generalized_nn')) {",
    "    bw <- npregbw(xdat = tx, ydat = y, regtype = regtype, bwtype = bwtype, bandwidth.compute = FALSE, bws = if (identical(bwtype, 'fixed')) 0.55 else 9L)",
    "    fit.train <- npreg(bws = bw, txdat = tx, tydat = y, gradients = identical(regtype, 'll'))",
    "    fit.eval <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = identical(regtype, 'll'))",
    "    mean.train.matrix <- npreghat(bws = bw, txdat = tx, output = 'matrix')",
    "    mean.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = 'apply')",
    "    mean.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix')",
    "    mean.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply')",
    "    stopifnot(isTRUE(all.equal(drop(mean.train.matrix %*% y), as.vector(fit.train$mean), tolerance = 1e-10)))",
    "    stopifnot(isTRUE(all.equal(as.vector(mean.train.apply), as.vector(fit.train$mean), tolerance = 1e-10)))",
    "    stopifnot(isTRUE(all.equal(drop(mean.eval.matrix %*% y), as.vector(fit.eval$mean), tolerance = 1e-10)))",
    "    stopifnot(isTRUE(all.equal(as.vector(mean.eval.apply), as.vector(fit.eval$mean), tolerance = 1e-10)))",
    "    if (identical(regtype, 'll')) {",
    "      grad.train.matrix <- npreghat(bws = bw, txdat = tx, output = 'matrix', s = 1L)",
    "      grad.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = 'apply', s = 1L)",
    "      grad.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = 1L)",
    "      grad.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = 1L)",
    "      stopifnot(isTRUE(all.equal(drop(grad.train.matrix %*% y), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)))",
    "      stopifnot(isTRUE(all.equal(as.vector(grad.train.apply), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)))",
    "      stopifnot(isTRUE(all.equal(drop(grad.eval.matrix %*% y), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)))",
    "      stopifnot(isTRUE(all.equal(as.vector(grad.eval.apply), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)))",
    "    }",
    "  }",
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

test_that("npRmpi fixed and generalized lc/ll matrix hats match apply and npreg", {
  npreghat_nonadaptive_lcll_matrix_case()
})

npreghat_fixed_lp_matrix_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_FIXED_LP_MATRIX_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "set.seed(20260317)",
    "n <- 96L",
    "x <- sort(rnorm(n))",
    "y <- x^3 - 0.4 * x + rnorm(n, sd = 0.05)",
    "tx <- data.frame(x = x)",
    "ex <- data.frame(x = seq(min(x), max(x), length.out = 31L))",
    "for (degree in 2:3) {",
    "  bw <- npregbw(xdat = tx, ydat = y, regtype = 'lp', degree = degree, basis = 'glp', bernstein.basis = FALSE, bwtype = 'fixed', bandwidth.compute = FALSE, bws = 0.55)",
    "  fit.train <- npreg(bws = bw, txdat = tx, tydat = y, gradients = TRUE)",
    "  fit.eval <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = TRUE)",
    "  mean.train.matrix <- npreghat(bws = bw, txdat = tx, output = 'matrix')",
    "  mean.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = 'apply')",
    "  mean.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix')",
    "  mean.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply')",
    "  grad.train.matrix <- npreghat(bws = bw, txdat = tx, output = 'matrix', s = 1L)",
    "  grad.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = 'apply', s = 1L)",
    "  grad.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = 1L)",
    "  grad.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = 'apply', s = 1L)",
    "  stopifnot(isTRUE(all.equal(drop(mean.train.matrix %*% y), as.vector(fit.train$mean), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(as.vector(mean.train.apply), as.vector(fit.train$mean), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(drop(mean.eval.matrix %*% y), as.vector(fit.eval$mean), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(as.vector(mean.eval.apply), as.vector(fit.eval$mean), tolerance = 1e-10)))",
    "  stopifnot(isTRUE(all.equal(drop(grad.train.matrix %*% y), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(grad.train.apply), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(drop(grad.eval.matrix %*% y), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)))",
    "  stopifnot(isTRUE(all.equal(as.vector(grad.eval.apply), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)))",
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

test_that("npRmpi fixed higher-degree lp matrix hats match apply and npreg", {
  npreghat_fixed_lp_matrix_case()
})
