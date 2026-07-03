local_npreg_quiet <- function(expr) {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  force(expr)
}

expect_omit_rows <- function(x, expected) {
  if (!length(expected)) {
    expect_length(x, 1L)
    expect_true(is.na(x))
  } else {
    expect_identical(as.integer(x), as.integer(expected))
  }
}

test_that("npreg formula bandwidth honors explicit data override", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  local_npreg_quiet({
    d1 <- data.frame(x = seq(0, 1, length.out = 40))
    d1$y <- 1 + d1$x
    d2 <- data.frame(x = seq(0, 1, length.out = 12))
    d2$y <- 10 - 2 * d2$x

    bw <- npRmpi::npregbw(y ~ x, data = d1, bws = 0.35, bandwidth.compute = FALSE)
    d1$y <- 100 + d1$x

    fit <- npRmpi::npreg(bws = bw, data = d2)
    direct <- npRmpi::npreg(txdat = d2["x"], tydat = d2$y, bws = bw)

    expect_identical(fit$ntrain, nrow(d2))
    expect_equal(fit$mean, direct$mean, tolerance = 1e-12)
    expect_equal(fit$merr, direct$merr, tolerance = 1e-12)
  })
})

test_that("npreg stores training and evaluation omission metadata separately", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  local_npreg_quiet({
    tx <- data.frame(x = c(0, 0.2, NA, 0.5, 0.9, 1.1))
    ty <- c(1, 2, 3, 4, 5, NA)
    ex <- data.frame(x = c(0.1, NA, 0.7, 0.8))
    ey <- c(1, 2, NA, 4)
    bw <- npRmpi::npregbw(xdat = tx, ydat = ty, bws = 0.4, bandwidth.compute = FALSE)

    fit_eval <- npRmpi::npreg(bws = bw, txdat = tx, tydat = ty, exdat = ex, eydat = ey)
    expect_omit_rows(fit_eval$train.rows.omit, c(3L, 6L))
    expect_identical(fit_eval$train.nobs.omit, 2L)
    expect_omit_rows(fit_eval$eval.rows.omit, c(2L, 3L))
    expect_identical(fit_eval$eval.nobs.omit, 2L)
    expect_omit_rows(fit_eval$rows.omit, c(2L, 3L))

    fit_train <- npRmpi::npreg(bws = bw, txdat = tx, tydat = ty)
    expect_omit_rows(fit_train$train.rows.omit, c(3L, 6L))
    expect_identical(fit_train$train.nobs.omit, 2L)
    expect_omit_rows(fit_train$eval.rows.omit, integer(0))
    expect_equal(fit_train$eval.nobs.omit, 0)
    expect_omit_rows(fit_train$rows.omit, c(3L, 6L))
  })
})

test_that("npreg formula route reports original training and evaluation omissions", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  local_npreg_quiet({
    old.na <- getOption("na.action")
    on.exit(options(na.action = old.na), add = TRUE)
    options(na.action = "na.exclude")

    dat <- data.frame(
      y = c(1, NA, 3, 4, 5, 6),
      x = c(0, 0.2, 0.4, NA, 0.8, 1)
    )
    newdat <- data.frame(
      y = c(1, 2, NA, 4),
      x = c(0.1, NA, 0.7, 0.8)
    )

    bw <- npRmpi::npregbw(y ~ x, data = dat, bws = 0.4, bandwidth.compute = FALSE)
    fit <- npRmpi::npreg(bws = bw, data = dat, newdata = newdat, y.eval = TRUE)

    expect_omit_rows(fit$train.rows.omit, c(2L, 4L))
    expect_identical(fit$train.nobs.omit, 2L)
    expect_omit_rows(fit$eval.rows.omit, c(2L, 3L))
    expect_identical(fit$eval.nobs.omit, 2L)
    expect_omit_rows(fit$rows.omit, c(2L, 3L))
  })
})

test_that("npreg residuals reuse the training fit result for training-data fits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  local_npreg_quiet({
    x <- data.frame(x = seq(-1, 1, length.out = 50))
    y <- sin(x$x) + 0.1 * x$x
    bw <- npRmpi::npregbw(xdat = x, ydat = y, bws = 0.3, bandwidth.compute = FALSE)

    fit <- npRmpi::npreg(bws = bw, txdat = x, tydat = y, residuals = TRUE)
    nores <- npRmpi::npreg(bws = bw, txdat = x, tydat = y, residuals = FALSE)

    expect_equal(fit$resid, y - fit$mean, tolerance = 1e-12)
    expect_equal(fit$mean, nores$mean, tolerance = 1e-12)

    fn.body <- paste(deparse(body(getFromNamespace("npreg.rbandwidth", "npRmpi")),
                             width.cutoff = 500L), collapse = " ")
    expect_match(fn.body, "compute\\.resid\\.from\\.fit")
    expect_match(fn.body, "resid <- as\\.double\\(resid\\.response\\) - myout\\$mean")
    expect_no_match(fn.body, "myopti\\$do_grad")
    expect_no_match(fn.body, "do_grad = gradients")
  })
})

test_that("npreg factor response and external-evaluation GOF sentinels remain stable", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  local_npreg_quiet({
    tx <- data.frame(x = seq(0, 1, length.out = 36))
    ty <- factor(rep(c("a", "b", "c"), length.out = 36))
    ex <- data.frame(x = seq(0.05, 0.95, length.out = 9))
    ey <- factor(rep(c("c", "a", "b"), length.out = 9), levels = levels(ty))
    bw.factor <- npRmpi::npregbw(xdat = tx, ydat = ty, bws = 0.4, bandwidth.compute = FALSE)
    fit.factor <- npRmpi::npreg(bws = bw.factor, txdat = tx, tydat = ty, exdat = ex, eydat = ey)
    expect_equal(length(fit.factor$mean), nrow(ex))
    expect_equal(length(fit.factor$merr), nrow(ex))
    expect_null(fit.factor$xtra)
    expect_true(all(is.finite(fit.factor$mean)))

    y <- 1 + tx$x + rnorm(nrow(tx), sd = 0.01)
    bw.numeric <- npRmpi::npregbw(xdat = tx, ydat = y, bws = 0.35, bandwidth.compute = FALSE)
    fit.no.ey <- npRmpi::npreg(bws = bw.numeric, txdat = tx, tydat = y, exdat = ex)
    expect_equal(length(fit.no.ey$mean), nrow(ex))
    expect_equal(length(fit.no.ey$merr), nrow(ex))
    expect_null(fit.no.ey$xtra)
    expect_true(all(is.finite(fit.no.ey$mean)))
    expect_true(all(is.finite(fit.no.ey$merr)))
  })
})

test_that("npreg compact kernels remain finite for far external evaluation", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  local_npreg_quiet({
    tx <- data.frame(x = seq(0, 1, length.out = 40))
    ty <- sin(tx$x)
    ex <- data.frame(x = c(10, 11, 12))
    bw <- npRmpi::npregbw(
      xdat = tx,
      ydat = ty,
      bws = 0.05,
      bandwidth.compute = FALSE,
      ckertype = "uniform"
    )
    mean.fit <- npRmpi::npreg(bws = bw, txdat = tx, tydat = ty, exdat = ex)
    grad.fit <- npRmpi::npreg(bws = bw, txdat = tx, tydat = ty, exdat = ex, gradients = TRUE)

    expect_true(all(is.finite(mean.fit$mean)))
    expect_true(all(is.finite(mean.fit$merr)))
    expect_true(all(is.finite(grad.fit$grad)))
    expect_true(all(is.finite(grad.fit$gerr)))
  })
})
