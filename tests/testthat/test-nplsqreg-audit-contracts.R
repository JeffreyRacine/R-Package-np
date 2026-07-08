skip_slow_nplsqreg_refined_tau_search <- function() {
  skip_on_cran()
}

test_that("npRmpi nplsqreg direct routes omit incomplete training rows and realign output", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260703)
  n <- 24L
  dat <- data.frame(
    x = seq(0.05, 0.95, length.out = n),
    y = sin(seq(0.05, 0.95, length.out = n) * 2 * pi)
  )
  dat$x[5L] <- NA_real_

  bw <- nplsqregbw(
    xdat = dat["x"], ydat = dat$y,
    bws = 0.25, bandwidth.compute = FALSE,
    scale = rep(1, n), tau = 0.5
  )
  expect_identical(bw$rows.omit, 5L)
  expect_identical(bw$nobs.omit, 1L)
  expect_equal(nrow(bw$xdat), n - 1L)

  fit <- nplsqreg(
    txdat = dat["x"], tydat = dat$y,
    bws = 0.25, bandwidth.compute = FALSE,
    scale = rep(1, n), tau = c(0.25, 0.5),
    gradients = TRUE, residuals = TRUE
  )
  expect_equal(dim(fit$quantile), c(n, 2L))
  expect_true(all(is.na(fit$quantile[5L, ])))
  expect_true(all(is.na(fit$quanterr[5L, ])))
  expect_true(all(is.na(fit$quantgrad[5L, , ])))
  expect_true(all(is.na(fit$quantgerr[5L, , ])))
  expect_true(is.na(fit$resid[5L]))
  expect_identical(fit$rows.omit, 5L)
  expect_identical(fit$nobs.omit, 1L)
})

test_that("npRmpi nplsqreg formula routes preserve training and evaluation row shape", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260703)
  n <- 24L
  dat <- data.frame(
    x = seq(0.05, 0.95, length.out = n),
    y = sin(seq(0.05, 0.95, length.out = n) * 2 * pi)
  )
  dat.train.na <- dat
  dat.train.na$y[6L] <- NA_real_
  complete.n <- sum(stats::complete.cases(dat.train.na))

  fit.train <- nplsqreg(
    y ~ x, data = dat.train.na,
    scale = rep(1, complete.n),
    tau = c(0.25, 0.5),
    gradients = TRUE, residuals = TRUE,
    nmulti = 1L, itmax = 20L,
    na.action = na.omit
  )
  expect_equal(dim(fit.train$quantile), c(n, 2L))
  expect_true(all(is.na(fit.train$quantile[6L, ])))
  expect_true(all(is.na(fit.train$quantgrad[6L, , ])))
  expect_true(is.na(fit.train$resid[6L]))
  expect_identical(fit.train$rows.omit, 6L)
  expect_identical(fit.train$nobs.omit, 1L)

  eval.dat <- data.frame(x = c(0.2, NA, 0.8))
  fit.eval <- nplsqreg(
    y ~ x, data = dat, newdata = eval.dat,
    scale = rep(1, n),
    tau = c(0.25, 0.5),
    gradients = TRUE,
    nmulti = 1L, itmax = 20L,
    na.action = na.omit
  )
  expect_equal(dim(fit.eval$quantile), c(3L, 2L))
  expect_true(all(is.na(fit.eval$quantile[2L, ])))
  expect_true(all(is.na(fit.eval$quantgrad[2L, , ])))
  expect_identical(fit.eval$eval.rows.omit, 2L)
  expect_identical(fit.eval$eval.nobs.omit, 1L)
  expect_identical(fit.eval$nobs, 3L)
})

test_that("npRmpi nplsqreg exact bandwidth reuse rejects mismatched explicit training data", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260703)
  n <- 24L
  dat <- data.frame(
    x = seq(0.05, 0.95, length.out = n),
    y = sin(seq(0.05, 0.95, length.out = n) * 2 * pi)
  )
  fit <- nplsqreg(
    txdat = dat["x"], tydat = dat$y,
    bws = 0.25, bandwidth.compute = FALSE,
    scale = rep(1, n), tau = 0.5
  )
  expect_equal(fitted(nplsqreg(fit$bws)), fitted(fit), tolerance = 0)
  expect_equal(
    fitted(nplsqreg(bws = fit$bws, txdat = dat["x"], tydat = dat$y)),
    fitted(fit),
    tolerance = 0
  )
  changed <- dat
  changed$x <- rev(changed$x)
  expect_error(
    nplsqreg(bws = fit$bws, txdat = changed["x"], tydat = dat$y),
    "do not match"
  )
})

test_that("npRmpi nplsqreg fixed delta layout is preserved by the check-loss backend", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260703)
  n <- 24L
  dat <- data.frame(
    x = seq(0.05, 0.95, length.out = n),
    y = sin(seq(0.05, 0.95, length.out = n) * 2 * pi)
  )
  scale <- seq(1.0, 1.5, length.out = n)
  delta <- 0.37
  bw <- nplsqregbw(
    xdat = dat["x"], ydat = dat$y,
    bws = 0.25, bandwidth.compute = FALSE,
    scale = scale, delta = delta, tau = 0.5
  )
  expect_equal(bw$delta, delta, tolerance = 1e-14)
  expect_equal(bw$qdat, dat$y + scale * stats::qnorm(delta),
               tolerance = 1e-12)
})

test_that("npRmpi nplsqreg refined tau search keeps full controls for the anchor", {
  skip_slow_nplsqreg_refined_tau_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260703)
  n <- 24L
  dat <- data.frame(
    x = seq(0.05, 0.95, length.out = n),
    y = sin(seq(0.05, 0.95, length.out = n) * 2 * pi)
  )
  bw <- nplsqregbw(
    xdat = dat["x"], ydat = dat$y,
    scale = rep(1, n),
    tau = c(0.25, 0.5, 0.75),
    tau.search = "refined",
    nmulti = 3L, powell.remin = TRUE, itmax = 20L
  )
  expect_identical(bw$fit.order[1L], 2L)
  expect_identical(bw$tau.search.controls$nmulti, 1L)
  expect_identical(bw$tau.search.controls$powell.remin, FALSE)
  expect_identical(bw$tau.search.controls$anchor.controls$nmulti, 3L)
  expect_identical(bw$tau.search.controls$anchor.controls$powell.remin, TRUE)
})
