test_that("npindex lc eval objective localizes the large-h leave-one-out leaf", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  x1 <- runif(300L, -1, 1)
  x2 <- runif(300L, -1, 1)
  y <- x1 + rnorm(300L, sd = 0.1)
  xdat <- data.frame(x1 = x1, x2 = x2)
  xmat <- as.matrix(xdat)

  bw <- npindexbw(xdat = xdat, ydat = y, regtype = "lc", nmulti = 1L)
  spec <- npRmpi:::.npindex_resolve_spec(bw)
  param <- c(as.double(bw$beta[2L]), 1e10)

  actual <- npRmpi:::.npindexbw_eval_objective(param, xmat, y, bw, spec)

  index <- xmat %*% c(1, as.double(param[1L]))
  wmat <- cbind(y, 1.0)
  tww <- npRmpi:::.npRmpi_with_local_regression(
    npksum(
      txdat = index,
      tydat = wmat,
      weights = wmat,
      leave.one.out = TRUE,
      bandwidth.divide = TRUE,
      bws = c(1e10),
      bwtype = bw$type,
      ckertype = bw$ckertype,
      ckerorder = bw$ckerorder,
      ckerbound = bw$ckerbound,
      ckerlb = bw$ckerlb,
      ckerub = bw$ckerub
    )
  )$ksum
  fit.loo <- tww[1, 2, ] / npRmpi:::NZD(tww[2, 2, ])
  expected.objective <- mean((y - fit.loo)^2)
  expected.fast <- if (npRmpi:::.npindexbw_fast_eligible(h = 1e10, bws = bw, eval.index = index)) 1L else 0L

  expect_equal(as.numeric(actual$objective), as.numeric(expected.objective), tolerance = 1e-12)
  expect_identical(as.integer(actual$num.feval.fast), as.integer(expected.fast))
  expect_identical(as.integer(actual$num.feval.fast), 1L)
})

test_that("npindex lc object-fed huge-h fit stays finite and constant under MPI", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  x1 <- runif(300L, -1, 1)
  x2 <- runif(300L, -1, 1)
  y <- x1 + rnorm(300L, sd = 0.1)
  xdat <- data.frame(x1 = x1, x2 = x2)

  bw <- npindexbw(xdat = xdat, ydat = y, regtype = "lc", nmulti = 1L)
  bw$bw <- 1e10
  bw$bandwidth[[1L]] <- 1e10

  fit <- npindex(bws = bw, txdat = xdat, tydat = y)

  expect_true(all(is.finite(fit$mean)))
  expect_equal(max(fit$mean) - min(fit$mean), 0, tolerance = 1e-12)
  expect_true(is.finite(fit$MSE))
})
