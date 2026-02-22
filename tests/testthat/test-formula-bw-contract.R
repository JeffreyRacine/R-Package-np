test_that("formula npudensbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    x1 = runif(30),
    x2 = rnorm(30)
  )
  dat$x2[c(4, 17)] <- NA_real_
  keep <- dat$x1 > 0.2
  bw_formula <- npRmpi::npudensbw(
    ~ x1 + x2,
    data = dat,
    subset = keep,
    na.action = na.omit,
    bws = c(0.4, 0.6),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(~ x1 + x2, data = dat, subset = keep, na.action = na.omit)
  bw_default <- npRmpi::npudensbw(
    dat = mf,
    bws = c(0.4, 0.6),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npudistbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    x1 = runif(30),
    x2 = rnorm(30)
  )
  dat$x1[c(2, 11)] <- NA_real_
  keep <- dat$x2 > -0.5
  bw_formula <- npRmpi::npudistbw(
    ~ x1 + x2,
    data = dat,
    subset = keep,
    na.action = na.omit,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(~ x1 + x2, data = dat, subset = keep, na.action = na.omit)
  bw_default <- npRmpi::npudistbw(
    dat = mf,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npregbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(40),
    x1 = runif(40),
    x2 = rnorm(40)
  )
  dat$y[c(3, 9)] <- NA_real_
  keep <- dat$x1 < 0.85
  bw_formula <- npRmpi::npregbw(
    y ~ x1 + x2,
    data = dat,
    subset = keep,
    na.action = na.omit,
    bws = c(0.45, 0.55),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  mf <- model.frame(y ~ x1 + x2, data = dat, subset = keep, na.action = na.omit)
  bw_default <- npRmpi::npregbw(
    xdat = mf[, c("x1", "x2"), drop = FALSE],
    ydat = model.response(mf),
    bws = c(0.45, 0.55),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})
