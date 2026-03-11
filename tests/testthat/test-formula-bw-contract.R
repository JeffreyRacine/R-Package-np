test_that("formula npudensbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    x1 = runif(30),
    x2 = rnorm(30)
  )
  dat$x2[c(4, 17)] <- NA_real_
  bw_formula <- npRmpi::npudensbw(
    ~ x1 + x2,
    data = dat,
    subset = x1 > 0.2,
    na.action = na.omit,
    bws = c(0.4, 0.6),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(~ x1 + x2, data = dat, subset = x1 > 0.2, na.action = na.omit)
  bw_default <- npRmpi::npudensbw(
    dat = mf,
    bws = c(0.4, 0.6),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npudistbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    x1 = runif(30),
    x2 = rnorm(30)
  )
  dat$x1[c(2, 11)] <- NA_real_
  bw_formula <- npRmpi::npudistbw(
    ~ x1 + x2,
    data = dat,
    subset = x2 > -0.5,
    na.action = na.omit,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(~ x1 + x2, data = dat, subset = x2 > -0.5, na.action = na.omit)
  bw_default <- npRmpi::npudistbw(
    dat = mf,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npudistbw gdata path matches default interface", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    x1 = runif(30),
    x2 = rnorm(30)
  )
  gdat <- data.frame(
    x1 = runif(22),
    x2 = rnorm(22)
  )

  bw_formula <- npRmpi::npudistbw(
    ~ x1 + x2,
    data = dat,
    gdata = gdat,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(~ x1 + x2, data = dat)
  gmf <- model.frame(~ x1 + x2, data = gdat)
  bw_default <- npRmpi::npudistbw(
    dat = mf,
    gdat = gmf,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npregbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(40),
    x1 = runif(40),
    x2 = rnorm(40)
  )
  dat$y[c(3, 9)] <- NA_real_
  bw_formula <- npRmpi::npregbw(
    y ~ x1 + x2,
    data = dat,
    subset = x1 < 0.85,
    na.action = na.omit,
    bws = c(0.45, 0.55),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  mf <- model.frame(y ~ x1 + x2, data = dat, subset = x1 < 0.85, na.action = na.omit)
  bw_default <- npRmpi::npregbw(
    xdat = mf[, c("x1", "x2"), drop = FALSE],
    ydat = model.response(mf),
    bws = c(0.45, 0.55),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npcdistbw gdata path matches default interface", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(36),
    x = runif(36)
  )
  gdat <- data.frame(
    y = rnorm(18),
    x = runif(18)
  )

  bw_formula <- npRmpi::npcdistbw(
    y ~ x,
    data = dat,
    gdata = gdat,
    bws = c(0.35, 0.45),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(y ~ x, data = dat)
  gmf <- model.frame(y ~ x, data = gdat)
  bw_default <- npRmpi::npcdistbw(
    xdat = mf[, "x", drop = FALSE],
    ydat = mf[, "y", drop = FALSE],
    gydat = gmf[, "y", drop = FALSE],
    bws = c(0.35, 0.45),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npcdensbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(36),
    x = runif(36)
  )
  dat$y[c(5, 19)] <- NA_real_

  bw_formula <- npRmpi::npcdensbw(
    y ~ x,
    data = dat,
    subset = x > 0.25,
    na.action = na.omit,
    bws = c(0.35, 0.45),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(y ~ x, data = dat, subset = x > 0.25, na.action = na.omit)
  bw_default <- npRmpi::npcdensbw(
    xdat = mf[, "x", drop = FALSE],
    ydat = mf[, "y", drop = FALSE],
    bws = c(0.35, 0.45),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npindexbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(42),
    x1 = runif(42),
    x2 = rnorm(42)
  )
  dat$y[c(4, 21)] <- NA_real_

  bw_formula <- npRmpi::npindexbw(
    y ~ x1 + x2,
    data = dat,
    subset = x1 < 0.9,
    na.action = na.omit,
    bws = c(0.2, 0.3, 0.4),
    bandwidth.compute = FALSE,
    method = "ichimura",
    nmulti = 1
  )

  mf <- model.frame(y ~ x1 + x2, data = dat, subset = x1 < 0.9, na.action = na.omit)
  bw_default <- npRmpi::npindexbw(
    xdat = mf[, c("x1", "x2"), drop = FALSE],
    ydat = model.response(mf),
    bws = c(0.2, 0.3, 0.4),
    bandwidth.compute = FALSE,
    method = "ichimura",
    nmulti = 1
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npplregbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(44),
    x = runif(44),
    z = rnorm(44)
  )
  bws_mat <- matrix(c(0.25, 0.35), nrow = 2, ncol = 1)

  bw_formula <- npRmpi::npplregbw(
    y ~ z | x,
    data = dat,
    subset = x < 0.92,
    na.action = na.omit,
    bws = bws_mat,
    bandwidth.compute = FALSE
  )

  mf <- model.frame(y ~ x + z, data = dat, subset = x < 0.92, na.action = na.omit)
  bw_default <- npRmpi::npplregbw(
    xdat = mf[, "z", drop = FALSE],
    ydat = model.response(mf),
    zdat = mf[, "x", drop = FALSE],
    bws = bws_mat,
    bandwidth.compute = FALSE
  )

  bw_formula_vec <- vapply(bw_formula$bw, function(bwi) as.numeric(bwi$bw)[1L], numeric(1))
  bw_default_vec <- vapply(bw_default$bw, function(bwi) as.numeric(bwi$bw)[1L], numeric(1))
  expect_equal(bw_formula_vec, bw_default_vec)
})

test_that("formula npscoefbw matches default interface with subset/na.action", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(44),
    x = runif(44),
    z = rnorm(44)
  )
  dat$y[c(8, 29)] <- NA_real_

  bw_formula <- npRmpi::npscoefbw(
    y ~ x | z,
    data = dat,
    subset = x > 0.08,
    na.action = na.omit,
    bws = 0.4,
    bandwidth.compute = FALSE
  )

  mf <- model.frame(y ~ x + z, data = dat, subset = x > 0.08, na.action = na.omit)
  bw_default <- npRmpi::npscoefbw(
    xdat = mf[, "x", drop = FALSE],
    ydat = model.response(mf),
    zdat = mf[, "z", drop = FALSE],
    bws = 0.4,
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})
