test_that("npscoefbw exhaustive degree search matches manual profile minimum", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 28
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- (1 + zdat$z^2) * xdat$x + rnorm(n, sd = 0.08)

  bw0 <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree = 0L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  bw1 <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree = 1L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  auto <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(auto, "scbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, min(bw0$fval, bw1$fval) + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))

  manual <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  expect_null(manual$degree.search)
})

test_that("npscoefbw coordinate search can be exhaustively certified on a small grid", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 26
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(
    z1 = runif(n),
    z2 = runif(n)
  )
  y <- (1 + zdat$z1 + zdat$z2^2) * xdat$x + rnorm(n, sd = 0.08)

  exhaustive <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  coordinate <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = TRUE,
    degree.restarts = 1L,
    degree.max.cycles = 4L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(coordinate$degree.search$mode, "coordinate")
  expect_true(isTRUE(coordinate$degree.search$completed))
  expect_true(isTRUE(coordinate$degree.search$certified))
  expect_equal(as.integer(coordinate$degree), as.integer(exhaustive$degree))
  expect_equal(coordinate$fval, exhaustive$fval, tolerance = 1e-10)
})

test_that("npscoefbw automatic degree search enforces pilot guardrails", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  y <- (1 + zdat$z) * xdat$x + rnorm(n, sd = 0.08)

  expect_error(
    npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lc",
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "automatic degree search currently requires regtype='lp'"
  )

  expect_error(
    npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lp",
      bandwidth.compute = FALSE,
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 1L,
      bws = 0.2
    ),
    "bandwidth.compute=TRUE"
  )

  expect_error(
    npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lp",
      bernstein.basis = FALSE,
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 4L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "degree.max <= 3"
  )

  expect_error(
    npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lp",
      degree.select = "coordinate",
      degree.min = 0L,
      degree.max = 1L,
      cv.iterate = TRUE,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "cv.iterate=FALSE"
  )
})

test_that("npscoef forwards automatic LP degree search through npscoefbw", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  dat <- data.frame(
    x = runif(n),
    z = runif(n)
  )
  dat$y <- (1 + dat$z^2) * dat$x + rnorm(n, sd = 0.08)

  fit <- npscoef(
    y ~ x | z,
    data = dat,
    regtype = "lp",
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(fit, "smoothcoefficient")
  expect_s3_class(fit$bws, "scbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})
