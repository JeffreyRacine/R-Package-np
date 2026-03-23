test_that("npplregbw exhaustive degree search matches manual profile minimum", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 28
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- 1 + 0.75 * xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.08)

  bw0 <- npplregbw(
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
  bw1 <- npplregbw(
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
  auto <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(auto, "plbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, min(bw0$fval, bw1$fval) + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))

  manual <- npplregbw(
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

test_that("npplregbw coordinate search can be exhaustively certified on a small grid", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 26
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(
    z1 = runif(n),
    z2 = runif(n)
  )
  y <- 1 + 0.5 * xdat$x + sin(2 * pi * zdat$z1) + zdat$z2^2 + rnorm(n, sd = 0.08)

  exhaustive <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  coordinate <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "cell",
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

test_that("npplregbw automatic degree search enforces pilot guardrails", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = runif(n))
  y <- 1 + xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.08)

  expect_error(
    npplregbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lc",
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "automatic degree search currently requires regtype='lp'"
  )

  expect_error(
    npplregbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lp",
      bandwidth.compute = FALSE,
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bws = matrix(0.2, nrow = 2L, ncol = 1L)
    ),
    "bandwidth.compute=TRUE"
  )

  bw <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    bernstein.basis = FALSE,
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 4L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(bw, "plbandwidth")
  expect_false(isTRUE(bw$bernstein.basis))
  expect_lte(max(as.integer(bw$degree)), 4L)
})

test_that("npplreg forwards automatic LP degree search through npplregbw", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  dat <- data.frame(
    x = rnorm(n),
    z = runif(n)
  )
  dat$y <- 1 + 0.75 * dat$x + sin(2 * pi * dat$z) + rnorm(n, sd = 0.08)

  fit <- npplreg(
    y ~ x | z,
    data = dat,
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(fit, "plregression")
  expect_s3_class(fit$bws, "plbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})

test_that("npplregbw automatic degree search defaults to NOMAD plus Powell", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = runif(n))
  y <- 1 + 0.75 * xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.08)

  auto <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(auto, "plbandwidth")
  expect_false(is.null(auto$degree.search))
  expect_identical(auto$degree.search$mode, "nomad+powell")
  expect_equal(length(auto$degree), 1L)
  expect_true(is.finite(auto$fval))
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-8)
  expect_lte(auto$fval, auto$degree.search$baseline.fval + 1e-8)
  expect_equal(auto$fval, auto$degree.search$best.fval, tolerance = 1e-8)
})

test_that("npplreg explicit plbandwidth route preserves NOMAD child payload names", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  n <- 40L
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1^2 + rnorm(n, sd = 0.1)

  bw <- npplregbw(
    xdat = data.frame(x1 = x1),
    zdat = data.frame(x2 = x2),
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad+powell",
    degree.max = 2L,
    nmulti = 1L
  )

  expect_identical(names(bw$bw), c("yzbw", "x1"))

  fit <- npplreg(
    bws = bw,
    txdat = data.frame(x1 = x1),
    tzdat = data.frame(x2 = x2),
    tydat = y
  )

  expect_s3_class(fit, "plregression")
  expect_equal(nrow(fit$evalx), n)
})
