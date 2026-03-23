test_that("npplreg nomad shortcut matches the explicit partially linear preset", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  n <- 22
  dat <- data.frame(x = rnorm(n), z = sort(runif(n)))
  dat$y <- 1 + 0.75 * dat$x + sin(2 * pi * dat$z) + rnorm(n, sd = 0.08)

  bw_short <- np::npplregbw(
    y ~ x | z,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )
  bw_long <- np::npplregbw(
    y ~ x | z,
    data = dat,
    regtype = "lp",
    search.engine = "nomad+powell",
    degree.select = "coordinate",
    bernstein.basis = TRUE,
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = FALSE,
    bwtype = "fixed",
    nmulti = 1L
  )

  fit_short <- np::npplreg(
    y ~ x | z,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )
  fit_long <- np::npplreg(
    y ~ x | z,
    data = dat,
    regtype = "lp",
    search.engine = "nomad+powell",
    degree.select = "coordinate",
    bernstein.basis = TRUE,
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = FALSE,
    bwtype = "fixed",
    nmulti = 1L
  )

  expect_identical(as.integer(bw_short$degree), as.integer(bw_long$degree))
  expect_true(is.list(bw_short$nomad.shortcut))
  expect_equal(fitted(fit_short), fitted(fit_long), tolerance = 5e-4)
  expect_true(is.list(fit_short$bws$nomad.shortcut))
})

test_that("npscoef nomad shortcut matches the explicit smooth-coefficient preset", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  n <- 22
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- (1 + zdat$z^2) * xdat$x + rnorm(n, sd = 0.08)

  bw_short <- np::npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )
  bw_long <- np::npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    search.engine = "nomad+powell",
    degree.select = "coordinate",
    bernstein.basis = TRUE,
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = FALSE,
    bwtype = "fixed",
    nmulti = 1L
  )

  fit_short <- np::npscoef(
    txdat = xdat,
    tzdat = zdat,
    tydat = y,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L,
    errors = FALSE,
    betas = FALSE
  )
  fit_long <- np::npscoef(
    txdat = xdat,
    tzdat = zdat,
    tydat = y,
    regtype = "lp",
    search.engine = "nomad+powell",
    degree.select = "coordinate",
    bernstein.basis = TRUE,
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = FALSE,
    bwtype = "fixed",
    nmulti = 1L,
    errors = FALSE,
    betas = FALSE
  )

  expect_identical(as.integer(bw_short$degree), as.integer(bw_long$degree))
  expect_true(is.list(bw_short$nomad.shortcut))
  expect_equal(fitted(fit_short), fitted(fit_long), tolerance = 5e-4)
  expect_true(is.list(fit_short$bws$nomad.shortcut))
})

test_that("npindex nomad shortcut matches the explicit single-index preset", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  n <- 24
  xdat <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  index <- xdat$x1 + 0.75 * xdat$x2
  y <- sin(index) + 0.25 * index^2 + rnorm(n, sd = 0.05)

  bw_short <- np::npindexbw(
    xdat = xdat,
    ydat = y,
    method = "ichimura",
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )
  bw_long <- np::npindexbw(
    xdat = xdat,
    ydat = y,
    method = "ichimura",
    regtype = "lp",
    search.engine = "nomad+powell",
    degree.select = "coordinate",
    bernstein.basis = TRUE,
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = FALSE,
    bwtype = "fixed",
    nmulti = 1L
  )

  fit_short <- np::npindex(
    txdat = xdat,
    tydat = y,
    method = "ichimura",
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )
  fit_long <- np::npindex(
    txdat = xdat,
    tydat = y,
    method = "ichimura",
    regtype = "lp",
    search.engine = "nomad+powell",
    degree.select = "coordinate",
    bernstein.basis = TRUE,
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = FALSE,
    bwtype = "fixed",
    nmulti = 1L
  )

  expect_identical(as.integer(bw_short$degree), as.integer(bw_long$degree))
  expect_true(is.list(bw_short$nomad.shortcut))
  expect_equal(fitted(fit_short), fitted(fit_long), tolerance = 5e-4)
  expect_true(is.list(fit_short$bws$nomad.shortcut))

  dat <- data.frame(y = y, x1 = xdat$x1, x2 = xdat$x2)
  fit_formula <- np::npindex(
    y ~ x1 + x2,
    data = dat,
    method = "ichimura",
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )
  expect_true(is.list(fit_formula$bws$nomad.shortcut))
  expect_true(isTRUE(fit_formula$bws$nomad.shortcut$enabled))
})

test_that("semiparametric NOMAD shortcut routes accept nomad.nmulti", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260323)
  n <- 18
  dat <- data.frame(x = rnorm(n), z = sort(runif(n)))
  dat$y <- 1 + 0.5 * dat$x + sin(2 * pi * dat$z) + rnorm(n, sd = 0.05)

  fit_pl <- np::npplreg(
    y ~ x | z,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L,
    nomad.nmulti = 2L
  )
  fit_sc <- np::npscoef(
    txdat = data.frame(x = dat$x),
    tzdat = data.frame(z = dat$z),
    tydat = dat$y,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L,
    nomad.nmulti = 2L,
    errors = FALSE,
    betas = FALSE
  )
  fit_si <- np::npindex(
    txdat = data.frame(x1 = dat$x, x2 = dat$z),
    tydat = dat$y,
    method = "ichimura",
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L,
    nomad.nmulti = 2L
  )

  expect_true(is.list(fit_pl$bws$nomad.shortcut))
  expect_true(is.list(fit_sc$bws$nomad.shortcut))
  expect_true(is.list(fit_si$bws$nomad.shortcut))
})

test_that("semiparametric nomad shortcuts fail fast on incompatible settings", {
  xdat <- data.frame(x = runif(10))
  zdat <- data.frame(z = runif(10))
  y <- runif(10)

  expect_error(
    np::npplregbw(xdat = xdat, zdat = zdat, ydat = y, nomad = TRUE, regtype = "ll"),
    "nomad=TRUE requires regtype='lp'"
  )
  expect_error(
    np::npscoefbw(xdat = xdat, zdat = zdat, ydat = y, nomad = TRUE, bernstein.basis = FALSE),
    "requires bernstein.basis=TRUE"
  )
  expect_error(
    np::npindexbw(xdat = data.frame(x1 = runif(10), x2 = runif(10)), ydat = y, nomad = TRUE, bwtype = "adaptive_nn"),
    "requires bwtype='fixed'"
  )
  expect_error(
    np::npplregbw(xdat = xdat, zdat = zdat, ydat = y, nomad.nmulti = 1L),
    "nomad.nmulti is only supported when regtype='lp', automatic degree search is active, and search.engine is 'nomad' or 'nomad\\+powell'"
  )
})
