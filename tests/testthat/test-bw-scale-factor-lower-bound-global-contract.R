suppressPackageStartupMessages(library(npRmpi))

test_that("omitted generic scale-factor floor defaults to 0.1", {
  expect_equal(npRmpi:::npResolveScaleFactorLowerBound(NULL), 0.1, tolerance = 0)
  expect_equal(npRmpi:::npResolveScaleFactorLowerBound(0.01), 0.01, tolerance = 0)
})

test_that("scale.factor.search.lower is accepted across npRmpi bandwidth selectors", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260423)
  n <- 32L
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.25)

  expect_floor <- function(object, expected = 0.2) {
    expect_equal(object$scale.factor.search.lower, expected, tolerance = 0)
    invisible(object)
  }

  common <- list(
    bwscaling = TRUE,
    bwtype = "fixed",
    nmulti = 1L,
    itmax = 10L,
    scale.factor.search.lower = 0.2
  )

  expect_floor(do.call(npregbw, c(
    list(xdat = data.frame(x = x), ydat = y, bwmethod = "cv.ls",
         regtype = "lp", degree = 0L),
    common
  )))

  expect_floor(do.call(npudensbw, c(
    list(dat = data.frame(y = y), bwmethod = "cv.ls"),
    common
  )))

  expect_floor(do.call(npudistbw, c(
    list(dat = data.frame(y = y), bwmethod = "cv.cdf"),
    common
  )))

  expect_floor(do.call(npcdensbw, c(
    list(xdat = data.frame(x = x), ydat = y, bwmethod = "cv.ls",
         regtype = "lp", degree = 0L),
    common
  )))

  expect_floor(do.call(npcdistbw, c(
    list(xdat = data.frame(x = x), ydat = y, bwmethod = "cv.ls",
         regtype = "lp", degree = 0L),
    common
  )))

  pl <- npplregbw(
    xdat = data.frame(x = x),
    ydat = y,
    zdat = data.frame(z = z),
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1L,
    itmax = 10L,
    scale.factor.search.lower = 0.2
  )
  expect_equal(
    unname(vapply(pl$bw, function(bwi) bwi$scale.factor.search.lower, numeric(1L))),
    rep(0.2, length(pl$bw)),
    tolerance = 0
  )

  si <- npindexbw(
    xdat = data.frame(x = x, z = z),
    ydat = y,
    bws = c(0, 0, 0),
    bwmethod = "ichimura",
    bwtype = "fixed",
    nmulti = 1L,
    optim.maxit = 10L,
    scale.factor.search.lower = 0.2
  )
  expect_equal(si$scale.factor.search.lower, 0.2, tolerance = 0)

  sc <- npscoefbw(
    xdat = data.frame(x = x),
    ydat = y,
    zdat = data.frame(z = z),
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1L,
    optim.maxit = 10L,
    scale.factor.search.lower = 0.2
  )
  expect_equal(sc$scale.factor.search.lower, 0.2, tolerance = 0)
})
