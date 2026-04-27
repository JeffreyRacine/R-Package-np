suppressPackageStartupMessages(library(npRmpi))

scale_floor_fixture <- function(n = 20L) {
  set.seed(8675309)
  data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )
}

expect_bad_hbc_error <- function(expr) {
  expect_error(
    expr,
    regexp = "scale\\.factor\\.init\\.upper.*max\\('scale\\.factor\\.init\\.lower', 'scale\\.factor\\.search\\.lower'\\)"
  )
}

test_that("continuous search starts reject hbc below the effective lower endpoint", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- scale_floor_fixture()

  expect_bad_hbc_error(
    npregbw(
      y ~ x,
      data = dat,
      scale.factor.search.lower = 1,
      scale.factor.init.lower = 0.1,
      scale.factor.init.upper = 0.5
    )
  )

  expect_bad_hbc_error(
    npudensbw(
      ~ x,
      data = dat,
      scale.factor.search.lower = 1,
      scale.factor.init.lower = 0.1,
      scale.factor.init.upper = 0.5
    )
  )

  expect_bad_hbc_error(
    npudistbw(
      ~ x,
      data = dat,
      scale.factor.search.lower = 1,
      scale.factor.init.lower = 0.1,
      scale.factor.init.upper = 0.5
    )
  )

  expect_bad_hbc_error(
    npcdensbw(
      y ~ x,
      data = dat,
      scale.factor.search.lower = 1,
      scale.factor.init.lower = 0.1,
      scale.factor.init.upper = 0.5
    )
  )

  expect_bad_hbc_error(
    npcdistbw(
      y ~ x,
      data = dat,
      scale.factor.search.lower = 1,
      scale.factor.init.lower = 0.1,
      scale.factor.init.upper = 0.5
    )
  )

  dat$z <- rnorm(nrow(dat))

  expect_bad_hbc_error(
    npscoefbw(
      y ~ x | z,
      data = dat,
      scale.factor.search.lower = 1,
      scale.factor.init.lower = 0.1,
      scale.factor.init.upper = 0.5
    )
  )

  expect_bad_hbc_error(
    npindexbw(
      y ~ x + z,
      data = dat,
      scale.factor.search.lower = 1,
      scale.factor.init.lower = 0.1,
      scale.factor.init.upper = 0.5
    )
  )
})

test_that("explicit bandwidth objects are not clamped by the search floor", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- scale_floor_fixture()
  tiny <- 1e-8

  reg <- npregbw(
    xdat = data.frame(x = dat$x),
    ydat = dat$y,
    bws = tiny,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.search.lower = 1
  )
  expect_equal(reg$bw[1L], tiny, tolerance = 0)

  dens <- npudensbw(
    dat = data.frame(x = dat$x),
    bws = tiny,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.search.lower = 1
  )
  expect_equal(dens$bw[1L], tiny, tolerance = 0)

  dist <- npudistbw(
    dat = data.frame(x = dat$x),
    bws = tiny,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.search.lower = 1
  )
  expect_equal(dist$bw[1L], tiny, tolerance = 0)

  cdens <- npcdensbw(
    xdat = data.frame(x = dat$x),
    ydat = data.frame(y = dat$y),
    bws = c(tiny, tiny),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.search.lower = 1
  )
  expect_equal(cdens$ybw[1L], tiny, tolerance = 0)
  expect_equal(cdens$xbw[1L], tiny, tolerance = 0)

  cdist <- npcdistbw(
    xdat = data.frame(x = dat$x),
    ydat = data.frame(y = dat$y),
    bws = c(tiny, tiny),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.search.lower = 1
  )
  expect_equal(cdist$ybw[1L], tiny, tolerance = 0)
  expect_equal(cdist$xbw[1L], tiny, tolerance = 0)
})
