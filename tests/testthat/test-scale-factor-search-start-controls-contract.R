library(np)

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
    regexp = "hbc\\.init.*max\\('lbc\\.init', 'scale\\.factor\\.lower\\.bound'\\)"
  )
}

test_that("continuous search starts reject hbc below the effective lower endpoint", {
  dat <- scale_floor_fixture()

  expect_bad_hbc_error(
    npregbw(
      y ~ x,
      data = dat,
      scale.factor.lower.bound = 1,
      lbc.init = 0.1,
      hbc.init = 0.5
    )
  )

  expect_bad_hbc_error(
    npudensbw(
      ~ x,
      data = dat,
      scale.factor.lower.bound = 1,
      lbc.init = 0.1,
      hbc.init = 0.5
    )
  )

  expect_bad_hbc_error(
    npudistbw(
      ~ x,
      data = dat,
      scale.factor.lower.bound = 1,
      lbc.init = 0.1,
      hbc.init = 0.5
    )
  )

  expect_bad_hbc_error(
    npcdensbw(
      y ~ x,
      data = dat,
      scale.factor.lower.bound = 1,
      lbc.init = 0.1,
      hbc.init = 0.5
    )
  )

  expect_bad_hbc_error(
    npcdistbw(
      y ~ x,
      data = dat,
      scale.factor.lower.bound = 1,
      lbc.init = 0.1,
      hbc.init = 0.5
    )
  )

  dat$z <- rnorm(nrow(dat))

  expect_bad_hbc_error(
    npscoefbw(
      y ~ x | z,
      data = dat,
      scale.factor.lower.bound = 1,
      lbc.init = 0.1,
      hbc.init = 0.5
    )
  )

  expect_bad_hbc_error(
    npindexbw(
      y ~ x + z,
      data = dat,
      scale.factor.lower.bound = 1,
      lbc.init = 0.1,
      hbc.init = 0.5
    )
  )
})

test_that("explicit bandwidth objects are not clamped by the search floor", {
  dat <- scale_floor_fixture()
  tiny <- 1e-8

  reg <- npregbw(
    xdat = data.frame(x = dat$x),
    ydat = dat$y,
    bws = tiny,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.lower.bound = 1
  )
  expect_equal(reg$bw[1L], tiny, tolerance = 0)

  dens <- npudensbw(
    dat = data.frame(x = dat$x),
    bws = tiny,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.lower.bound = 1
  )
  expect_equal(dens$bw[1L], tiny, tolerance = 0)

  dist <- npudistbw(
    dat = data.frame(x = dat$x),
    bws = tiny,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.lower.bound = 1
  )
  expect_equal(dist$bw[1L], tiny, tolerance = 0)

  cdens <- npcdensbw(
    xdat = data.frame(x = dat$x),
    ydat = data.frame(y = dat$y),
    bws = c(tiny, tiny),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.lower.bound = 1
  )
  expect_equal(cdens$ybw[1L], tiny, tolerance = 0)
  expect_equal(cdens$xbw[1L], tiny, tolerance = 0)

  cdist <- npcdistbw(
    xdat = data.frame(x = dat$x),
    ydat = data.frame(y = dat$y),
    bws = c(tiny, tiny),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    scale.factor.lower.bound = 1
  )
  expect_equal(cdist$ybw[1L], tiny, tolerance = 0)
  expect_equal(cdist$xbw[1L], tiny, tolerance = 0)
})
