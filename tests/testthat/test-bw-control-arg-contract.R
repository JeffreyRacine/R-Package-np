test_that("bandwidth wrappers reject invalid scalar control flags", {
  dat <- data.frame(x = c(0.1, 0.2))
  bws <- list()

  expect_error(
    npudensbw.bandwidth(dat = dat, bws = bws, bandwidth.compute = c(TRUE, FALSE)),
    "'bandwidth.compute' must be TRUE or FALSE"
  )
  expect_error(
    npudensbw.bandwidth(dat = dat, bws = bws, nmulti = -1),
    "'nmulti' must be a non-negative integer"
  )
  expect_error(
    npudensbw.bandwidth(dat = dat, bws = bws, itmax = 0),
    "'itmax' must be a positive integer"
  )
  expect_error(
    npudensbw.bandwidth(dat = dat, bws = bws, ftol = 0),
    "'ftol' must be a positive finite numeric scalar"
  )
  expect_error(
    npudensbw.bandwidth(dat = dat, bws = bws, penalty.multiplier = 0),
    "'penalty.multiplier' must be a positive finite numeric scalar"
  )

  expect_error(
    npudistbw.dbandwidth(dat = dat, bws = bws, do.full.integral = c(TRUE, FALSE)),
    "'do.full.integral' must be TRUE or FALSE"
  )
  expect_error(
    npudistbw.dbandwidth(dat = dat, bws = bws, ngrid = 0),
    "'ngrid' must be a positive integer"
  )
  expect_error(
    npudistbw.dbandwidth(dat = dat, bws = bws, memfac = 0),
    "'memfac' must be a positive finite numeric scalar"
  )

  expect_error(
    npregbw.rbandwidth(xdat = dat, ydat = dat$x, bws = bws, bandwidth.compute = c(TRUE, FALSE)),
    "'bandwidth.compute' must be TRUE or FALSE"
  )
  expect_error(
    npregbw.rbandwidth(xdat = dat, ydat = dat$x, bws = bws, itmax = 0),
    "'itmax' must be a positive integer"
  )

  expect_error(
    npcdensbw.conbandwidth(xdat = dat, ydat = dat, bws = bws, memfac = 0),
    "'memfac' must be a positive finite numeric scalar"
  )
  expect_error(
    npcdensbw.conbandwidth(xdat = dat, ydat = dat, bws = bws, nmulti = -1),
    "'nmulti' must be a non-negative integer"
  )

  expect_error(
    npcdistbw.condbandwidth(xdat = dat, ydat = dat, bws = bws, do.full.integral = c(TRUE, FALSE)),
    "'do.full.integral' must be TRUE or FALSE"
  )
  expect_error(
    npcdistbw.condbandwidth(xdat = dat, ydat = dat, bws = bws, ngrid = 0),
    "'ngrid' must be a positive integer"
  )

  expect_error(
    npindexbw.sibandwidth(xdat = dat, ydat = dat$x, bws = bws, bandwidth.compute = c(TRUE, FALSE)),
    "'bandwidth.compute' must be TRUE or FALSE"
  )
  expect_error(
    npindexbw.sibandwidth(xdat = dat, ydat = dat$x, bws = bws, optim.maxit = 0),
    "'optim.maxit' must be a positive integer"
  )

  expect_error(
    npscoefbw.scbandwidth(xdat = dat, ydat = dat$x, bws = bws, cv.iterate = c(TRUE, FALSE)),
    "'cv.iterate' must be TRUE or FALSE"
  )
  expect_error(
    npscoefbw.scbandwidth(xdat = dat, ydat = dat$x, bws = bws, backfit.tol = 0),
    "'backfit.tol' must be a positive finite numeric scalar"
  )
})
