test_that("bandwidth wrappers reject invalid scalar control flags", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- data.frame(x = c(0.1, 0.2))
  bws <- list()
  ns_fun <- function(name) get(name, envir = asNamespace("npRmpi"))
  npudensbw_bandwidth <- ns_fun("npudensbw.bandwidth")
  npudistbw_dbandwidth <- ns_fun("npudistbw.dbandwidth")
  npregbw_rbandwidth <- ns_fun("npregbw.rbandwidth")
  npcdensbw_conbandwidth <- ns_fun("npcdensbw.conbandwidth")
  npcdistbw_condbandwidth <- ns_fun("npcdistbw.condbandwidth")
  npindexbw_sibandwidth <- ns_fun("npindexbw.sibandwidth")
  npscoefbw_scbandwidth <- ns_fun("npscoefbw.scbandwidth")
  npplregbw_default <- ns_fun("npplregbw.default")
  npplregbw_plbandwidth <- ns_fun("npplregbw.plbandwidth")
  cdens.bws <- npcdensbw(xdat = dat, ydat = dat, bws = c(1, 1),
                         bandwidth.compute = FALSE)
  cdist.bws <- npcdistbw(xdat = dat, ydat = dat, bws = c(1, 1),
                         bandwidth.compute = FALSE)

  expect_error(
    npudensbw_bandwidth(dat = dat, bws = bws, bandwidth.compute = c(TRUE, FALSE)),
    "'bandwidth.compute' must be TRUE or FALSE"
  )
  expect_error(
    npudensbw_bandwidth(dat = dat, bws = bws, nmulti = -1),
    "'nmulti' must be a positive integer"
  )
  expect_error(
    npudensbw_bandwidth(dat = dat, bws = bws, itmax = 0),
    "'itmax' must be a positive integer"
  )
  expect_error(
    npudensbw_bandwidth(dat = dat, bws = bws, ftol = 0),
    "'ftol' must be a positive finite numeric scalar"
  )
  expect_error(
    npudensbw_bandwidth(dat = dat, bws = bws, penalty.multiplier = 0),
    "'penalty.multiplier' must be a positive finite numeric scalar"
  )

  expect_error(
    npudistbw_dbandwidth(dat = dat, bws = bws, do.full.integral = c(TRUE, FALSE)),
    "'do.full.integral' must be TRUE or FALSE"
  )
  expect_error(
    npudistbw_dbandwidth(dat = dat, bws = bws, ngrid = 0),
    "'ngrid' must be a positive integer"
  )
  expect_error(
    npudistbw_dbandwidth(dat = dat, bws = bws, memfac = 0),
    "'memfac' must be a positive finite numeric scalar"
  )

  expect_error(
    npregbw_rbandwidth(xdat = dat, ydat = dat$x, bws = bws, bandwidth.compute = c(TRUE, FALSE)),
    "'bandwidth.compute' must be TRUE or FALSE"
  )
  expect_error(
    npregbw_rbandwidth(xdat = dat, ydat = dat$x, bws = bws, itmax = 0),
    "'itmax' must be a positive integer"
  )

  expect_error(
    npcdensbw_conbandwidth(xdat = dat, ydat = dat, bws = cdens.bws, memfac = 0),
    "'memfac' must be a positive finite numeric scalar"
  )
  expect_error(
    npcdensbw_conbandwidth(xdat = dat, ydat = dat, bws = cdens.bws, nmulti = -1),
    "'nmulti' must be a positive integer"
  )

  expect_error(
    npcdistbw_condbandwidth(xdat = dat, ydat = dat, bws = cdist.bws, do.full.integral = c(TRUE, FALSE)),
    "'do.full.integral' must be TRUE or FALSE"
  )
  expect_error(
    npcdistbw_condbandwidth(xdat = dat, ydat = dat, bws = cdist.bws, ngrid = 0),
    "'ngrid' must be a positive integer"
  )

  expect_error(
    npindexbw_sibandwidth(xdat = dat, ydat = dat$x, bws = bws, bandwidth.compute = c(TRUE, FALSE)),
    "'bandwidth.compute' must be TRUE or FALSE"
  )
  expect_error(
    npindexbw_sibandwidth(xdat = dat, ydat = dat$x, bws = bws, optim.maxit = 0),
    "'optim.maxit' must be a positive integer"
  )
  expect_error(
    npindexbw_sibandwidth(xdat = dat, ydat = dat$x, bws = bws, scale.factor.init.lower = 0.6, scale.factor.init.upper = 0.5),
    "'scale.factor.init.upper' must be greater than or equal to max"
  )

  expect_error(
    npscoefbw_scbandwidth(xdat = dat, ydat = dat$x, bws = bws, cv.iterate = c(TRUE, FALSE)),
    "'cv.iterate' must be TRUE or FALSE"
  )
  expect_error(
    npscoefbw_scbandwidth(xdat = dat, ydat = dat$x, bws = bws, backfit.tol = 0),
    "'backfit.tol' must be a positive finite numeric scalar"
  )
  expect_error(
    npscoefbw_scbandwidth(xdat = dat, ydat = dat$x, bws = bws, scale.factor.init.lower = 0.8, scale.factor.init.upper = 0.6),
    "'scale.factor.init.upper' must be greater than or equal to max"
  )
  expect_error(
    npscoefbw_scbandwidth(xdat = dat, ydat = dat$x, bws = bws, hbd.init = 2.1),
    "categorical start factors must be less than or equal to 2"
  )

  expect_error(
    npplregbw_default(xdat = dat, ydat = dat$x, zdat = dat,
                      bws = matrix(0, nrow = 2, ncol = 1),
                      bandwidth.compute = c(TRUE, FALSE)),
    "'bandwidth.compute' must be TRUE or FALSE"
  )
  expect_error(
    npplregbw_plbandwidth(xdat = dat, ydat = dat$x, zdat = dat, bws = bws, nmulti = -1),
    "'nmulti' must be a positive integer"
  )
})
