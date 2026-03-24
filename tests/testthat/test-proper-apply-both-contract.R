library(npRmpi)

test_that("npcdens proper apply='both' combines fitted and predict slice behavior", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260324)
  x <- runif(40, -1, 1)
  y <- sin(2 * pi * x) + rnorm(40, sd = 0.18)
  nd <- rbind(
    data.frame(y = c(-0.75, -0.2, 0.4), x = rep(-0.35, 3L)),
    data.frame(y = c(-0.1, 0.3, 0.75, 1.0), x = rep(0.45, 4L))
  )

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  ctrl.both <- list(mode = "slice", apply = "both", slice.grid.size = 21L, slice.extend.factor = 0)
  ctrl.fit <- modifyList(ctrl.both, list(apply = "fitted"))

  fit.raw <- npcdens(bws = bw, txdat = data.frame(x = x), tydat = data.frame(y = y))
  fit.fitted <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE,
    proper.control = ctrl.fit
  )
  fit.both <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE,
    proper.control = ctrl.both
  )

  pred.base <- predict(fit.raw, newdata = nd, proper = FALSE)
  pred.ref <- predict(fit.raw, newdata = nd, proper = TRUE, proper.control = ctrl.both)
  pred.both <- predict(fit.both, newdata = nd, proper.control = ctrl.both)

  expect_true(isTRUE(fit.both$proper.requested))
  expect_true(isTRUE(fit.both$proper.applied))
  expect_identical(fit.both$proper.info$route, "slice")
  expect_equal(fit.both$condens.raw, fit.raw$condens, tolerance = 1e-12)
  expect_equal(fit.both$condens, fit.fitted$condens, tolerance = 1e-10)
  expect_equal(as.numeric(pred.both), as.numeric(pred.ref), tolerance = 1e-10)
  expect_true(any(abs(as.numeric(pred.both) - as.numeric(pred.base)) > 1e-8))
})

test_that("npcdist proper apply='both' combines fitted and predict slice behavior", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260324)
  x <- runif(40, -1, 1)
  y <- sin(2 * pi * x) + rnorm(40, sd = 0.18)
  nd <- rbind(
    data.frame(y = c(-0.75, -0.2, 0.4), x = rep(-0.35, 3L)),
    data.frame(y = c(-0.1, 0.3, 0.75, 1.0), x = rep(0.45, 4L))
  )

  bw <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  ctrl.both <- list(mode = "slice", apply = "both", slice.grid.size = 21L, slice.extend.factor = 0)
  ctrl.fit <- modifyList(ctrl.both, list(apply = "fitted"))

  fit.raw <- npcdist(bws = bw, txdat = data.frame(x = x), tydat = data.frame(y = y))
  fit.fitted <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE,
    proper.control = ctrl.fit
  )
  fit.both <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE,
    proper.control = ctrl.both
  )

  pred.base <- predict(fit.raw, newdata = nd, proper = FALSE)
  pred.ref <- predict(fit.raw, newdata = nd, proper = TRUE, proper.control = ctrl.both)
  pred.both <- predict(fit.both, newdata = nd, proper.control = ctrl.both)

  expect_true(isTRUE(fit.both$proper.requested))
  expect_true(isTRUE(fit.both$proper.applied))
  expect_identical(fit.both$proper.info$route, "slice")
  expect_equal(fit.both$condist.raw, fit.raw$condist, tolerance = 1e-12)
  expect_equal(fit.both$condist, fit.fitted$condist, tolerance = 1e-10)
  expect_equal(as.numeric(pred.both), as.numeric(pred.ref), tolerance = 1e-10)
})
