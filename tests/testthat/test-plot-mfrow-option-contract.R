test_that("plot.par.mfrow=FALSE is honored for npcdens auto and manual plots", {
  skip_if_not(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old.opts <- options(plot.par.mfrow = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(101)
  n <- 40
  x <- runif(n)
  y <- rbeta(n, 1, 1)

  f.auto <- npcdens(
    y ~ x,
    cykerbound = "range",
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "cell",
    degree.max = 1,
    degree.verify = FALSE,
    bwmethod = "cv.ml",
    nmulti = 1
  )

  f.manual <- npcdens(
    y ~ x,
    cykerbound = "range",
    regtype = "lp",
    degree = 1,
    bwmethod = "cv.ml",
    nmulti = 1
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  par(mfrow = c(2, 2))
  plot(1:10)
  plot(1:10)
  suppressWarnings(
    plot(f.auto, theta = 70, phi = 10, view = "fixed", zlim = c(0, 2))
  )
  expect_identical(par("mfrow"), c(2L, 2L))
  expect_identical(par("mfg"), c(2L, 1L, 2L, 2L))

  suppressWarnings(
    plot(f.manual, theta = 70, phi = 10, view = "fixed", zlim = c(0, 2))
  )
  expect_identical(par("mfrow"), c(2L, 2L))
  expect_identical(par("mfg"), c(2L, 2L, 2L, 2L))
})
