test_that("npcdensbw NOMAD degree search backend improves over the baseline", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(14)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  bw <- npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(bw, "conbandwidth")
  expect_identical(bw$degree.search$mode, "nomad")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_gte(bw$degree.search$n.unique, 1L)
  expect_gte(bw$degree.search$best.fval, bw$degree.search$baseline.fval - 1e-10)
})

test_that("npcdensbw beta degree search uses the canonical LP shadow", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(2026080114L)
  n <- 14L
  x <- data.frame(x = sort(runif(n, 0.04, 0.96)))
  y <- data.frame(y = pmin(
    0.96, pmax(0.04, 0.25 + 0.55 * x$x + rnorm(n, sd = 0.06))
  ))

  bw <- npRmpi::npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L,
    cxkertype = "beta", cxkerorder = 2L,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykertype = "beta", cykerorder = 2L,
    cykerbound = "fixed", cykerlb = 0, cykerub = 1
  )

  expect_s3_class(bw, "conbandwidth")
  expect_identical(bw$degree.search$mode, "nomad")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_true(all(is.finite(c(bw$xbw, bw$ybw, bw$fval[1L]))))
})

test_that("npcdensbw automatic degree search defaults to NOMAD plus Powell", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(14)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  bw <- npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_true(is.finite(bw$nomad.time))
  expect_true(is.finite(bw$powell.time))
  expect_equal(as.double(bw$total.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
  expect_equal(as.double(bw$degree.search$optim.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
})
