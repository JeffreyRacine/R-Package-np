test_that("npregbw NOMAD degree search backend improves over the baseline", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(16)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  bw <- npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 2L
  )

  expect_s3_class(bw, "rbandwidth")
  expect_identical(bw$degree.search$mode, "nomad")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_gte(bw$degree.search$n.unique, 1L)
  expect_lte(bw$degree.search$best.fval, bw$degree.search$baseline.fval + 1e-10)
  expect_identical(length(bw$degree.search$restart.degree.starts), 2L)
  expect_identical(length(bw$degree.search$restart.results), 2L)
  expect_identical(bw$degree.search$restart.results[[1L]]$restart, 1L)
  expect_identical(length(bw$nomad.restart.fval), 2L)
  expect_identical(length(bw$nomad.restart.results), 2L)
  expect_identical(length(bw$nomad.restart.degree.starts), 2L)
  restart.objectives <- unlist(lapply(
    bw$degree.search$restart.results,
    function(x) if (is.null(x$objective)) NA_real_ else as.numeric(x$objective[1L])
  ))
  restart.objectives <- restart.objectives[is.finite(restart.objectives)]
  expect_true(length(restart.objectives) >= 1L)
  expect_lte(bw$degree.search$best.fval, min(restart.objectives) + 1e-6)
  expect_identical(as.numeric(bw$nomad.restart.fval), as.numeric(unlist(lapply(
    bw$nomad.restart.results,
    function(x) if (is.null(x$objective)) NA_real_ else as.numeric(x$objective[1L])
  ))))
  expect_identical(bw$nomad.best.restart, which.min(bw$nomad.restart.fval))
})
