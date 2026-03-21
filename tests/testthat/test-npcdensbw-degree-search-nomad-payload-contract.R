test_that("npcdensbw direct nomad payload preserves CV metadata", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260320)
  dat <- data.frame(x = runif(24))
  dat$y <- rbeta(nrow(dat), 1, 1)

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
    nmulti = 1L,
    cxkerbound = "range",
    cykerbound = "range"
  )

  expect_identical(bw$method, "cv.ls")
  expect_identical(bw$pmethod, "Least Squares Cross-Validation")
  expect_true(length(bw$ifval) == 1L && is.na(bw$ifval))
  expect_true(length(bw$fval.history) == 1L && is.na(bw$fval.history))
  expect_true(length(bw$eval.history) == 1L && is.na(bw$eval.history))
  expect_true(length(bw$invalid.history) == 1L && is.na(bw$invalid.history))
  expect_equal(
    bw$fval,
    npRmpi:::.npcdensbw_eval_only(data.frame(x = dat$x), data.frame(y = dat$y), bw)$objective,
    tolerance = 1e-12
  )

  printed <- paste(capture.output(npRmpi:::print.conbandwidth(bw)), collapse = "\n")
  expect_false(grepl("Manual", printed, fixed = TRUE))
  expect_match(printed, "achieved on multistart 1", fixed = TRUE)
})
