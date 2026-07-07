test_that("npregbw automatic degree search defaults to NOMAD plus Powell", {
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

test_that("npregbw native NOMAD route has an explicit crs availability guard", {
  require_body <- paste(
    deparse(body(get(".npregbw_nomad_native_require_crs", envir = asNamespace("npRmpi"), inherits = FALSE))),
    collapse = "\n"
  )
  runner_body <- paste(
    deparse(body(get(".npregbw_run_fixed_degree_mads", envir = asNamespace("npRmpi"), inherits = FALSE))),
    collapse = "\n"
  )

  expect_match(require_body, "requireNamespace\\(\"crs\",\\s*quietly = TRUE\\)")
  expect_true(grepl("packageVersion(\"crs\") < \"0.15.44\"", require_body, fixed = TRUE))
  expect_true(grepl("native npreg NOMAD route requires crs >= 0.15-44", require_body, fixed = TRUE))
  expect_true(grepl(".npregbw_nomad_native_require_crs()", runner_body, fixed = TRUE))
})
