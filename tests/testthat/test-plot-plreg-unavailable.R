test_that("unreliable partially linear rows preserve curves and consistent band availability", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(101)
  n <- 60L
  z <- data.frame(z = runif(n))
  x <- data.frame(x = rnorm(n) + z$z)
  y <- x$x + sin(2 * pi * z$z) + rnorm(n, sd = .4)
  bw <- npplregbw(xdat = x, ydat = y, zdat = z,
    bws = matrix(c(.2, .25), 2L), bandwidth.compute = FALSE)
  reference <- plot(bw, errors = "asymptotic", band = "all",
                    perspective = FALSE, neval = 7L, output = "data")
  influence <- getFromNamespace(".np_plot_plreg_influence", "npRmpi")
  local_mocked_bindings(.np_plot_plreg_influence = function(bws, xdat, ydat,
      zdat, exdat, ezdat, ...) {
    out <- influence(bws, xdat, ydat, zdat, exdat, ezdat, ...)
    # Force the documented discrepancy only on the varying-z panel. This
    # tests the availability contract independently of platform-specific ridging.
    if (length(unique(ezdat[[1L]])) > 1L)
      out$mean[2L] <- out$mean[2L] + out$se[2L]
    out
  }, .package = "npRmpi")
  warnings <- character()
  out <- withCallingHandlers(plot(bw, errors = "asymptotic", band = "all",
      perspective = FALSE, neval = 7L, output = "data"), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("numerically unstable at 1 of 7", warnings, fixed = TRUE)))
  expect_identical(out[[1L]]$band.errors, reference[[1L]]$band.errors)
  for (j in seq_along(out)) expect_identical(out[[j]]$mean, reference[[j]]$mean)
  bands <- out[[2L]]$band.errors
  expect_true(all(is.na(bands$simultaneous)))
  for (family in c("pointwise", "bonferroni")) {
    expect_true(all(is.na(bands[[family]][2L, ])))
    expect_identical(bands[[family]][-2L, ], reference[[2L]]$band.errors[[family]][-2L, ])
  }
  for (band in c("pointwise", "bonferroni", "simultaneous")) {
    single <- suppressWarnings(plot(bw, errors = "asymptotic", band = band,
                                    perspective = FALSE, neval = 7L, output = "data"))
    expect_equal(single[[2L]]$merr, bands[[band]], ignore_attr = TRUE)
    expect_identical(single[[2L]]$mean, reference[[2L]]$mean)
  }
})
