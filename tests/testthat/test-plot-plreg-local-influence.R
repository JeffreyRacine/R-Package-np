test_that("partially linear plot influence keeps LC row blocks on the master", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(43032)
  n <- 120L
  z <- data.frame(z = runif(n))
  x <- data.frame(x = rnorm(n) + z$z)
  y <- x$x + sin(2 * pi * z$z) + rnorm(n, sd = .5)
  bw <- npplregbw(xdat = x, ydat = y, zdat = z,
    bws = matrix(c(.15, .2), 2L, 1L), bandwidth.compute = FALSE)
  influence <- getFromNamespace(".np_plot_plreg_influence", "npRmpi")
  called <- 0L
  local_mocked_bindings(.np_plot_plreg_influence = function(...) {
    called <<- called + 1L
    expect_true(isTRUE(getOption("npRmpi.local.regression.mode")))
    expect_true(isTRUE(getOption("npRmpi.autodispatch.disable")))
    influence(...)
  }, .package = "npRmpi")
  state.names <- c("npRmpi.local.regression.mode", "npRmpi.autodispatch.disable",
                   "npRmpi.autodispatch.context")
  before <- options()[state.names]
  out <- plot(bw, errors = "asymptotic", output = "data", neval = 7L)
  expect_gt(called, 0L)
  expect_true(all(is.finite(unlist(lapply(out, function(panel) panel$merr)))))
  expect_identical(options()[state.names], before)
  local_mocked_bindings(.np_plot_plreg_influence = function(...) {
    stop("intentional influence failure")
  }, .package = "npRmpi")
  expect_error(plot(bw, errors = "asymptotic", output = "data", neval = 7L),
               "intentional influence failure")
  expect_identical(options()[state.names], before)
})
