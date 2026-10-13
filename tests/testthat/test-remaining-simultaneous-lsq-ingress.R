test_that("LSQ plot rejects duplicate tau before fitting or seed changes", {
  # This ingress fixture deliberately has no usable fit: invalid tau must be
  # rejected before the plotting engine can inspect estimator state.
  object <- structure(list(tau = c(.3, .7), bws = list()),
                      class = "lsqregression")
  calls <- 0L
  local_mocked_bindings(.np_plot_condbandwidth_engine = function(...) {
    calls <<- calls + 1L
    stop("unexpected plot engine entry")
  }, .package = "npRmpi")
  set.seed(173)
  seed <- .Random.seed
  for (gradient in c(FALSE, TRUE)) {
    for (tau in list(c(.3, .3), c(.7, .3, .7))) {
      expect_error(plot(object, tau = tau, gradients = gradient, output = "data"),
                   "duplicate 'tau' values", fixed = TRUE)
    }
  }
  expect_identical(calls, 0L)
  expect_identical(.Random.seed, seed)
})
