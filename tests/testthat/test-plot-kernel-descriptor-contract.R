make_plot_kernel_bandwidth <- function(xdat,
                                       bandwidth,
                                       ckertype,
                                       ckerorder = 2L,
                                       ckerbound = "none",
                                       ckerlb = NA_real_,
                                       ckerub = NA_real_) {
  npRmpi:::kbandwidth(
    bw = bandwidth,
    xdati = npRmpi:::untangle(xdat),
    xnames = names(xdat),
    bwscaling = FALSE,
    ckertype = ckertype,
    ckerorder = ckerorder,
    ckerbound = ckerbound,
    ckerlb = ckerlb,
    ckerub = ckerub
  )
}

test_that("direct plot weights use canonical kernel code and descriptor owners", {
  helper_text <- paste(
    deparse(body(npRmpi:::.np_plot_kernel_weights_direct)),
    collapse = "\n"
  )

  expect_match(
    helper_text,
    "kerneval = npContinuousKernelCode\\(bws\\)"
  )
  expect_match(
    helper_text,
    "npContinuousKernelDescriptorOptions\\(bws\\)"
  )
  expect_match(
    helper_text,
    "divide.returned.kernel.weights = FALSE",
    fixed = TRUE
  )
})

test_that("direct plot weights match independent order-two kernel formulas", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)

  training <- data.frame(x = c(0.05, 0.20, 0.48, 0.72, 0.94))
  evaluation <- data.frame(x = c(0.16, 0.51, 0.86))
  h <- 0.24
  standardized <- outer(training$x, evaluation$x, "-") / h
  helper <- npRmpi:::.np_plot_kernel_weights_direct

  expected <- list(
    gaussian = dnorm(standardized),
    epanechnikov = 3 / (4 * sqrt(5)) *
      (1 - standardized^2 / 5) *
      (abs(standardized) <= sqrt(5)),
    uniform = 0.5 * (abs(standardized) <= 1)
  )

  for (family in names(expected)) {
    bws <- suppressWarnings(make_plot_kernel_bandwidth(
      xdat = training,
      bandwidth = h,
      ckertype = family
    ))
    actual <- helper(
      bws = bws,
      txdat = training,
      exdat = evaluation,
      operator = "normal",
      bandwidth.divide = TRUE
    )
    expect_equal(actual, expected[[family]], tolerance = 2e-14)
  }

  training.beta <- data.frame(x = c(0, 0.08, 0.31, 0.67, 0.93, 1))
  evaluation.beta <- data.frame(x = c(0, 0.22, 0.59, 1))
  h <- 0.14
  tau <- (1 / h)^2
  bws <- make_plot_kernel_bandwidth(
    xdat = training.beta,
    bandwidth = h,
    ckertype = "beta",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  actual <- npRmpi:::.np_plot_kernel_weights_direct(
    bws = bws,
    txdat = training.beta,
    exdat = evaluation.beta,
    operator = "normal",
    bandwidth.divide = FALSE
  )
  expected.beta <- vapply(evaluation.beta$x, function(target) {
    dbeta(
      training.beta$x,
      shape1 = 1 + target * tau,
      shape2 = 1 + (1 - target) * tau
    )
  }, numeric(nrow(training.beta)))

  expect_equal(actual, expected.beta, tolerance = 2e-12)
})
