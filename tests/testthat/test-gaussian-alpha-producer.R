# Small producer contract only; coordinator/topology matrices live in the
# retained forensic campaign. This does not assert general seed-only replay.
test_that("Gaussian alpha reductions share one uninterrupted normal stream", {
  kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(kind)), add = TRUE)
  quiet <- function(code) withCallingHandlers(force(code), warning = function(w) {
    if (grepl("below the recommended Monte Carlo precision count",
              conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
  alpha <- c(.05, .1)
  for (normal in c("Inversion", "Box-Muller")) for (pending in 0:1) {
    reset <- function() {
      RNGkind("Mersenne-Twister", normal, "Rejection")
      set.seed(6108)
      if (pending) invisible(rnorm(1))
    }
    reset()
    maxima <- apply(abs(matrix(rnorm(3L * 13L), 3L)), 2L, max)
    raw <- unname(quantile(maxima, 1 - alpha, type = 1L))
    continuation <- list(seed = .Random.seed, normal = rnorm(5), uniform = runif(3))
    reset()
    out <- quiet(.np_plot_gaussian_calibrate_alphas(diag(3), alpha,
                                                    draws = 13L, block = 5L))
    expect_identical(vapply(out, `[[`, 0, "critical.raw"), raw)
    expect_identical(list(seed = .Random.seed, normal = rnorm(5), uniform = runif(3)),
                     continuation)
    reset()
    scalar <- quiet(.np_plot_gaussian_calibrate(diag(3), .05,
                                                draws = 13L, block = 5L))
    expect_identical(scalar, out[[1L]])
  }
  expect_error(.np_plot_gaussian_calibrate(diag(2), alpha = alpha), "alpha")
  expect_error(.np_plot_gaussian_calibrate_alphas(diag(2), c(.1, .1)), "unique")
})

test_that("signed Gaussian maxima use one uninterrupted normal sample", {
  kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(kind)), add = TRUE)
  quiet <- function(code) withCallingHandlers(force(code), warning = function(w) {
    if (grepl("below the recommended Monte Carlo precision count",
              conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
  for (normal in c("Inversion", "Box-Muller")) for (pending in 0:1) {
    reset <- function() {
      RNGkind("Mersenne-Twister", normal, "Rejection"); set.seed(314)
      if (pending) invisible(rnorm(1))
    }
    for (m in 2:3) {
      reset()
      expected <- apply(-matrix(rnorm(m * 13L), m), 2L, max)
      continuation <- list(seed = .Random.seed, normal = rnorm(5), uniform = runif(3))
      reset()
      out <- quiet(.np_udens_scale_gaussian(diag(m), 13L, 1024^2))
      expect_identical(out$maxima, expected)
      expect_identical(list(seed = .Random.seed, normal = rnorm(5), uniform = runif(3)),
                       continuation)
    }
    reset()
    expected <- -colSums(matrix(rnorm(26L), 2L)) / sqrt(2)
    reset()
    out <- quiet(.np_udens_scale_gaussian(matrix(1, 2L, 1L), 13L, 1024^2))
    expect_equal(out$maxima, expected, tolerance = 1e-13)
    expect_identical(out$rank, 1L)
  }
})
