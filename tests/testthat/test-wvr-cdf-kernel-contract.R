test_that("Wang-van Ryzin integrals equal cumulative PMF mass", {
  old <- options(np.messages = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  dat <- data.frame(z = ordered(1:5, levels = 1:5))
  for (lambda in c(0, 0.2, 0.73, 1 - 1e-8)) {
    pmf <- function(z, center) {
      ifelse(z == center, 1 - lambda,
             0.5 * (1 - lambda) * lambda^abs(z - center))
    }
    expected <- vapply(1:5, function(threshold) {
      vapply(1:5, function(center) {
        anchor <- 0.5 + 0.5 * pmf(center, center)
        if (threshold == center) return(anchor)
        if (threshold > center)
          return(anchor + sum(pmf(seq.int(center + 1L, threshold), center)))
        anchor - sum(pmf(seq.int(threshold + 1L, center), center))
      }, numeric(1L))
    }, numeric(5L))
    actual <- npksum(txdat = dat, exdat = dat, bws = lambda,
                    okertype = "wangvanryzin", operator = "integral",
                    return.kernel.weights = TRUE)$kw
    expect_equal(actual, expected, tolerance = 3e-14)
    # Near one, subtracting two CDFs near 1/2 needs an absolute error bound.
    expect_lt(max(abs(actual[, -1L] - actual[, -5L] -
      outer(1:5, 2:5, function(center, z) pmf(z, center)))), 3e-14)
  }
})
