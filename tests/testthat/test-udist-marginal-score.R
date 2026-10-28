test_that("CDF score envelope uses the signed range and preserves nesting", {
  for (q in 1:2) {
    range <- .np_udist_marginal_range(q)
    expect_lt(range[1], 0)
    expect_equal(unname(sum(range)), 1)
    center <- c(-.01, .5, 1.01)
    alpha <- c(.05, .1)
    widths <- list(pointwise = matrix(.01, 3, 2),
      simultaneous = matrix(.02, 3, 2), bonferroni = matrix(.03, 3, 2))
    bands <- lapply(widths, function(w) {
      z <- array(0, c(3, 2, 2)); z[, 1, ] <- center-w; z[, 2, ] <- center+w; z
    })
    answer <- list(center = center, intervals = bands, critical = widths, settings = list())
    z <- .np_udist_marginal_envelope(answer, widths, 100, q, alpha)
    expect_identical(z$center, center)
    expect_identical(z$reference.bands, bands)
    expect_true(all(z$intervals$pointwise[, 1, ] <= bands$pointwise[, 1, ]))
    expect_true(all(z$intervals$pointwise[, 2, ] >= bands$pointwise[, 2, ]))
    expect_true(all(z$intervals$simultaneous[, 1, ] <= z$intervals$pointwise[, 1, ]))
    expect_true(all(z$intervals$simultaneous[, 2, ] >= z$intervals$pointwise[, 2, ]))
    expect_true(all(z$intervals$bonferroni[, 1, ] <= z$intervals$simultaneous[, 1, ]))
    expect_true(all(z$intervals$bonferroni[, 2, ] >= z$intervals$simultaneous[, 2, ]))
    widths <- lapply(widths, function(w) w[1, , drop = FALSE])
    answer$center <- center[1]
    answer$intervals <- lapply(bands, function(b) b[1, , , drop = FALSE])
    answer$critical <- widths
    one <- .np_udist_marginal_envelope(answer, widths, 100, q, alpha)
    expect_identical(one$intervals$pointwise, one$intervals$simultaneous)
    expect_identical(one$intervals$pointwise, one$intervals$bonferroni)
  }
})
