test_that("the uniform Gaussian maximum respects shared innovations across targets", {
  F <- matrix(c(1, .7, -.2, .3, -.8, .6), nrow = 3)
  B <- 271L
  # The owner receives a square covariance factor, including any zero columns.
  factor <- cbind(F, 0)
  expected <- .np_with_seed(121L, {
    innovations <- matrix(rnorm(nrow(factor)*B), nrow(factor))
    apply(abs(factor %*% innovations), 2L, max)
  })
  actual <- .np_with_seed(121L, .np_udens_gaussian_absolute_maximum(factor, B, function(...) NULL))
  expect_equal(actual, expected, tolerance = 1e-14)
})
test_that("the joint envelope compares two marginal maximum quantiles", {
  a <- diag(c(1, .1))
  b <- diag(c(.1, 1))
  B <- 501L
  maxima <- .np_with_seed(27L, .np_udens_gaussian_absolute_maximum(a, B, function(...) NULL, b))
  z <- .np_with_seed(27L, matrix(rnorm(2L * B), 2L))
  reference <- cbind(apply(abs(a %*% z), 2L, max), apply(abs(b %*% z), 2L, max))
  expect_equal(unname(maxima), reference, tolerance = 1e-14)
  expect_lt(max(apply(maxima, 2L, quantile, .95)), unname(quantile(apply(maxima, 1L, max), .95)))
  g <- rep(seq(-1, 1, length.out = 11), each = 2)
  factor <- .np_udens_uniform_empirical_factor(seq(-2, 2, length.out = 9), .43, g)
  expect_equal(factor, t(factor), tolerance = 1e-12)
  expect_equal(factor[1, ], factor[2, ], tolerance = 1e-12)
})
