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
