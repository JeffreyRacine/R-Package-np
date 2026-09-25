test_that("conditioning advice retains existing warning and refusal severity", {
  set.seed(5)
  z <- runif(60, -1, 1)
  check <- function(x, degree = 2L, bernstein = FALSE) {
    npCheckRegressionDesignCondition(REGTYPE_LP, data.frame(x = x),
      degree = degree, bernstein.basis = bernstein)
  }
  expect_silent(check(z))
  expect_silent(check(1985 + z, degree = 3L, bernstein = TRUE))
  expect_warning(check(1985 + 35*z), "bernstein.basis=TRUE", fixed = TRUE)
  expect_error(check(1985 + z), "severely ill-conditioned")
  expect_error(check(10000 + z), "numerically rank deficient")
  expect_error(check(10000 + z), "centring/scaling", fixed = TRUE)
  expect_error(check(rep(1, 60), bernstein = TRUE), "rank deficient")
})
