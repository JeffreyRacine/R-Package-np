test_that("private D8 rows use the complete corrected observation covariance", {
  x <- data.frame(x = c(-1.8, -.4, .1, .7, 1.4))
  grid <- data.frame(x = c(-.9, .2, 1.1))
  state <- list(x = x, n = nrow(x), q = 1L, r = 0L, s = 0L, h = .6)
  actual <- .np_udens_d8_rows(state, grid)
  u <- outer(grid$x, x$x, "-") / .6
  literal <- dnorm(u) / .6 * (35/16 - 35*u^2/16 + 7*u^4/16 - u^6/48)
  expected <- rowMeans(literal)
  count.covariance <- diag(nrow(x)) - matrix(1/nrow(x), nrow(x), nrow(x))
  expect_equal(actual$center, expected, tolerance = 1e-13)
  expect_equal(tcrossprod(actual$factor),
    literal %*% count.covariance %*% t(literal) / nrow(x)^2,
    tolerance = 1e-13)
  expect_equal(actual$se^2, rowSums(actual$factor^2), tolerance = 1e-13)
  expect_equal(rowSums(actual$factor), rep(0, nrow(grid)), tolerance = 1e-13)
})

test_that("D8 has total-degree radial moments and coordinate symmetry", {
  for (q in 1:3) {
    a <- 1 + q/2
    coefficients <- c(a*(a+1)*(a+2)/6, -(a+1)*(a+2)/2, (a+2)/2, -1/6)
    for (k in 0:3)
      expect_equal(sum(coefficients * exp(lgamma(q/2+k+0:3)-lgamma(q/2))),
        as.numeric(k == 0), tolerance = 1e-10)
  }
  x <- data.frame(x = c(-1, .3, 1), y = c(.2, -1, .6))
  grid <- data.frame(x = c(-.1, .5), y = c(.4, -.2))
  state <- list(x = x, n = 3L, q = 2L, r = 0L, s = 0L, h = c(.5, .8))
  actual <- .np_udens_d8_rows(state, grid)
  state$x <- x[2:1]; state$h <- state$h[2:1]
  expect_equal(.np_udens_d8_rows(state, grid[2:1]), actual, tolerance = 1e-13)
})

test_that("D8 finite extreme distances do not produce zero-times-infinity", {
  state <- list(x = data.frame(x = c(-1, 0, 1)), n = 3L, q = 1L,
                r = 0L, s = 0L, h = .5)
  actual <- .np_udens_d8_rows(state, data.frame(x = c(-1e300, 1e300)))
  expect_identical(actual$center, c(0, 0))
  expect_true(all(actual$factor == 0))
  expect_identical(actual$se, c(0, 0))
})
