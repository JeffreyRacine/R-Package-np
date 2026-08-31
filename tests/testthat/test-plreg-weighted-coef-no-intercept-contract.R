test_that("partial-linear weighted ridge has no privileged coordinate", {
  solve.coef <- getFromNamespace(".np_plreg_weighted_coef", "np")
  X <- cbind(
    x1 = c(-2, -1, 0, 1, 2, 3),
    x2 = c(1, 0, 2, -1, 3, 1)
  )
  y <- c(-1.0, 0.2, 1.7, -0.4, 3.1, 2.2)
  w <- c(2, 0, 3, 1, 4, 2)
  ridge <- 0.25

  gram <- crossprod(X, X * w)
  rhs <- drop(crossprod(X, y * w))
  oracle <- unname(drop(solve(gram + diag(ridge, ncol(X)), rhs)))
  actual <- solve.coef(X = X, y = y, w = w, ridge = ridge)
  expect_equal(actual, oracle, tolerance = 2e-15)

  permutation <- c(2L, 1L)
  permuted <- solve.coef(
    X = X[, permutation, drop = FALSE], y = y, w = w, ridge = ridge
  )
  expect_equal(permuted[permutation], actual, tolerance = 2e-15)
})

test_that("partial-linear weighted solve starts unregularized and then uses 1/n", {
  solve.coef <- getFromNamespace(".np_plreg_weighted_coef", "np")
  X <- cbind(
    x1 = c(-2, -1, 0, 1, 2, 3),
    x2 = c(1, 0, 2, -1, 3, 1)
  )
  y <- c(-1.0, 0.2, 1.7, -0.4, 3.1, 2.2)
  w <- c(2, 0, 3, 1, 4, 2)
  oracle <- unname(drop(solve(crossprod(X, X * w), crossprod(X, y * w))))
  expect_equal(solve.coef(X, y, w, ridge = 0), oracle, tolerance = 2e-15)
  expect_equal(solve.coef(X, y, w), oracle, tolerance = 2e-15)

  singular.X <- cbind(x1 = X[, 1L], x2 = 2 * X[, 1L])
  ridge <- 1 / nrow(singular.X)
  singular.oracle <- unname(drop(solve(
    crossprod(singular.X, singular.X * w) + diag(ridge, 2L),
    crossprod(singular.X, y * w)
  )))
  singular.actual <- solve.coef(singular.X, y, w)
  expect_true(all(is.finite(singular.actual)))
  expect_equal(singular.actual, singular.oracle, tolerance = 2e-14)
})
