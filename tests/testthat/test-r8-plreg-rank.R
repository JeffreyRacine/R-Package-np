test_that("partial-linear rank validation reuses QR and respects formation scale", {
  ns <- asNamespace("np")
  check <- get(".np_plreg_check_residualized_rank", ns)
  formation <- get(".np_plreg_residual_formation_error", ns)
  set.seed(681)
  X <- matrix(rnorm(96L), 32L, 3L)
  bad <- cbind(X[, 1L], X[, 1L])
  for (perm in list(1:2, 2:1)) {
    expect_error(check(qr(bad[, perm], tol = .Machine$double.eps), 2L, "test"),
                 "rank deficient after smoothing", fixed = TRUE)
  }
  expect_error(check(qr(cbind(X[, 1L], 0), tol = .Machine$double.eps), 2L, "test"),
               "rank deficient after smoothing", fixed = TRUE)
  expect_error(check(qr(cbind(X[, 1:2], X[, 1L]+2*X[, 2L]),
                           tol = .Machine$double.eps), 3L, "test"),
               "rank deficient after smoothing", fixed = TRUE)
  for (s in list(c(1, 1, 1), c(.001, 1, 1000))) {
    q <- qr(sweep(X, 2L, s, `*`), tol = .Machine$double.eps)
    saved <- q
    expect_silent(check(q, 3L, "test"))
    expect_identical(q, saved)
  }
  # Direct factor boundary: no BLAS-dependent QR ambiguity at the threshold.
  t <- 32*.Machine$double.eps
  for (a in c(.5, 1, 2)) {
    R <- matrix(0, 32L, 2L); R[1L, ] <- 1; R[2L, 2L] <- a*t
    q <- list(qr = R, rank = 2L, pivot = 1:2)
    if (a <= 1) expect_error(check(q, 2L, "test"), "rank deficient", fixed = TRUE)
    else expect_silent(check(q, 2L, "test"))
  }
  x <- rep(1, 32L); xhat <- x - .Machine$double.eps
  e <- formation(x, xhat)
  expect_error(check(qr(matrix(x-xhat)), 1L, "test", e),
               "rank deficient", fixed = TRUE)
  expect_silent(check(qr(matrix(rep(1e-10, 32L))), 1L, "test", e))
})
