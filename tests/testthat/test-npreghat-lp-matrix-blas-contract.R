.np_test_lp_hat_matrix_reference <- function(kw, W.train, W.eval) {
  ntrain <- nrow(W.train)
  neval <- nrow(W.eval)
  H <- matrix(NA_real_, nrow = neval, ncol = ntrain)
  eps <- 1.0 / max(1L, ntrain)

  for (j in seq_len(neval)) {
    w <- kw[, j]
    A.base <- crossprod(W.train, W.train * w)
    rhs <- W.eval[j, ]
    solved <- tryCatch(solve(A.base, rhs), error = function(e) NULL)

    if (is.null(solved) || !all(is.finite(solved))) {
      A.try <- A.base
      nepsilon <- 0.0
      repeat {
        diag(A.try) <- diag(A.try) + eps
        nepsilon <- nepsilon + eps
        solved <- tryCatch(solve(A.try, rhs), error = function(e) NULL)
        if (!is.null(solved) && all(is.finite(solved)))
          break
      }
      denom <- A.try[1L, 1L]
      if (!is.finite(denom) || abs(denom) < .Machine$double.xmin)
        denom <- .Machine$double.xmin
      solved[1L] <- solved[1L] * (1.0 + nepsilon / denom)
    }

    H[j, ] <- w * drop(W.train %*% solved)
  }
  H
}

test_that("compiled LP hat matrix preserves the BLAS-backed R loop exactly", {
  old <- getOption("matprod")
  on.exit(options(matprod = old), add = TRUE)
  options(matprod = "default")

  set.seed(2026072299L)
  n <- 97L
  x1 <- runif(n)
  x2 <- runif(n)
  W.train <- cbind(1.0, x1, x2)
  W.eval <- W.train
  d1 <- outer(x1, x1, "-") / 0.24
  d2 <- outer(x2, x2, "-") / 0.29
  kw <- exp(-0.5 * (d1^2 + d2^2))

  reference <- .np_test_lp_hat_matrix_reference(kw, W.train, W.eval)
  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "npRmpi"
  )
  expect_identical(compiled, reference)
})

test_that("compiled LP hat matrix preserves the incumbent ridge sequence", {
  old <- getOption("matprod")
  on.exit(options(matprod = old), add = TRUE)
  options(matprod = "default")

  set.seed(2026072300L)
  n <- 41L
  x1 <- runif(n)
  x2 <- runif(n)
  W.train <- cbind(1.0, x1, x1^2, x2, x2^2, x1 * x2)
  W.eval <- W.train
  kw <- matrix(0.0, nrow = n, ncol = n)
  diag(kw) <- 1.0

  reference <- .np_test_lp_hat_matrix_reference(kw, W.train, W.eval)
  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "npRmpi"
  )
  expect_identical(compiled, reference)
})
