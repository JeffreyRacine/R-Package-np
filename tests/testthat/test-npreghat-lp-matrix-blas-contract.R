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
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(tx$x1) - cos(tx$x2) + rnorm(n, sd = 0.1)
  bw <- npregbw(
    xdat = tx, ydat = y, regtype = "lp", degree = c(1L, 1L),
    degree.select = "manual", basis = "glp", bernstein.basis = FALSE,
    bwmethod = "cv.ls", bwtype = "fixed", ckertype = "gaussian",
    ckerorder = 2L, bws = c(0.24, 0.29), bandwidth.compute = FALSE
  )
  kw <- suppressWarnings(np:::.np_kernel_weights_direct(
    bws = bw, txdat = tx, bandwidth.divide = TRUE, kernel.pow = 1.0,
    int.do.tree = np:::.npreg_fit_tree_code(
      bw, ncon = bw$ncon, ncat = bw$nuno + bw$nord
    )
  ))
  W.train <- np:::W.lp(xdat = tx, degree = c(1L, 1L), basis = "glp",
                       bernstein.basis = FALSE)
  W.eval <- np:::W.lp(xdat = tx, degree = c(1L, 1L), basis = "glp",
                      bernstein.basis = FALSE)

  reference <- .np_test_lp_hat_matrix_reference(kw, W.train, W.eval)
  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "np"
  )
  expect_identical(compiled, reference)
  expect_identical(as.double(npreghat(bws = bw, txdat = tx)),
                   as.double(reference))
})

test_that("compiled LP hat matrix preserves the incumbent ridge sequence", {
  old <- getOption("matprod")
  on.exit(options(matprod = old), add = TRUE)
  options(matprod = "default")

  set.seed(2026072300L)
  n <- 41L
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(tx$x1) + rnorm(n, sd = 0.1)
  bw <- suppressWarnings(npregbw(
    xdat = tx, ydat = y, regtype = "lp", degree = c(2L, 2L),
    degree.select = "manual", basis = "glp", bernstein.basis = FALSE,
    bwmethod = "cv.ls", bwtype = "fixed", ckertype = "uniform",
    ckerorder = 2L, bws = c(0.002, 0.002), bandwidth.compute = FALSE
  ))
  kw <- suppressWarnings(np:::.np_kernel_weights_direct(
    bws = bw, txdat = tx, bandwidth.divide = TRUE, kernel.pow = 1.0,
    int.do.tree = np:::.npreg_fit_tree_code(
      bw, ncon = bw$ncon, ncat = bw$nuno + bw$nord
    )
  ))
  W.train <- np:::W.lp(xdat = tx, degree = c(2L, 2L), basis = "glp",
                       bernstein.basis = FALSE)
  W.eval <- np:::W.lp(xdat = tx, degree = c(2L, 2L), basis = "glp",
                      bernstein.basis = FALSE)

  reference <- .np_test_lp_hat_matrix_reference(kw, W.train, W.eval)
  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "np"
  )
  expect_identical(compiled, reference)
})
