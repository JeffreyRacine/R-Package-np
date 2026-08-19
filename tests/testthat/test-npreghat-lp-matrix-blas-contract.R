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

test_that("compiled LP hat matrix agrees with the public canonical owner", {
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
  kw <- suppressWarnings(npRmpi:::.np_kernel_weights_direct(
    bws = bw, txdat = tx, bandwidth.divide = TRUE, kernel.pow = 1.0,
    int.do.tree = npRmpi:::.npreg_fit_tree_code(
      bw, ncon = bw$ncon, ncat = bw$nuno + bw$nord
    )
  ))
  W.train <- npRmpi:::W.lp(xdat = tx, degree = c(1L, 1L), basis = "glp",
                           bernstein.basis = FALSE)
  W.eval <- npRmpi:::W.lp(xdat = tx, degree = c(1L, 1L), basis = "glp",
                          bernstein.basis = FALSE)

  reference <- .np_test_lp_hat_matrix_reference(kw, W.train, W.eval)
  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "npRmpi"
  )
  expect_equal(compiled, reference, tolerance = 1e-12)
  expect_equal(as.double(npreghat(bws = bw, txdat = tx)),
               as.double(reference), tolerance = 1e-12)
})

test_that("compiled and public LP hats share the canonical ridge policy", {
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
  kw <- suppressWarnings(npRmpi:::.np_kernel_weights_direct(
    bws = bw, txdat = tx, bandwidth.divide = TRUE, kernel.pow = 1.0,
    int.do.tree = npRmpi:::.npreg_fit_tree_code(
      bw, ncon = bw$ncon, ncat = bw$nuno + bw$nord
    )
  ))
  W.train <- npRmpi:::W.lp(xdat = tx, degree = c(2L, 2L), basis = "glp",
                           bernstein.basis = FALSE)
  W.eval <- npRmpi:::W.lp(xdat = tx, degree = c(2L, 2L), basis = "glp",
                          bernstein.basis = FALSE)

  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "npRmpi"
  )
  public <- npreghat(bws = bw, txdat = tx, output = "matrix")
  fit <- npreg(
    bws = bw, txdat = tx, tydat = y,
    warn.glp.gradient = FALSE
  )$mean

  expect_gt(sum(attr(public, "ridge.used", exact = TRUE) > 0), 0L)
  expect_equal(as.vector(compiled), as.vector(public), tolerance = 1e-10)
  expect_equal(as.vector(compiled %*% y), as.vector(fit), tolerance = 1e-8)
})

test_that("width-one scalar hats retain signed higher-order kernel weights", {
  set.seed(2026072504L)
  n <- 97L
  neval <- 19L
  x <- runif(n)
  xe <- seq(0.05, 0.95, length.out = neval)
  distance <- outer(x, xe, "-") / 0.19
  kw <- exp(-0.5 * distance^2) * (1.5 - distance^2)
  W.train <- matrix(1.0, nrow = n, ncol = 1L)
  W.eval <- matrix(1.0, nrow = neval, ncol = 1L)

  expect_true(any(kw < 0.0))
  reference <- .np_test_lp_hat_matrix_reference(kw, W.train, W.eval)
  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "npRmpi"
  )
  expect_equal(compiled, reference, tolerance = 5e-15)
})

test_that("width-one scalar hats own degenerate and non-finite systems", {
  scalar_hat <- function(weights, basis.eval = 1.0) {
    .Call(
      "C_np_reghat_lp_matrix_fast",
      matrix(as.double(weights), ncol = 1L),
      matrix(1.0, nrow = length(weights), ncol = 1L),
      matrix(as.double(basis.eval), nrow = 1L, ncol = 1L),
      PACKAGE = "npRmpi"
    )
  }

  expect_identical(scalar_hat(c(0.0, 0.0)), matrix(c(0.0, 0.0), 1L))
  expect_identical(scalar_hat(c(1.0, -1.0)), matrix(c(4.0, -4.0), 1L))
  expect_identical(
    scalar_hat(c(1.0, -1.0 + 2^-50)),
    matrix(c(2^50, -2^50 + 1.0), 1L)
  )
  expect_error(
    scalar_hat(c(1.0, NaN)),
    "non-finite system",
    fixed = TRUE
  )
  expect_error(
    scalar_hat(c(1.0, 1.0), basis.eval = Inf),
    "non-finite system",
    fixed = TRUE
  )
})
