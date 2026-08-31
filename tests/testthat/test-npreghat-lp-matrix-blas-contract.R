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
      denom <- A.base[1L, 1L]
      if (!is.finite(denom) || denom == 0.0)
        stop("undefined pristine ridge intercept anchor")
      solved[1L] <- solved[1L] + nepsilon * (solved[1L] / denom)
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

  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "np"
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

test_that("canonical LP ridge uses the signed pristine intercept anchor", {
  solve_response <- function(anchor) {
    .Call(
      "C_np_lp_batch_project",
      matrix(c(anchor, 0, 0), nrow = 1L),
      matrix(c(3 * anchor, 0), nrow = 1L),
      c(1, 0), 2, TRUE, PACKAGE = "np"
    )
  }
  influence <- function(anchor) {
    .Call(
      "C_np_reghat_lp_matrix_fast",
      matrix(rep(anchor / 2, 2L), nrow = 2L),
      rbind(c(1, 0), c(1, 0)),
      matrix(c(1, 0), nrow = 1L), PACKAGE = "np"
    )
  }

  for (anchor in c(1, 2^-900, -1)) {
    response <- solve_response(anchor)
    hat <- influence(anchor)
    expect_identical(response$status, 0L)
    expect_equal(response$values[[1L]], 3, tolerance = 4e-15)
    expect_equal(sum(hat), 1, tolerance = 4e-15)
  }

  zero.response <- solve_response(0)
  expect_identical(zero.response$status, 1L)
  expect_null(zero.response$values)
  expect_error(influence(0), "invalid wider-LP compiled hat-matrix input")
})

test_that("width-one scalar hats retain signed higher-order kernel weights", {
  old <- getOption("matprod")
  on.exit(options(matprod = old), add = TRUE)

  set.seed(2026072504L)
  n <- 97L
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  ex <- tx[seq_len(19L), , drop = FALSE]
  y <- sin(tx$x1) - cos(tx$x2) + rnorm(n, sd = 0.1)
  bw <- npregbw(
    xdat = tx, ydat = y, regtype = "lp", degree = c(0L, 0L),
    degree.select = "manual", basis = "glp", bernstein.basis = FALSE,
    bwmethod = "cv.ls", bwtype = "fixed", ckertype = "gaussian",
    ckerorder = 4L, bws = c(0.24, 0.29), bandwidth.compute = FALSE
  )
  kw <- suppressWarnings(np:::.np_kernel_weights_direct(
    bws = bw, txdat = tx, exdat = ex, bandwidth.divide = TRUE,
    kernel.pow = 1.0,
    int.do.tree = np:::.npreg_fit_tree_code(
      bw, ncon = bw$ncon, ncat = bw$nuno + bw$nord
    )
  ))
  W.train <- np:::W.lp(
    xdat = tx, degree = c(0L, 0L), basis = "glp",
    bernstein.basis = FALSE
  )
  W.eval <- np:::W.lp(
    xdat = ex, degree = c(0L, 0L), basis = "glp",
    bernstein.basis = FALSE
  )

  expect_true(any(kw < 0.0))
  reference <- .np_test_lp_hat_matrix_reference(kw, W.train, W.eval)
  compiled <- .Call(
    "C_np_reghat_lp_matrix_fast",
    as.matrix(kw), as.matrix(W.train), as.matrix(W.eval), PACKAGE = "np"
  )
  expect_equal(compiled, reference, tolerance = 5e-15)

  options(matprod = "internal")
  internal <- npreghat(
    bws = bw, txdat = tx, exdat = ex, output = "matrix", s = 0L
  )
  options(matprod = "default")
  default <- npreghat(
    bws = bw, txdat = tx, exdat = ex, output = "matrix", s = 0L
  )
  expect_identical(dim(internal), dim(compiled))
  expect_identical(dim(default), dim(compiled))
  expect_identical(as.double(internal), as.double(default))
  expect_equal(as.double(default), as.double(compiled), tolerance = 5e-15)
})

test_that("width-one scalar hats own degenerate and non-finite systems", {
  scalar_hat <- function(weights, basis.eval = 1.0) {
    .Call(
      "C_np_reghat_lp_matrix_fast",
      matrix(as.double(weights), ncol = 1L),
      matrix(1.0, nrow = length(weights), ncol = 1L),
      matrix(as.double(basis.eval), nrow = 1L, ncol = 1L),
      PACKAGE = "np"
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
