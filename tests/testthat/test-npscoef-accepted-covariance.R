.dc2_covariance_test_helper <- function() {
  owner <- getFromNamespace(".np_scoef_fit_internal", "npRmpi")
  definitions <- Filter(function(expr)
    is.call(expr) && identical(expr[[1L]], quote(`<-`)) &&
      identical(expr[[2L]], quote(accepted_moment_covariance)),
    as.list(body(owner)))
  stopifnot(length(definitions) == 1L)
  eval(definitions[[1L]][[3L]], envir = environment(owner))
}

test_that("smooth-coefficient covariance preserves unregularized SPD arithmetic", {
  covariance <- .dc2_covariance_test_helper()
  a <- matrix(c(2, .4, .4, 1), 2L)
  s <- matrix(c(.8, .2, .2, .5), 2L)
  saved.a <- a
  saved.s <- s
  cm <- chol2inv(chol(a + 0 * diag(rep(1.0, 2L))))
  expect_identical(covariance(a, 0, s), cm %*% s %*% cm)
  expect_identical(a, saved.a)
  expect_identical(s, saved.s)
})

test_that("smooth-coefficient covariance includes both accepted RHS corrections", {
  covariance <- .dc2_covariance_test_helper()
  a <- matrix(c(2, .8, .8, 1.5), 2L)
  s <- matrix(c(.7, .2, .2, .9), 2L)
  saved.s <- s
  ridge <- .6
  cmat <- solve(a + diag(ridge, 2L))
  dmat <- diag(c(1 + ridge / a[1L, 1L], 1))
  g <- cmat %*% dmat
  expect_false(isTRUE(all.equal(g, t(g), tolerance = 1e-15)))
  expect_equal(covariance(a, ridge, s), g %*% s %*% t(g), tolerance = 1e-14)
  expect_identical(s, saved.s)
  expect_gt(abs((g %*% s %*% t(g))[1L, 2L]), .001)
  for (r in c(0, .25)) {
    expect_equal(covariance(matrix(2, 1L), r, matrix(.8, 1L)),
                 matrix(.8 / 4, 1L), tolerance = 1e-14)
  }
})

test_that("accepted indefinite covariance changes factorization not ridge", {
  covariance <- .dc2_covariance_test_helper()
  a <- matrix(c(2, .3, .3, -.5), 2L)
  s <- matrix(c(.8, .1, .1, .5), 2L)
  for (r in c(0, .1)) {
    g <- solve(a + diag(r, 2L)) %*% diag(c(1 + r / a[1L, 1L], 1))
    expected <- g %*% s %*% t(g)
    expect_equal(covariance(a, r, s), expected, tolerance = 1e-14)
    expect_gt(min(eigen(expected, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
  expect_null(covariance(matrix(1, 2L, 2L), 0, diag(2L)))
  # Covariance retains solve.default's condition check, not tol=0.
  expect_null(covariance(diag(c(1, -.Machine$double.eps^2)), 0, diag(2L)))
})

test_that("public signed-kernel smooth-coefficient errors use the point map", {
  skip_if_not(isTRUE(getOption("npRmpi.mpi.initialized", FALSE)),
              "public proof requires the test harness MPI pool")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(-1, 0, 0, 0, 0, 1))
  z <- data.frame(z = c(-2, -.2, -.1, .1, .2, 2))
  y <- c(.4, -.3, .8, -.2, .1, .7)
  b <- npscoefbw(xdat = x, ydat = y, zdat = z, bws = 1,
                 bandwidth.compute = FALSE, regtype = "lc",
                 ckertype = "gaussian", ckerorder = 4)
  fit <- npscoef(bws = b, txdat = x, tydat = y, tzdat = z,
                 exdat = data.frame(x = .5), ezdat = data.frame(z = 0),
                 se = TRUE, betas = TRUE, iterate = FALSE)
  design <- cbind(1, x$x)
  kernel <- function(u) (1.5 - .5 * u^2) * dnorm(u)
  train.weights <- outer(z$z, z$z, function(a, b) kernel(a - b))
  train.beta <- vapply(seq_along(y), function(j)
    solve(crossprod(design, design * train.weights[, j]),
          crossprod(design, y * train.weights[, j])), numeric(2L))
  residual <- y - rowSums(design * t(train.beta))
  w <- kernel(z$z)
  a <- crossprod(design, design * w)
  score.meat <- crossprod(design, design * (w * residual)^2)
  cm <- solve(a)
  expected <- cm %*% score.meat %*% t(cm)
  point <- solve(a, crossprod(design, y * w))
  expect_equal(as.double(fit$beta), as.double(point), tolerance = 2e-13)
  expect_equal(as.double(fit$merr), sqrt(drop(c(1, .5) %*% expected %*% c(1, .5))),
               tolerance = 2e-13)
  expect_equal(as.double(fit$gerr), sqrt(expected[2L, 2L]), tolerance = 2e-13)
  expect_lt(min(eigen(a, symmetric = TRUE, only.values = TRUE)$values), 0)
})
