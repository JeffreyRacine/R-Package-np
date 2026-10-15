# Primitive quotient and finite case-mass oracles for the LC permutation owner.
test_that("density continuous gradients retain literal explanatory factors", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  set.seed(731)
  x <- data.frame(x = rnorm(120), f = factor(rep(letters[1:3], 40)))
  y <- data.frame(y = sin(x$x) + rnorm(120, sd = .7))
  ex <- x[c(11, 31, 51), , drop = FALSE]
  ey <- y[c(11, 31, 51), , drop = FALSE]
  for (cdf in c(FALSE, TRUE)) {
    b <- if (cdf) npcdistbw(xdat = x, ydat = y, bws = c(.5, .6, .663),
      bandwidth.compute = FALSE) else npcdensbw(xdat = x, ydat = y,
      bws = c(.5, .6, .663), bandwidth.compute = FALSE)
    actual <- if (cdf) npcdist(bws = b, txdat = x, tydat = y,
      exdat = ex, eydat = ey, gradients = TRUE) else npcdens(bws = b,
      txdat = x, tydat = y, exdat = ex, eydat = ey, gradients = TRUE)
    joint <- .np_plot_density_conditional(b, x, y, ex, ey, cdf,
      gradient = 1L)
    expect_equal(joint$estimate, gradients(actual)[, 1L], tolerance = 1e-12)
    # Independent literal unordered-kernel weights: no effective-factor
    # adjustment is part of the conditional estimator's requested target.
    raw <- outer(as.integer(ex$f), as.integer(x$f),
      function(a, b) ifelse(a == b, 1-.663, .663/2))
    u <- outer(ex$x, x$x, "-")/.6
    K <- dnorm(u)/.6 * raw
    D <- -u*dnorm(u)/.6^2 * raw
    v <- outer(ey$y, y$y, "-")/.5
    Y <- if (cdf) pnorm(sqrt(2)*0.7071067810*v) else dnorm(v)/.5
    target <- function(mass) {
      s <- drop(K %*% mass); ds <- drop(D %*% mass)
      (drop((D*Y) %*% mass)*s - drop((K*Y) %*% mass)*ds)/s^2
    }
    h <- 1e-4
    oracle <- vapply(seq_len(120), function(i) {
      plus <- minus <- rep(1, 120)
      plus[i] <- 1+h; minus[i] <- 1-h
      (target(plus)-target(minus))/(2*h)
    }, numeric(3))
    oracle <- oracle-rowMeans(oracle)
    expect_equal(joint$A, oracle, tolerance = 2e-10)
    corr <- function(z) {
      V <- tcrossprod(z)
      V/sqrt(outer(diag(V), diag(V)))
    }
    expect_equal(corr(joint$A), corr(oracle), tolerance = 2e-7)
    expect_error(.np_plot_density_joint(b, x, ex, rep(.01, 3), y, ey,
      cdf, gradient = 1L, estimate = joint$estimate+.1), "does not reproduce")
  }
})

