# Exact affine/constant-response identities, independent of the fitted solver.
test_that("conditional LP CDF reproduces a saturated response", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  data("Engel95", package = "npRmpi", envir = environment())
  b <- npcdistbw(food ~ logexp, data = Engel95, regtype = "lp",
                degree = 3L, bwmethod = "normal-reference")
  ex <- data.frame(logexp = c(min(Engel95$logexp), median(Engel95$logexp),
                             max(Engel95$logexp)))
  ey <- data.frame(food = rep(max(Engel95$food) + 100, 3L))
  z <- npcdist(bws = b, exdat = ex, eydat = ey, gradients = TRUE, se = TRUE)
  expect_equal(as.numeric(fitted(z)), rep(1, 3L), tolerance = 0)
  expect_equal(as.numeric(gradients(z)), rep(0, 3L), tolerance = 0)
  expect_equal(as.numeric(se(z)), rep(0, 3L), tolerance = 0)
})

test_that("conditional uncertainty retains the literal response operator", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(873)
  x <- data.frame(x = rnorm(120))
  y <- data.frame(y = factor(rep(letters[1:3], 40)))
  ex <- data.frame(x = c(-.5, 0, .5)); ey <- y[1:3, , drop = FALSE]
  # A positive-variance oracle independent of the plotting/kernel helpers.
  # lambda=.25 is outside both near-bound gates for every accepted tolerance
  # in (0,.5), including an enclosing session's cached option values.
  literal <- npcdensbw(xdat = x, ydat = y, bws = c(.25, .4),
                       bandwidth.compute = FALSE)
  literal.fit <- npcdens(bws = literal, exdat = ex, eydat = ey,
                        gradients = TRUE, se = TRUE)
  Y.literal <- ifelse(outer(as.character(ey$y), as.character(y$y), "=="),
                      .75, .25/2)
  u.literal <- outer(ex$x, x$x, "-")/.4
  K.literal <- dnorm(u.literal)/.4
  den.literal <- rowSums(K.literal)
  mean.literal <- rowSums(K.literal*Y.literal)/den.literal
  variance.literal <- rowSums(K.literal*(Y.literal-mean.literal)^2)/
    (2*sqrt(pi)*.4*den.literal^2)
  expected.literal <- sqrt(variance.literal)/(.4*sqrt(2))
  expect_true(all(expected.literal > .05))
  expect_equal(as.numeric(gradients(literal.fit, se = TRUE)), expected.literal,
               tolerance = 2e-12)

  b <- npcdensbw(xdat = x, ydat = y, bws = c(2/3-1e-6, .4),
                 bandwidth.compute = FALSE)
  z <- npcdens(bws = b, exdat = ex, eydat = ey, gradients = TRUE, se = TRUE)
  # Literal nominal-response weights retain small differences near the bound.
  # Subtract a row constant before quotient differentiation: its derivative and
  # case-mass influence are exactly zero, avoiding cancellation in the oracle.
  lambda <- 2/3-1e-6
  Y <- ifelse(outer(as.character(ey$y), as.character(y$y), "=="),
              1-lambda, lambda/2)
  Z <- Y-Y[, 1L]
  u <- outer(ex$x, x$x, "-")/.4
  K <- dnorm(u)/.4; D <- -u*dnorm(u)/.4^2
  den <- rowSums(K); mean <- rowSums(K*Z)/den
  variance <- rowSums(K*(Z-mean)^2)/(2*sqrt(pi)*.4*den^2)
  expected <- sqrt(variance)/(.4*sqrt(2))
  expect_true(all(expected > 0))
  expect_equal(as.numeric(gradients(z, se = TRUE)), expected, tolerance = 2e-12)
  H <- D/den-(K/den)*(rowSums(D)/den)
  expected.gradient <- rowSums(H*Z)
  expected.A <- H*(Z-mean)-(K/den)*expected.gradient
  expected.A <- expected.A-rowMeans(expected.A)
  expect_true(all(rowSums(expected.A^2) > 0))
  joint <- .np_plot_density_conditional(b, x, y, ex, ey, gradient = 1L)
  expect_equal(joint$estimate, expected.gradient, tolerance = 2e-12)
  expect_equal(joint$A, expected.A, tolerance = 2e-12)
  # The public quotient subtracts uncentered O(1) terms. Its absolute error
  # scales with their magnitudes, not with the near-zero derivative. gamma
  # bounds summation roundoff; the factor eight covers the quotient chain.
  eps <- .Machine$double.eps
  gamma <- (nrow(x)+8)*eps/(1-(nrow(x)+8)*eps)
  scale <- rowSums(abs(D*Y))/den +
    abs(rowSums(K*Y)/den)*rowSums(abs(D))/den
  bound <- 8*gamma*scale
  expect_true(all(bound < abs(expected.gradient)*1e-4))
  expect_true(all(abs(gradients(z)[, 1L]-expected.gradient) <= bound))

})

# Nominal Y coordinates use the normal kernel even in a CDF product.
test_that("conditional CDF response factors retain their native operator", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(982)
  x <- data.frame(x = rnorm(80))
  y <- data.frame(y = factor(rep(letters[1:2], 40)))
  b <- npcdistbw(xdat = x, ydat = y, bws = c(.5-1e-6, .4),
                 bandwidth.compute = FALSE)
  ex <- x[1:3, , drop = FALSE]; ey <- y[1:3, , drop = FALSE]
  z <- npcdist(bws = b, exdat = ex, eydat = ey, gradients = TRUE, se = TRUE)
  joint <- .np_plot_density_conditional(b, x, y, ex, ey, cdf = TRUE, gradient = 1L)
  # A nominal response still uses its literal unordered kernel in a CDF.
  lambda <- .5-1e-6
  Y <- ifelse(outer(as.character(ey$y), as.character(y$y), "=="),
              1-lambda, lambda)
  Z <- Y-Y[, 1L]
  u <- outer(ex$x, x$x, "-")/.4
  K <- dnorm(u)/.4; D <- -u*dnorm(u)/.4^2
  den <- rowSums(K)
  H <- D/den-(K/den)*(rowSums(D)/den)
  expected.gradient <- rowSums(H*Z)
  expect_equal(joint$estimate, expected.gradient, tolerance = 2e-12)
  # Same uncentered quotient forward-error allowance as the density fixture.
  eps <- .Machine$double.eps
  gamma <- (nrow(x)+8)*eps/(1-(nrow(x)+8)*eps)
  scale <- rowSums(abs(D*Y))/den +
    abs(rowSums(K*Y)/den)*rowSums(abs(D))/den
  bound <- 8*gamma*scale
  expect_true(all(bound < abs(expected.gradient)*1e-4))
  expect_true(all(abs(gradients(z)[, 1L]-expected.gradient) <= bound))
  p <- plot(z, errors = "asymptotic", band = "all", gradients = TRUE,
            neval = 5L, output = "data")
  expect_length(p[[1]]$band.errors[[1]]$simultaneous, 10L)
  expect_true(all(is.finite(p[[1]]$band.errors[[1]]$simultaneous)))
})
