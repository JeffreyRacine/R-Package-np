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

test_that("conditional uncertainty retains the effective response operator", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(873)
  x <- data.frame(x = rnorm(120))
  y <- data.frame(y = factor(rep(letters[1:3], 40)))
  ex <- data.frame(x = c(-.5, 0, .5)); ey <- y[1:3, , drop = FALSE]
  b <- npcdensbw(xdat = x, ydat = y, bws = c(2/3-1e-6, .4),
                 bandwidth.compute = FALSE)
  z <- npcdens(bws = b, exdat = ex, eydat = ey, gradients = TRUE, se = TRUE)
  # Ask the existing tile for its cached effective operator; an enclosing
  # test session may have intentionally disabled the gate before loading.
  kb <- .npcdhat_make_ybw(b, y)
  Y <- .np_plot_density_factor_weights(
    .np_plot_density_matrix(kb, y, ey, "normal"), kb, y, ey, "normal",
    ordered = FALSE, constant.response = TRUE)
  u <- outer(ex$x, x$x, "-")/.4
  K <- dnorm(u)/.4; D <- -u*dnorm(u)/.4^2
  den <- rowSums(K); mean <- rowSums(K*Y)/den
  constant <- apply(Y, 1L, function(row) all(row == row[1L]))
  variance <- rowSums(K*(Y-mean)^2)/(2*sqrt(pi)*.4*den^2)
  variance[constant] <- 0
  expected <- sqrt(variance)/(.4*sqrt(2))
  expect_equal(as.numeric(gradients(z, se = TRUE)), expected, tolerance = 2e-12)
  joint <- .np_plot_density_conditional(b, x, y, ex, ey, gradient = 1L)
  expect_equal(joint$estimate, gradients(z)[, 1L], tolerance = 2e-12)
  expect_equal(joint$A[constant, , drop = FALSE],
               matrix(0, sum(constant), nrow(x)), tolerance = 0)
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
  expect_equal(joint$estimate, gradients(z)[, 1L], tolerance = 2e-12)
  p <- plot(z, errors = "asymptotic", band = "all", gradients = TRUE,
            neval = 5L, output = "data")
  expect_length(p[[1]]$band.errors[[1]]$simultaneous, 10L)
  expect_true(all(is.finite(p[[1]]$band.errors[[1]]$simultaneous)))
})
