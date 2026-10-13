# Small public-contract sentinels. Broad LP/NN/kernel/renderer matrices and timing
# evidence live in the external integration campaign, not the default test lane.
test_that("density joint calibration retains independently computed dependence", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  x <- data.frame(x = seq(-3, 3, length.out = 120L))
  ex <- data.frame(x = seq(-.4, .4, length.out = 7L))
  bw <- npudensbw(dat = x, bws = .65, bandwidth.compute = FALSE)
  # Primitive Gaussian contributions provide a reference independent of the
  # package kernel, influence, and covariance constructors.
  K <- outer(ex$x, x$x, function(a, b) dnorm((a-b)/.65)/.65)/nrow(x)
  A <- K - rowMeans(K)
  ref <- tcrossprod(A)
  ref <- ref / sqrt(outer(diag(ref), diag(ref)))
  owner <- .np_plot_density_unconditional(bw, x, ex)
  expect_equal(owner$estimate, rowSums(K), tolerance = 1e-12)
  expect_equal(tcrossprod(owner$A), tcrossprod(A), tolerance = 1e-12)
  # Independent Cholesky simulation checks the resulting critical value. The
  # tolerance allows Monte Carlo variation, but excludes pointwise, Bonferroni
  # and independent-coordinate substitutions for this correlated grid.
  reference <- .np_with_seed(904L, {
    z <- t(chol(ref)) %*% matrix(rnorm(7L * 99999L), 7L)
    unname(quantile(apply(abs(z), 2L, max), .95, type = 1L))
  })
  actual <- .np_with_seed(903L, .np_plot_density_joint(bw, x, ex,
    se = rep(1, 7L), estimate = rowSums(K)))
  expect_lte(abs(actual$critical - reference), .06)
  expect_gt(actual$critical, qnorm(.975) + .1)
  expect_lt(actual$critical, qnorm(1-.05/14) - .1)
})

test_that("density Gaussian plots preserve marginal bands, seed and grid size", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  set.seed(834)
  x <- data.frame(x = rnorm(80))
  y <- data.frame(y = .3 * x$x + rnorm(80))
  bw <- npcdensbw(xdat = x, ydat = y, bws = c(.7, .7),
                  bandwidth.compute = FALSE, regtype = "lp", degree = 2L)
  fit <- npcdens(bws = bw, txdat = x, tydat = y,
                 exdat = x[1:3, , drop = FALSE], eydat = y[1:3, , drop = FALSE])
  rng <- .Random.seed
  out <- plot(fit, errors = "asymptotic", band = "all", output = "data", neval = 5)
  expect_identical(.Random.seed, rng)
  obj <- out[[1L]]
  expect_equal(obj$band.calibration$draws, 9999L)
  expect_equal(obj$band.calibration$grid.size, 25L)
  b <- obj$band.errors
  expect_true(all(b$pointwise <= b$simultaneous))
  expect_true(all(b$simultaneous <= b$bonferroni))
  # This catches the former conditional-surface nrow(vector) multiplicity bug.
  expect_equal(b$bonferroni / b$pointwise,
               matrix(qnorm(1 - .05 / 50) / qnorm(.975), 25, 2))
  single <- plot(fit, errors = "asymptotic", band = "bonferroni",
                 output = "data", neval = 5)
  expect_equal(single[[1L]]$conderr, b$bonferroni)
  gradient <- plot(fit, errors = "asymptotic", band = "all", output = "data",
                   neval = 5, perspective = FALSE, gradients = TRUE)
  expect_length(gradient[[1L]]$band.calibration, 1L)
  expect_true(all(is.finite(gradient[[1L]]$band.errors[[1L]]$simultaneous)))
  expect_error(plot(fit, errors = "asymptotic", band = "pointwise",
                    gaussian.draws = 9999), "gaussian.draws requires")
})

test_that("unconditional NN Gaussian plots and unavailable grids keep their shape", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  set.seed(835)
  x <- data.frame(x = rnorm(80))
  bw <- npudistbw(dat = x, bws = 25, bwtype = "adaptive_nn",
                  bandwidth.compute = FALSE)
  fit <- npudist(bws = bw, tdat = x)
  out <- plot(fit, errors = "asymptotic", band = "all", output = "data", neval = 7)
  expect_equal(out[[1L]]$band.calibration$grid.size, 7L)
  expect_true(all(is.finite(out[[1L]]$band.errors$simultaneous)))
  bw <- npudensbw(dat = x, bws = .7, ckerorder = 4,
                  bandwidth.compute = FALSE)
  fit <- npudens(bws = bw, tdat = x)
  expect_warning(out <- plot(fit, errors = "asymptotic", band = "simultaneous",
                              output = "data", neval = 7), "kernel specification")
  expect_equal(nrow(out[[1L]]$derr), 7L)
  expect_true(all(is.na(out[[1L]]$derr)))
})
