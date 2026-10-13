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

test_that("density plot factors follow the incumbent level and contrast owners", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.largelambda = TRUE)
  on.exit(options(old))
  set.seed(1060)
  x <- data.frame(x = rnorm(100), f = factor(rep(letters[1:3], length.out = 100)))
  ex <- x[c(15, 30, 45), , drop = FALSE]
  b <- npudensbw(dat = x, bws = c(.45, .66), bandwidth.compute = FALSE)
  actual <- npudens(bws = b, tdat = x, edat = ex)
  joint <- .np_plot_density_unconditional(b, x, ex)
  expect_equal(joint$estimate, fitted(actual), tolerance = 1e-12)
  y <- data.frame(y = sin(x$x))
  ey <- y[c(15, 30, 45), , drop = FALSE]
  for (degree in c(0L, 2L)) {
    b <- npcdistbw(xdat = x, ydat = y, bws = c(.4, .45, .66),
      bandwidth.compute = FALSE, regtype = "lp", degree = degree)
    actual <- npcdist(bws = b, txdat = x, tydat = y, exdat = ex, eydat = ey,
      gradients = TRUE)
    joint <- .np_plot_density_conditional(b, x, y, ex, ey, cdf = TRUE, gradient = 2L)
    expect_equal(joint$estimate, gradients(actual)[, 2L], tolerance = 1e-10)
  }
})

test_that("centered density LP influence preserves constant response equivariance", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  W <- outer(seq(.1, .9, length.out = 60), 0:2, "^")
  w <- rep(0, 60)
  w[c(12, 45)] <- 1
  identity <- .np_inid_lp_design_identity(W)
  response <- cbind(seq_len(60)/60, rep(1, 60))
  a <- .np_plot_density_lp_score(W, w, response, c(1, .5, .25), 60, identity)
  shifted <- .np_plot_density_lp_score(W, w, response + 1, c(1, .5, .25), 60, identity)
  expect_true(all(a$ridge > 0))
  expect_equal(shifted$estimate, a$estimate + 1, tolerance = 1e-12)
  expect_equal(shifted$score, a$score, tolerance = 1e-12)
  expect_equal(a$score[2L, ], rep(0, 60), tolerance = 0)
})

test_that("density reconstruction allows rounding without inflating zero SEs", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 80))
  ex <- data.frame(x = c(-20, 20))
  bw <- npudistbw(dat = x, bws = .5, bandwidth.compute = FALSE)
  target <- .np_plot_density_unconditional(bw, x, ex, cdf = TRUE)$estimate
  actual <- .np_plot_density_joint(bw, x, ex, se = c(0, 0), cdf = TRUE,
    estimate = target + 1e-13)
  expect_true(is.finite(actual$critical))
  expect_error(.np_plot_density_joint(bw, x, ex, se = c(0, 0), cdf = TRUE,
    estimate = target + 1e-7), "does not reproduce")
  scaled <- .np_plot_density_rescale(matrix(0, 2, 80), c(0, 0))
  expect_identical(scaled$se, c(0, 0))
  expect_identical(scaled$factor, matrix(0, 2, 80))
})


test_that("ordered distribution influence follows the native integral owner", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.largelambda = TRUE)
  on.exit(options(old))
  set.seed(516)
  x <- data.frame(x = rnorm(140), o = ordered(rep(1:3, length.out = 140)))
  ex <- x[1:7, , drop = FALSE]
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    h <- if (type == "fixed") .4 else 40
    bw <- npudistbw(dat = x, bws = c(h, .995), bwtype = type,
      okertype = "liracine", bandwidth.compute = FALSE)
    actual <- npudist(bws = bw, tdat = x, edat = ex)
    joint <- .np_plot_density_unconditional(bw, x, ex, cdf = TRUE)
    expect_equal(joint$estimate, fitted(actual), tolerance = 1e-12)
  }
})
