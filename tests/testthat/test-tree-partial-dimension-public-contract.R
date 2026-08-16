test_that("exact partial-dimension tree plan preserves public descendants", {
  old_options <- options(np.messages = FALSE, np.largeh = FALSE,
                         np.macMseries.accelerate = FALSE)
  on.exit(options(old_options), add = TRUE)

  set.seed(20260819L)
  n <- 127L
  x <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
  y <- sin(2 * pi * x$x3) + 0.2 * rnorm(n)
  ranges <- vapply(x, function(value) diff(range(value)), numeric(1L))
  bandwidths <- 1.25 * ranges / sqrt(5)
  bandwidths[[3L]] <- 0.16 * ranges[[3L]] / sqrt(5)

  reg_bw <- npregbw(
    xdat = x, ydat = y, bws = bandwidths,
    bandwidth.compute = FALSE, bwscaling = FALSE,
    bwtype = "fixed", regtype = "lc", ckertype = "epanechnikov"
  )
  eval_only <- getFromNamespace(".npregbw_eval_only", "np")
  options(np.tree = FALSE)
  reg_objective_dense <- eval_only(x, y, reg_bw)
  reg_fit_dense <- npreg(bws = reg_bw, se = TRUE, gradients = FALSE)
  options(np.tree = TRUE)
  reg_objective_tree <- eval_only(x, y, reg_bw)
  reg_fit_tree <- npreg(bws = reg_bw, se = TRUE, gradients = FALSE)

  expect_equal(reg_objective_tree$objective,
               reg_objective_dense$objective, tolerance = 2e-11)
  for (field in c("mean", "merr", "resid"))
    expect_equal(reg_fit_tree[[field]], reg_fit_dense[[field]],
                 tolerance = 2e-11, info = paste("npreg", field))

  density_bw <- npudensbw(
    dat = x, bws = bandwidths,
    bandwidth.compute = FALSE, bwscaling = FALSE,
    bwtype = "fixed", ckertype = "epanechnikov"
  )
  options(np.tree = FALSE)
  density_dense <- npudens(bws = density_bw)
  options(np.tree = TRUE)
  density_tree <- npudens(bws = density_bw)
  for (field in c("dens", "derr", "log_likelihood"))
    expect_equal(density_tree[[field]], density_dense[[field]],
                 tolerance = 2e-11, info = paste("npudens", field))

  distribution_bw <- npudistbw(
    dat = x, bws = bandwidths,
    bandwidth.compute = FALSE, bwscaling = FALSE,
    bwtype = "fixed", ckertype = "epanechnikov"
  )
  options(np.tree = FALSE)
  distribution_dense <- npudist(bws = distribution_bw)
  options(np.tree = TRUE)
  distribution_tree <- npudist(bws = distribution_bw)
  for (field in c("dist", "derr"))
    expect_equal(distribution_tree[[field]], distribution_dense[[field]],
                 tolerance = 2e-11, info = paste("npudist", field))
})
