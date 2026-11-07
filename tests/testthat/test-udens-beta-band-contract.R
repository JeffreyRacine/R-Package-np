# Fast default-lane contract: admission, signed inversion, the containing hull
# and the dedicated numerical owner. Full coverage/MPI/resource work is external.
test_that("fixed beta band admission preserves the incumbent owners", {
  b <- structure(list(type = "fixed", ckertype = "beta", ckerorder = 2L,
    ckerbound = "fixed", ndim = 1L, ncon = 1L), class = "bandwidth")
  expect_true(.np_udens_beta_request(b, NULL, "bias-corrected", "asymptotic"))
  expect_true(.np_udens_beta_request(b, NULL, "bias-corrected", "bootstrap"))
  expect_false(.np_udens_beta_request(b, list(version = 1L), "bias-corrected", "asymptotic"))
  expect_false(.np_udens_beta_request(b, NULL, "estimate", "asymptotic"))
  expect_false(.np_udens_beta_request(b, NULL, "bias-corrected", "none"))
  for (kernel in c("gaussian", "uniform", "epanechnikov")) {
    b$ckertype <- kernel
    expect_false(.np_udens_beta_request(b, NULL, "bias-corrected", "asymptotic"))
  }
})

test_that("beta residual roots use signed inversion and a single containing hull", {
  roots <- rbind(c(-8, -3, -3, -.2, .5, 1, 2, 2, 4),
                 c(-2, -1, 0, 0, .2, 1, 3, 5, 9))
  rhat <- c(.1, -.2); scale <- c(.2, .4)
  region <- list(lower = c(.3, .8), upper = c(.7, 1.4))
  center <- c(.55, 1.0); alpha <- .8
  k <- floor(ncol(roots)*alpha/4) + 1L
  upper <- ncol(roots) + 1L - k
  ka <- ceiling(ncol(roots)*(1-alpha/2))
  for (band in c("pointwise", "simultaneous")) {
    lower.root <- if (band == "pointwise") apply(roots, 1, function(x) sort(x)[k]) else
      sort(apply(roots, 2, min))[k]
    upper.root <- if (band == "pointwise") apply(roots, 1, function(x) sort(x)[upper]) else
      sort(apply(roots, 2, max))[upper]
    got <- .np_udens_beta_sd_reduce(roots, rhat, scale, region, center, alpha,
      band, "asymmetric-sd")
    expect_equal(got$lower, rhat + region$lower - scale*upper.root, tolerance = 1e-14)
    expect_equal(got$upper, rhat + region$upper - scale*lower.root, tolerance = 1e-14)
    expect_identical(got$construction, "signed-residual-inversion")
    absolute <- if (band == "pointwise") apply(abs(roots), 1, function(x) sort(x)[ka]) else
      sort(apply(abs(roots), 2, max))[ka]
    raw.lower <- rhat + region$lower - scale*absolute
    raw.upper <- rhat + region$upper + scale*absolute
    radius <- pmax(abs(raw.lower-center), abs(raw.upper-center))
    got <- .np_udens_beta_sd_reduce(roots, rhat, scale, region, center, alpha,
      band, "symmetric-sd")
    expect_equal(got$lower, center-radius, tolerance = 1e-14)
    expect_equal(got$upper, center+radius, tolerance = 1e-14)
    expect_true(all(got$lower <= raw.lower & got$upper >= raw.upper))
    expect_identical(got$construction, "residual-absolute-containing-hull")
  }
  roots[1, 1] <- NA_real_
  expect_error(.np_udens_beta_sd_reduce(roots, rhat, scale, region, center, alpha,
    "pointwise", "asymmetric-sd"), "invalid beta plot SD reduction")
})

test_that("plot beta quadrature has the beta probability orientation", {
  for (shape in list(c(.5, .5), c(.9, 1), c(2, 5))) {
    rule <- .np_udens_beta_rule(shape, 16L)
    expect_true(all(rule[, 1] > 0 & rule[, 1] < 1 & rule[, 2] > 0))
    for (k in 0:3)
      expect_equal(sum(rule[, 2]*rule[, 1]^k),
        exp(lbeta(shape[1]+k, shape[2])-lbeta(shape[1], shape[2])), tolerance = 2e-13)
  }
  expect_error(.np_udens_beta_rule(c(0, 1), 16L), "invalid beta plot quadrature geometry", fixed = TRUE)
})

test_that("invalid beta requests fail before RNG and numerical preparation", {
  x <- matrix(seq(.05, .95, length.out = 10L), ncol = 1L)
  expect_error(.np_udens_beta_hidden(x, .1, matrix(c(0, .5), ncol = 1L), .05,
    "asymptotic", "pointwise"), "invalid beta plot hidden provider request")
  expect_error(.np_udens_beta_random_request("multinomial", 10L, 0L, 42L,
    c("Mersenne-Twister", "Inversion", "Rejection")), "invalid beta plot random stream")
})


test_that("beta guide likelihood increments retain Newton descent near stationarity", {
  # Independent 80-digit beta-moment series value. Direct subtraction of
  # near-equal normalizers gives the wrong sign at this ordinary near-uniform fit.
  theta <- c(-0.0019176123663135151, -0.006830269895892126, 0.017261954149090944, 0.011411625904602433, -0.080214918194572979)
  trial <- c(-0.0019174896110850639, -0.0068301485055371142, 0.017262081881624272, 0.011411752008505976, -0.08021493442337238)
  target <- c(-0.99773426984724944, -1.0056384015921216, -0.99049873160507795, -0.99973334755532695, -0.0088476067617440698)
  delta <- trial - theta
  value <- .np_udens_beta_joint_objective_delta(theta, delta, target)
  expect_lt(abs(value - (-1.088271278671157349e-14)), 1e-20)
  expect_lt(value, 0)
  expect_identical(.np_udens_beta_joint_objective_delta(theta, rep(0, 5L), target), 0)
  reverse <- .np_udens_beta_joint_objective_delta(trial, -delta, target)
  expect_lt(abs(reverse + value), 1e-20)
  expect_gt(reverse, 0)
})
