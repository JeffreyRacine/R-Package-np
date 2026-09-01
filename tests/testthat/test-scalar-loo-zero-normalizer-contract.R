.scalar_loo_regression_bandwidth <- function(response, value) {
  x <- data.frame(x = seq(0, 1, length.out = length(response)))
  npregbw(
    xdat = x, ydat = response, bws = value,
    bandwidth.compute = FALSE, bwtype = "fixed", bwscaling = FALSE,
    ckertype = "epanechnikov", ckerorder = 2L,
    bwmethod = "cv.ls", regtype = "lc"
  )
}

.scalar_loo_regression_core <- function(response, bandwidth, objective,
                                        penalty) {
  getFromNamespace(".npregbw_call_fixed_degree_core", "np")(
    xdat = data.frame(x = seq(0, 1, length.out = length(response))),
    ydat = response, bws = bandwidth,
    invalid.penalty = penalty, penalty.multiplier = 10,
    eval.only = TRUE, objective = objective
  )
}

.scalar_loo_check_core <- function(response, value, penalty) {
  x <- data.frame(x = seq(0, 1, length.out = length(response)))
  scale <- rep(1, length(response))
  bounds <- c(1e-4, 1 - 1e-4)
  template <- nplsqregbw(
    xdat = x, ydat = response, scale = scale, tau = 0.5,
    delta = 0.5, delta.bounds = bounds, bws = value,
    bandwidth.compute = FALSE, bwtype = "fixed", bwscaling = FALSE,
    ckertype = "epanechnikov", ckerorder = 2L,
    bwmethod = "cv.ls", regtype = "lc", invalid.penalty = "dbmax"
  )$reg.bws
  getFromNamespace(".nplsqreg_call_fixed_degree_core", "np")(
    xdat = x, ydat = response, scale = scale, tau = 0.5,
    bws = template, delta = 0.5, delta.bounds = bounds,
    opt.args = list(invalid.penalty = penalty, penalty.multiplier = 10),
    bandwidth.compute = FALSE
  )
}

.scalar_loo_observe <- function(value) {
  response <- sin(seq_len(8L))
  binary <- rep(c(0, 1), 4L)
  cvls.bw <- .scalar_loo_regression_bandwidth(response, value)
  cvks.bw <- .scalar_loo_regression_bandwidth(binary, value)
  run <- function(penalty) list(
    CVLS = .scalar_loo_regression_core(response, cvls.bw, "ls", penalty),
    CVCHECK = .scalar_loo_check_core(response, value, penalty),
    CVKS = .scalar_loo_regression_core(binary, cvks.bw, "ks", penalty)
  )
  list(raw = run("dbmax"), mapped = run("baseline"))
}

test_that("scalar normalized LOO rejects an exact zero normalizer", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  observed <- .scalar_loo_observe(1e-8)
  for (mode in names(observed$raw)) {
    expect_identical(as.double(observed$raw[[mode]]$objective),
                     .Machine$double.xmax, info = mode)
    expect_identical(as.double(observed$raw[[mode]]$invalid.history[[1L]]),
                     1, info = mode)
    expect_true(is.finite(observed$mapped[[mode]]$objective), info = mode)
    expect_lt(observed$mapped[[mode]]$objective, .Machine$double.xmax)
  }
})

test_that("scalar normalized LOO preserves supported rows", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  observed <- .scalar_loo_observe(0.35)
  for (mode in names(observed$raw)) {
    expect_identical(as.double(observed$raw[[mode]]$invalid.history[[1L]]),
                     0, info = mode)
    expect_identical(as.double(observed$mapped[[mode]]$objective),
                     as.double(observed$raw[[mode]]$objective), info = mode)
  }
})
