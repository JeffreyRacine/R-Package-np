npscoefbw_default <- getFromNamespace("npscoefbw.default", "np")
npscoefbw_scbandwidth <- getFromNamespace("npscoefbw.scbandwidth", "np")
npscoefbw_start_controls <- getFromNamespace(".npscoefbw_start_controls", "np")
npscoef_default_start_bandwidth <- getFromNamespace(".npscoef_default_start_bandwidth", "np")
npscoef_random_start_bandwidth <- getFromNamespace(".npscoef_random_start_bandwidth", "np")
npscoef_nn_candidate_bandwidth <- getFromNamespace(".npscoef_nn_candidate_bandwidth", "np")
npscoef_candidate_is_admissible <- getFromNamespace(".npscoef_candidate_is_admissible", "np")
np_round_half_to_even <- getFromNamespace(".np_round_half_to_even", "np")
validate_bandwidth <- getFromNamespace("validateBandwidthTF", "np")
npscoef_raw_objective_valid <- getFromNamespace(".npscoefbw_raw_objective_valid", "np")
npscoef_assert_training_radius <- getFromNamespace(".npscoef_nn_assert_training_radius", "np")
semihat_regbw_args <- getFromNamespace(".np_semihat_make_regbw_args", "np")

test_that("npscoefbw surfaces fixed-start controls as formal arguments", {
  expect_true(all(c("scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init", "lbd.init", "hbd.init", "dfac.init") %in%
                    names(formals(npscoefbw_default))))
  expect_true(all(c("scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init", "lbd.init", "hbd.init", "dfac.init") %in%
                    names(formals(npscoefbw_scbandwidth))))
})

test_that("npscoefbw fixed start helpers replay scale-factor defaults", {
  param <- c(1.25, 0.6, 0.45)
  icon <- c(TRUE, FALSE, TRUE)
  iord <- c(FALSE, TRUE, FALSE)
  iuno <- c(FALSE, FALSE, FALSE)
  controls <- npscoefbw_start_controls()

  expect_equal(
    npscoef_default_start_bandwidth(
      param = param,
      bwtype = "fixed",
      nobs = 80L,
      start.controls = controls,
      icon = icon,
      iord = iord,
      iuno = iuno
    ),
    c(0.5, 1.0, 0.5) * param
  )

  set.seed(20260405)
  got <- npscoef_random_start_bandwidth(
    param = param,
    bwtype = "fixed",
    nobs = 80L,
    start.controls = controls,
    icon = icon,
    iord = iord,
    iuno = iuno
  )
  set.seed(20260405)
  expected <- c(
    runif(1, min = 0.1, max = 2.0),
    runif(1, min = 0.5, max = 1.5),
    runif(1, min = 0.1, max = 2.0)
  ) * param
  expect_equal(got, expected)
})

test_that("npscoefbw fixed start helpers split continuous and categorical controls", {
  param <- c(1.25, 0.6, 0.45)
  icon <- c(TRUE, FALSE, TRUE)
  iord <- c(FALSE, TRUE, FALSE)
  iuno <- c(FALSE, FALSE, FALSE)
  controls <- npscoefbw_start_controls(
    scale.factor.init.lower = 0.7,
    scale.factor.init.upper = 0.9,
    scale.factor.init = 1.2,
    lbd.init = 1.1,
    hbd.init = 1.3,
    dfac.init = 0.8
  )

  expect_equal(
    npscoef_default_start_bandwidth(
      param = param,
      bwtype = "fixed",
      nobs = 80L,
      start.controls = controls,
      icon = icon,
      iord = iord,
      iuno = iuno
    ),
    c(1.2, 0.8, 1.2) * param
  )

  set.seed(11)
  got <- npscoef_random_start_bandwidth(
    param = param,
    bwtype = "fixed",
    nobs = 80L,
    start.controls = controls,
    icon = icon,
    iord = iord,
    iuno = iuno
  )
  set.seed(11)
  expected <- c(
    runif(1, min = 0.7, max = 0.9),
    runif(1, min = 1.1, max = 1.3),
    runif(1, min = 0.7, max = 0.9)
  ) * param
  expect_equal(got, expected)
})

test_that("npscoefbw nearest-neighbor start helpers ignore fixed-start controls", {
  param <- c(8, 8, 8)
  controls <- npscoefbw_start_controls(
    scale.factor.init.lower = 0.7,
    scale.factor.init.upper = 0.9,
    scale.factor.init = 1.2,
    lbd.init = 1.1,
    hbd.init = 1.3,
    dfac.init = 0.8
  )

  expect_equal(
    npscoef_default_start_bandwidth(
      param = param,
      bwtype = "generalized_nn",
      nobs = 80L,
      start.controls = controls
    ),
    npscoef_default_start_bandwidth(param = param, bwtype = "generalized_nn", nobs = 80L)
  )

  set.seed(20260405)
  got <- npscoef_random_start_bandwidth(
    param = param,
    bwtype = "generalized_nn",
    nobs = 80L,
    start.controls = controls
  )
  set.seed(20260405)
  expect_equal(got, npscoef_random_start_bandwidth(param = param, bwtype = "generalized_nn", nobs = 80L))
})

test_that("npscoefbw nearest-neighbor helpers preserve categorical coordinates", {
  param <- c(0.25, 8, 0.4)
  icon <- c(FALSE, TRUE, FALSE)
  iord <- c(FALSE, FALSE, TRUE)
  iuno <- c(TRUE, FALSE, FALSE)
  controls <- npscoefbw_start_controls()

  expect_equal(
    npscoef_nn_candidate_bandwidth(
      param = c(0.25, 11.6, 0.4),
      bwtype = "generalized_nn",
      nobs = 80L,
      icon = icon
    ),
    c(0.25, 12, 0.4)
  )
  expect_equal(
    npscoef_default_start_bandwidth(
      param = param,
      bwtype = "generalized_nn",
      nobs = 80L,
      start.controls = controls,
      icon = icon,
      iord = iord,
      iuno = iuno
    ),
    c(0.25, 9, 0.4)
  )

  set.seed(20260902L)
  got <- npscoef_random_start_bandwidth(
    param = param,
    bwtype = "adaptive_nn",
    nobs = 80L,
    start.controls = controls,
    icon = icon,
    iord = iord,
    iuno = iuno
  )
  set.seed(20260902L)
  continuous <- np_round_half_to_even(runif(1L, min = 2, max = 79))
  categorical <- c(0.25, 0.4) *
    runif(2L, min = controls$lbd.init, max = controls$hbd.init)
  expect_equal(got, c(categorical[1L], continuous, categorical[2L]))

  expect_false(npscoef_candidate_is_admissible(
    param = c(0.6, 12, 0.4),
    bwtype = "generalized_nn",
    nobs = 80L,
    upper = c(0.5, Inf, 0.5),
    icon = icon
  ))

  zdat <- data.frame(
    u = factor(rep(c("a", "b"), 40L)),
    x = seq_len(80L),
    o = ordered(rep(c("low", "high"), each = 40L))
  )
  manual <- npscoefbw(
    xdat = data.frame(w = seq_len(80L)),
    ydat = seq_len(80L),
    zdat = zdat,
    bws = c(0.25, 12, 0.4),
    bandwidth.compute = FALSE,
    bwtype = "generalized_nn"
  )
  expect_equal(manual$bw, c(0.25, 12, 0.4))
  helper.args <- semihat_regbw_args(
    source = manual,
    xdat = zdat,
    ydat = seq_len(80L),
    bw = manual$bw
  )
  expect_length(helper.args$ckerlb, ncol(zdat))
  expect_length(helper.args$ckerub, ncol(zdat))
})

test_that("npscoefbw NN radius admission uses exact training multiplicities", {
  zdat <- data.frame(x = c(0, 0, 0, 1, 2), u = factor(c("a", "b", "a", "b", "a")))
  invalid <- list(type = "generalized_nn", icon = c(TRUE, FALSE), bw = c(2, 0.25))
  valid <- invalid
  valid$bw[[1L]] <- 3

  condition <- tryCatch(
    npscoef_assert_training_radius(invalid, zdat, "synthetic npscoef owner"),
    np_nn_candidate_invalid = identity
  )
  expect_s3_class(condition, "np_nn_candidate_invalid")
  expect_identical(condition$owner, "synthetic npscoef owner")
  expect_true(npscoef_assert_training_radius(valid, zdat, "synthetic npscoef owner"))

  penalty <- sqrt(.Machine$double.xmax)
  expect_true(npscoef_raw_objective_valid(1, penalty))
  expect_false(npscoef_raw_objective_valid(penalty, penalty))
  expect_false(npscoef_raw_objective_valid(Inf, penalty))
})

test_that("npscoefbw mixed generalized-NN search crosses a zero-radius plateau", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  data("wage1", package = "np")

  bw <- npscoefbw(
    lwage ~ numdep + smsa |
      married + female + nonwhite + educ + exper + tenure,
    data = wage1,
    bwtype = "generalized_nn",
    regtype = "lc",
    nmulti = 1L,
    optim.maxattempts = 1L,
    optim.maxit = 120L,
    powell.remin = FALSE,
    random.seed = 2026090201L
  )

  expect_true(is.finite(bw$fval) && bw$fval < sqrt(.Machine$double.xmax))
  expect_true(all(bw$bw[bw$icon] >= 2L))
  expect_identical(bw$bw[bw$icon], as.double(as.integer(bw$bw[bw$icon])))
  expect_true(validate_bandwidth(bw))
})
