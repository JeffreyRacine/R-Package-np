npscoefbw_default <- getFromNamespace("npscoefbw.default", "npRmpi")
npscoefbw_scbandwidth <- getFromNamespace("npscoefbw.scbandwidth", "npRmpi")
npscoefbw_start_controls <- getFromNamespace(".npscoefbw_start_controls", "npRmpi")
npscoef_default_start_bandwidth <- getFromNamespace(".npscoef_default_start_bandwidth", "npRmpi")
npscoef_random_start_bandwidth <- getFromNamespace(".npscoef_random_start_bandwidth", "npRmpi")

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
