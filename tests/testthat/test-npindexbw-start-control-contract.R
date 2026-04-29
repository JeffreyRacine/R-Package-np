npindexbw_default <- getFromNamespace("npindexbw.default", "npRmpi")
npindexbw_sibandwidth <- getFromNamespace("npindexbw.sibandwidth", "npRmpi")
npindexbw_h_start_controls <- getFromNamespace(".npindexbw_h_start_controls", "npRmpi")
npindex_default_start_bandwidth <- getFromNamespace(".npindex_default_start_bandwidth", "npRmpi")
npindex_random_start_bandwidth <- getFromNamespace(".npindex_random_start_bandwidth", "npRmpi")
npindex_nomad_fixed_start_setup <- getFromNamespace(".npindexbw_nomad_fixed_start_setup", "npRmpi")

test_that("npindexbw surfaces fixed-h start controls as formal arguments", {
  expect_true(all(c("scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init") %in% names(formals(npindexbw_default))))
  expect_true(all(c("scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init") %in% names(formals(npindexbw_sibandwidth))))
})

test_that("npindexbw fixed-h helper replays scale-factor defaults", {
  fit <- c(-1.5, -0.25, 0.1, 0.8, 1.7)
  controls <- npindexbw_h_start_controls()
  h.scale <- getFromNamespace("EssDee", "npRmpi")(fit) * 80^(-1 / 5)

  expect_equal(
    npindex_default_start_bandwidth(fit = fit, bwtype = "fixed", nobs = 80L, start.controls = controls),
    0.5 * h.scale
  )

  set.seed(20260405)
  got <- npindex_random_start_bandwidth(
    fit = fit,
    bwtype = "fixed",
    nobs = 80L,
    start.controls = controls
  )
  set.seed(20260405)
  expect_equal(got, runif(1, min = 0.1, max = 2.0) * h.scale)
})

test_that("npindexbw fixed-h helper honors explicit overrides", {
  fit <- c(-1.5, -0.25, 0.1, 0.8, 1.7)
  controls <- npindexbw_h_start_controls(scale.factor.init.lower = 0.7, scale.factor.init.upper = 0.9, scale.factor.init = 1.2)

  expect_equal(
    npindex_default_start_bandwidth(fit = fit, bwtype = "fixed", nobs = 80L, start.controls = controls),
    1.2 * getFromNamespace("EssDee", "npRmpi")(fit) * 80^(-1 / 5)
  )

  set.seed(11)
  got <- npindex_random_start_bandwidth(
    fit = fit,
    bwtype = "fixed",
    nobs = 80L,
    start.controls = controls
  )
  set.seed(11)
  expect_equal(got, runif(1, min = 0.7, max = 0.9) * getFromNamespace("EssDee", "npRmpi")(fit) * 80^(-1 / 5))
})

test_that("npindexbw nearest-neighbor h starts ignore fixed-h controls", {
  fit <- c(-1.5, -0.25, 0.1, 0.8, 1.7)
  controls <- npindexbw_h_start_controls(scale.factor.init.lower = 0.7, scale.factor.init.upper = 0.9, scale.factor.init = 1.2)

  expect_equal(
    npindex_default_start_bandwidth(fit = fit, bwtype = "generalized_nn", nobs = 80L, start.controls = controls),
    npindex_default_start_bandwidth(fit = fit, bwtype = "generalized_nn", nobs = 80L)
  )

  set.seed(20260405)
  got <- npindex_random_start_bandwidth(
    fit = fit,
    bwtype = "generalized_nn",
    nobs = 80L,
    start.controls = controls
  )
  set.seed(20260405)
  expect_equal(got, npindex_random_start_bandwidth(fit = fit, bwtype = "generalized_nn", nobs = 80L))
})

test_that("npindexbw fixed NOMAD start setup preserves beta starts and replays default h starts", {
  set.seed(42)
  n <- 120L
  x1 <- runif(n, -1, 1)
  x2 <- runif(n, -1, 1)
  y <- x1 - x2 + rnorm(n)
  xmat <- cbind(x1 = x1, x2 = x2)
  degree.search <- list(start.degree = 1L, lower = 0L, upper = 3L, basis = "glp", nobs = n, start.user = FALSE)
  baseline.bws <- list(beta = c(1, 0), bw = 0)

  legacy <- npindex_nomad_fixed_start_setup(
    xmat = xmat,
    ydat = y,
    baseline.bws = baseline.bws,
    degree.search = degree.search,
    nmulti = 4L,
    random.seed = 42L
  )
  replay <- npindex_nomad_fixed_start_setup(
    xmat = xmat,
    ydat = y,
    baseline.bws = baseline.bws,
    degree.search = degree.search,
    nmulti = 4L,
    random.seed = 42L,
    h.start.controls = npindexbw_h_start_controls()
  )

  expect_equal(replay$start_matrix.raw, legacy$start_matrix.raw)
  expect_equal(replay$start_matrix.point, legacy$start_matrix.point)

  custom <- npindex_nomad_fixed_start_setup(
    xmat = xmat,
    ydat = y,
    baseline.bws = baseline.bws,
    degree.search = degree.search,
    nmulti = 4L,
    random.seed = 42L,
    h.start.controls = npindexbw_h_start_controls(scale.factor.init.lower = 0.7, scale.factor.init.upper = 0.9, scale.factor.init = 1.2)
  )

  expect_equal(custom$start_matrix.raw[, legacy$beta.free, drop = FALSE],
               legacy$start_matrix.raw[, legacy$beta.free, drop = FALSE])
  expect_false(isTRUE(all.equal(custom$start_matrix.raw[, legacy$h.col], legacy$start_matrix.raw[, legacy$h.col])))
})
