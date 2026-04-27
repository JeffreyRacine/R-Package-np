library(npRmpi)

test_that("old continuous scale-factor search names fail with rename diagnostics", {
  dat <- data.frame(x = rnorm(12), y = rnorm(12))

  expect_error(
    npregbw(y ~ x, data = dat, cfac.init = 0.5),
    "'cfac.init' has been renamed to 'scale.factor.init'",
    fixed = TRUE
  )
  expect_error(
    npudensbw(~ x, data = dat, lbc.init = 0.1),
    "'lbc.init' has been renamed to 'scale.factor.init.lower'",
    fixed = TRUE
  )
  expect_error(
    npcdensbw(y ~ x, data = dat, hbc.init = 2),
    "'hbc.init' has been renamed to 'scale.factor.init.upper'",
    fixed = TRUE
  )
  expect_error(
    npindexbw(y ~ x + I(x^2), data = dat, scale.factor.lower.bound = 0.2),
    "'scale.factor.lower.bound' has been renamed to 'scale.factor.search.lower'",
    fixed = TRUE
  )
})

test_that("scale-factor search lower metadata uses canonical field with legacy read fallback", {
  canonical <- list(scale.factor.search.lower = 0.25)
  legacy <- list(scale.factor.lower.bound = 0.25)
  conflict <- list(scale.factor.search.lower = 0.25, scale.factor.lower.bound = 0.5)

  expect_equal(npRmpi:::npGetScaleFactorSearchLower(canonical), 0.25, tolerance = 0)
  expect_equal(npRmpi:::npGetScaleFactorSearchLower(legacy), 0.25, tolerance = 0)
  expect_error(
    npRmpi:::npGetScaleFactorSearchLower(conflict),
    "conflicts with legacy 'scale.factor.lower.bound'",
    fixed = TRUE
  )

  object <- npRmpi:::npSetScaleFactorSearchLower(list(scale.factor.lower.bound = 0.1), 0.3)
  expect_equal(object$scale.factor.search.lower, 0.3, tolerance = 0)
  expect_null(object$scale.factor.lower.bound)
})

test_that("protected non-continuous initialization names remain separate", {
  npscoefbw_default <- getFromNamespace("npscoefbw.default", "npRmpi")
  npudensbw_default <- getFromNamespace("npudensbw.default", "npRmpi")

  expect_true(all(c("dfac.init", "lbd.init", "hbd.init") %in%
                    names(formals(npscoefbw_default))))
  expect_true(all(c("dfac.init", "lbd.init", "hbd.init",
                    "initc.dir", "initd.dir",
                    "scale.init.categorical.sample") %in%
                    names(formals(npudensbw_default))))
})
