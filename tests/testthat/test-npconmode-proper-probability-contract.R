test_that("npconmode proper helper enforces binary complement probabilities", {
  helper <- getFromNamespace(".npConmodeProperProbabilities", "npRmpi")
  raw <- matrix(c(-0.2, 1.2,
                  0.4, 0.8,
                  0.3, 0.3), ncol = 2L, byrow = TRUE)

  out <- helper(raw, levels = c("0", "1"), proper = TRUE)

  expect_true(isTRUE(out$proper.requested))
  expect_true(isTRUE(out$proper.applied))
  expect_equal(rowSums(out$probabilities), rep(1, nrow(raw)), tolerance = 1e-12)
  expect_true(all(out$probabilities >= -1e-12))
  expect_equal(out$probabilities[, 2L], 1 - out$probabilities[, 1L], tolerance = 1e-12)
  expect_identical(out$proper.info$reason, "projected")
})

test_that("npconmode proper defaults follow the canonical regression type", {
  effective <- getFromNamespace(".npConmodeEffectiveProper", "npRmpi")

  expect_false(effective(list(regtype = "lc"), NULL))
  expect_true(effective(list(regtype = "ll"), NULL))
  expect_true(effective(list(regtype = "lp"), NULL))
  expect_false(effective(list(regtype = "lp", regtype.engine = "lc"), NULL))
  expect_true(effective(list(regtype = "lc", regtype.engine = "lp"), NULL))
  expect_false(effective(list(regtype = "lp"), FALSE))
  expect_true(effective(list(regtype = "lc"), TRUE))
})

