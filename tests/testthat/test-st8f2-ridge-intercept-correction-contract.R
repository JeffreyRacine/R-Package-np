test_that("R ridge correction preserves signed pristine anchors", {
  correction <- getFromNamespace("npRidgeInterceptCorrection", "npRmpi")
  tiny <- 2^-900
  ridge <- 2^-26

  expect_identical(correction(0, 3, 0), 0)
  expect_equal(correction(ridge, 3, 1), 3 * ridge, tolerance = 0)
  expect_equal(correction(ridge, 3 * tiny, tiny), 3 * ridge, tolerance = 0)
  expect_equal(correction(ridge, -3 * tiny, -tiny), 3 * ridge, tolerance = 0)
  expect_true(is.na(correction(ridge, 0, 0)))
  expect_true(is.na(correction(ridge, 1, Inf)))
  expect_true(is.na(correction(-ridge, 1, 1)))
  expect_equal(correction(1e-300, 1e300, 1e-300), 1e300,
               tolerance = 2e-15)
  expect_error(correction(c(0, ridge), c(1, 2, 3), 1), "not conformable")
})

test_that("scalar zero-ridge correction retains input contracts", {
  correction <- getFromNamespace("npRidgeInterceptCorrection", "npRmpi")
  for (value in c(0, -0, 2^-900, NA_real_, NaN, Inf, -Inf)) {
    expect_identical(correction(0, value, value), 0)
    expect_identical(correction(-0, value, value), 0)
  }
  expect_identical(correction(0, c(a = 2), matrix(3)), 0)
  expect_identical(correction(0L, 2L, 3L), 0)
  expect_identical(correction(0, c(1, 2), c(3, 4)), c(0, 0))
  expect_identical(correction(c(0, 0), 1, 2), c(0, 0))
  expect_error(correction(0, numeric(), 1), "not conformable")
  expect_error(correction(0, 1, numeric()), "not conformable")
  expect_error(correction(0, c(1, 2), c(1, 2, 3)), "not conformable")
  expect_warning(correction(0, "bad", 1), "NAs introduced by coercion")
})

test_that("smooth-coefficient ridge solve reproduces constants at tiny scale", {
  solve.moment <- getFromNamespace(
    ".npscoefbw_nomad_solve_cv_moment_system", "npRmpi"
  )
  maxPenalty <- sqrt(.Machine$double.xmax)
  solve.anchor <- function(anchor) {
    solve.moment(
      tyw = matrix(c(3 * anchor, 0), nrow = 2L),
      tww = array(c(anchor, 0, 0, 0), dim = c(2L, 2L, 1L)),
      W.eval.design = matrix(c(1, 0), nrow = 1L),
      maxPenalty = maxPenalty,
      n.train = 2L
    )[, 1L]
  }

  for (anchor in c(1, 2^-900, -1, -2^-900))
    expect_equal(solve.anchor(anchor), c(3, 0), tolerance = 2e-15)
  expect_identical(solve.anchor(0), rep(maxPenalty, 2L))
})
