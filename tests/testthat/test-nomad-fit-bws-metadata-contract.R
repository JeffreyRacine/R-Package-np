test_that("nomad fit metadata restore keeps authoritative bandwidth telemetry", {
  authoritative <- list(
    method = "cv.ls",
    pmethod = "Least Squares Cross-Validation",
    fval = 0.25,
    ifval = 0.5,
    num.feval = 259,
    num.feval.fast = 43,
    fval.history = c(0.6, 0.25),
    eval.history = c(216, 43),
    invalid.history = c(0, 0),
    timing = 1.75,
    timing.profile = list(total = 1.75),
    degree.search = list(mode = "nomad+powell"),
    nomad.shortcut = list(enabled = TRUE),
    nomad.time = 1.2,
    powell.time = 0.3,
    total.time = 1.5
  )
  degraded <- list(
    method = "cv.ls",
    pmethod = "Manual",
    fval = 0.25,
    ifval = 0.25,
    num.feval = 216,
    num.feval.fast = 0,
    nomad.time = 1.2,
    powell.time = 0.3,
    total.time = 1.5
  )
  result <- list(
    fit.time = 2.0,
    bws = degraded
  )

  restored <- npRmpi:::.npRmpi_restore_nomad_fit_bws_metadata(result, authoritative)

  expect_identical(restored$bws$method, authoritative$method)
  expect_identical(restored$bws$pmethod, authoritative$pmethod)
  expect_equal(restored$bws$fval, authoritative$fval)
  expect_equal(restored$bws$ifval, authoritative$ifval)
  expect_equal(restored$bws$num.feval, authoritative$num.feval)
  expect_equal(restored$bws$num.feval.fast, authoritative$num.feval.fast)
  expect_identical(restored$bws$fval.history, authoritative$fval.history)
  expect_identical(restored$bws$eval.history, authoritative$eval.history)
  expect_identical(restored$bws$invalid.history, authoritative$invalid.history)
  expect_identical(restored$bws$timing, authoritative$timing)
  expect_identical(restored$bws$timing.profile, authoritative$timing.profile)
  expect_identical(restored$bws$degree.search, authoritative$degree.search)
  expect_identical(restored$bws$nomad.shortcut, authoritative$nomad.shortcut)
  expect_equal(restored$nomad.time, authoritative$nomad.time)
  expect_equal(restored$powell.time, authoritative$powell.time)
  expect_equal(restored$optim.time, authoritative$total.time)
  expect_equal(restored$total.time, authoritative$total.time + result$fit.time)
})
