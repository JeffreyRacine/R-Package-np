.regression_penalty_error <-
  "'penalty.multiplier' must be a finite numeric scalar greater than or equal to 1"

test_that("regression penalty multipliers have a strict public boundary", {
  validate <- getFromNamespace(
    "npValidateRegressionPenaltyMultiplier", "np"
  )

  for (value in c(1, 1 + .Machine$double.eps, 2, 10,
                  .Machine$double.xmax))
    expect_identical(validate(value), as.double(value))

  rejected <- list(
    1 - .Machine$double.eps / 2, 0, -1, NA_real_, NaN, Inf, -Inf,
    numeric(), c(1, 2), TRUE, "1", factor("1"), 1 + 0i
  )
  for (value in rejected)
    expect_error(validate(value), .regression_penalty_error, fixed = TRUE)

  x <- data.frame(x = seq(-1, 1, length.out = 8L))
  y <- c(-2, -1, 0, 0.5, 1, 2, 4, 8)
  expect_error(
    npregbw(
      xdat = x, ydat = y, bws = 1, bandwidth.compute = FALSE,
      bwtype = "fixed", bwscaling = FALSE, ckertype = "uniform",
      bwmethod = "cv.ls", regtype = "lc", penalty.multiplier = 0.5
    ),
    .regression_penalty_error,
    fixed = TRUE
  )
  expect_error(
    npregbw(
      y ~ x, data = data.frame(y = y, x = x$x),
      bws = 1, bandwidth.compute = FALSE,
      bwtype = "fixed", bwscaling = FALSE, ckertype = "uniform",
      bwmethod = "cv.ls", regtype = "lc", invalid.penalty = "dbmax",
      penalty.multiplier = 0.5
    ),
    .regression_penalty_error,
    fixed = TRUE
  )
  expect_error(
    nplsqregbw(
      xdat = x, ydat = y, scale = rep(1, length(y)), tau = 0.5,
      delta = 0.5, delta.bounds = c(0.1, 0.9), bws = 1,
      bandwidth.compute = FALSE, regtype = "lc", ckertype = "uniform",
      invalid.penalty = "baseline", penalty.multiplier = 0.5
    ),
    .regression_penalty_error,
    fixed = TRUE
  )
})

test_that("internal and native regression ingresses reject subunit multipliers", {
  x <- data.frame(x = seq(-1, 1, length.out = 8L))
  y <- c(-2, -1, 0, 0.5, 1, 2, 4, 8)
  bw <- npregbw(
    xdat = x, ydat = y, bws = 1, bandwidth.compute = FALSE,
    bwtype = "fixed", bwscaling = FALSE, ckertype = "uniform",
    bwmethod = "cv.ls", regtype = "lc"
  )
  prepare <- getFromNamespace(".npregbw_nomad_native_prepare_args", "np")
  evaluate <- getFromNamespace(".npregbw_eval_only", "np")

  expect_error(
    evaluate(x, y, bw, penalty.multiplier = 0.5),
    .regression_penalty_error,
    fixed = TRUE
  )
  expect_error(
    prepare(x, y, bw, penalty.multiplier = 0.5),
    .regression_penalty_error,
    fixed = TRUE
  )

  prep <- prepare(x, y, bw)
  expect_error(
    .Call(
      "C_np_regression_bw_eval",
      prep$runo, prep$rord, prep$rcon, prep$y, prep$mysd,
      prep$myopti, prep$myoptd, as.double(bw$bw), 1L,
      prep$penalty_mode, 0.5, prep$degree, prep$bernstein, prep$basis,
      prep$ckerlb, prep$ckerub, PACKAGE = "np"
    ),
    paste0(
      "C_np_regression_bw: penalty.multiplier must be finite ",
      "and greater than or equal to 1"
    ),
    fixed = TRUE
  )

  valid <- .Call(
    "C_np_regression_bw_eval",
    prep$runo, prep$rord, prep$rcon, prep$y, prep$mysd,
    prep$myopti, prep$myoptd, as.double(bw$bw), 1L,
    prep$penalty_mode, prep$penalty_multiplier,
    prep$degree, prep$bernstein, prep$basis, prep$ckerlb, prep$ckerub,
    PACKAGE = "np"
  )
  expect_true(is.finite(valid$fval[[1L]]))
})
