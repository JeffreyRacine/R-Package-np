test_that("npcdist NOMAD accounting is owner-level on fast-path fits", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)
  suppressPackageStartupMessages(library(np))

  set.seed(41)
  n <- 220L
  x <- runif(n)
  y <- pnorm(2 * sin(2 * pi * x) + rnorm(n))

  fit <- npcdist(
    y ~ x,
    nomad = TRUE,
    search.engine = "nomad",
    degree.max = 2L
  )

  expect_true(any(as.integer(fit$bws$degree.engine) > 0L))
  expect_gt(as.numeric(fit$bws$num.feval.fast[1L]), 0)
  expect_lte(
    as.numeric(fit$bws$num.feval.fast[1L]),
    as.numeric(fit$bws$num.feval[1L])
  )

  ev <- np:::.npcdistbw_eval_only(
    data.frame(x = x),
    data.frame(y = y),
    NULL,
    fit$bws,
    invalid.penalty = "baseline",
    penalty.multiplier = 10
  )

  expect_equal(as.numeric(ev$num.feval), 1)
  expect_equal(
    as.numeric(fit$bws$fval[1L]),
    as.numeric(ev$objective[1L]),
    tolerance = 1e-12
  )
})

test_that("npcdistbw generalized-NN degree endpoints equal their raw replay", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)
  data("wage1", package = "np")

  bw <- np::npcdistbw(
    lwage ~ married + female + nonwhite + educ + exper + tenure,
    data = wage1,
    bwtype = "generalized_nn",
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 1L,
    degree.start = rep.int(0L, 3L),
    degree.max.cycles = 2L,
    degree.verify = FALSE,
    nomad.opts = list(MAX_BB_EVAL = 120L),
    nmulti = 1L,
    powell.remin = FALSE
  )
  raw <- np:::.npcdistbw_eval_only(
    xdat = wage1[c("married", "female", "nonwhite", "educ", "exper", "tenure")],
    ydat = wage1["lwage"],
    gydat = NULL,
    bws = bw,
    invalid.penalty = "dbmax"
  )$objective

  expect_identical(as.double(bw$fval), as.double(raw))
  expect_true(np:::.np_nn_raw_objective_valid(raw))
})
