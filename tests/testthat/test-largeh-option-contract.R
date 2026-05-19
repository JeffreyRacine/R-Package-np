test_that("np.largeh toggles continuous large-bandwidth shortcuts", {
  old.options <- options(np.tree = FALSE, np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "np")
  eval_only <- getFromNamespace(".npregbw_eval_only", "np")

  set.seed(42)
  n <- 120L
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- xdat$x1 + xdat$x2 + rnorm(n)

  bws <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "ll",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    bandwidth.compute = FALSE,
    bws = c(1e8, 1e8)
  )

  options(np.largeh = TRUE)
  enabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
  options(np.largeh = FALSE)
  disabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
  options(np.largeh = TRUE)
  reenabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)

  expect_equal(enabled$objective, disabled$objective, tolerance = 1e-10)
  expect_equal(enabled$objective, reenabled$objective, tolerance = 1e-10)
  expect_equal(enabled$num.feval.fast, 1)
  expect_equal(disabled$num.feval.fast, 0)
  expect_equal(reenabled$num.feval.fast, 1)
})

test_that("np.largelambda toggles discrete upper-lambda shortcuts", {
  old.options <- options(np.tree = FALSE, np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "np")
  eval_only <- getFromNamespace(".npregbw_eval_only", "np")

  set.seed(43)
  n <- 120L
  xdat <- data.frame(z = factor(sample(letters[1:3], n, replace = TRUE)))
  ydat <- rnorm(n)

  bws <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "lc",
    bwmethod = "cv.ls",
    ukertype = "aitchisonaitken",
    bandwidth.compute = FALSE,
    bws = 2/3
  )

  options(np.largelambda = TRUE)
  enabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
  options(np.largelambda = FALSE)
  disabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
  options(np.largelambda = TRUE)
  reenabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)

  expect_equal(enabled$objective, disabled$objective, tolerance = 1e-10)
  expect_equal(enabled$objective, reenabled$objective, tolerance = 1e-10)
  expect_equal(enabled$num.feval.fast, 1)
  expect_equal(disabled$num.feval.fast, 0)
  expect_equal(reenabled$num.feval.fast, 1)
})
