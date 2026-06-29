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

test_that("np.largeh toggles one-step continuous bandwidth-search fast counts", {
  old.options <- options(np.messages = FALSE, np.tree = FALSE,
                         np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "np")

  set.seed(101)
  n <- 80L
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
  enabled <- npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    regtype = "ll",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    bandwidth.compute = TRUE,
    nmulti = 1L,
    itmax = 1L
  )

  options(np.largeh = FALSE)
  disabled <- npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    regtype = "ll",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    bandwidth.compute = TRUE,
    nmulti = 1L,
    itmax = 1L
  )

  expect_equal(enabled$fval, disabled$fval, tolerance = 1e-10)
  expect_gt(as.numeric(enabled$num.feval.fast[1L]), 0)
  expect_gt(as.numeric(enabled$num.feval.fast[1L]),
            as.numeric(disabled$num.feval.fast[1L]))
})

test_that("np.largeh and np.largelambda both gate mixed fast objective rows", {
  old.options <- options(np.messages = FALSE, np.tree = FALSE,
                         np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "np")
  npreg_eval_only <- getFromNamespace(".npregbw_eval_only", "np")
  npscoefbw <- getFromNamespace("npscoefbw", "np")
  npscoef_eval_only <- getFromNamespace(".npscoefbw_eval_only", "np")

  set.seed(107)
  n <- 80L
  xdat <- data.frame(
    x = runif(n),
    z = factor(sample(letters[1:3], n, replace = TRUE))
  )
  ydat <- rnorm(n)

  rbw <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "lc",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    ukertype = "aitchisonaitken",
    bandwidth.compute = FALSE,
    bws = c(1e8, 2 / 3)
  )

  options(np.largeh = TRUE, np.largelambda = TRUE)
  reg.both <- npreg_eval_only(xdat = xdat, ydat = ydat, bws = rbw)
  options(np.largeh = FALSE, np.largelambda = TRUE)
  reg.noh <- npreg_eval_only(xdat = xdat, ydat = ydat, bws = rbw)
  options(np.largeh = TRUE, np.largelambda = FALSE)
  reg.nolam <- npreg_eval_only(xdat = xdat, ydat = ydat, bws = rbw)
  options(np.largeh = TRUE, np.largelambda = TRUE)
  reg.both2 <- npreg_eval_only(xdat = xdat, ydat = ydat, bws = rbw)

  expect_equal(reg.both$objective, reg.noh$objective, tolerance = 1e-10)
  expect_equal(reg.both$objective, reg.nolam$objective, tolerance = 1e-10)
  expect_equal(reg.both$objective, reg.both2$objective, tolerance = 1e-10)
  expect_equal(as.numeric(reg.both$num.feval.fast[1L]), 1)
  expect_equal(as.numeric(reg.noh$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(reg.nolam$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(reg.both2$num.feval.fast[1L]), 1)

  set.seed(108)
  sxdat <- data.frame(x = runif(n))
  szdat <- data.frame(
    z = runif(n),
    g = factor(sample(letters[1:3], n, replace = TRUE))
  )
  sydat <- rnorm(n)

  sbw <- npscoefbw(
    xdat = sxdat,
    zdat = szdat,
    ydat = sydat,
    bwmethod = "cv.ls",
    regtype = "lc",
    bandwidth.compute = FALSE,
    bws = c(1e8, 2 / 3)
  )

  options(np.largeh = TRUE, np.largelambda = TRUE)
  sc.both <- npscoef_eval_only(xdat = sxdat, zdat = szdat, ydat = sydat, bws = sbw)
  options(np.largeh = FALSE, np.largelambda = TRUE)
  sc.noh <- npscoef_eval_only(xdat = sxdat, zdat = szdat, ydat = sydat, bws = sbw)
  options(np.largeh = TRUE, np.largelambda = FALSE)
  sc.nolam <- npscoef_eval_only(xdat = sxdat, zdat = szdat, ydat = sydat, bws = sbw)
  options(np.largeh = TRUE, np.largelambda = TRUE)
  sc.both2 <- npscoef_eval_only(xdat = sxdat, zdat = szdat, ydat = sydat, bws = sbw)

  expect_equal(sc.both$objective, sc.noh$objective, tolerance = 1e-10)
  expect_equal(sc.both$objective, sc.nolam$objective, tolerance = 1e-10)
  expect_equal(sc.both$objective, sc.both2$objective, tolerance = 1e-10)
  expect_equal(as.numeric(sc.both$num.feval.fast[1L]), 1)
  expect_equal(as.numeric(sc.noh$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(sc.nolam$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(sc.both2$num.feval.fast[1L]), 1)
})

test_that("np.largeh toggles conditional density eval-only fast counts", {
  old.options <- options(np.messages = FALSE, np.tree = FALSE,
                         np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npcdensbw <- getFromNamespace("npcdensbw", "np")
  eval_only <- getFromNamespace(".npcdensbw_eval_only", "np")

  set.seed(105)
  n <- 80L
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(y = runif(n))

  bws <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "lc",
    bwmethod = "cv.ml",
    cxkertype = "epanechnikov",
    cykertype = "epanechnikov",
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
  expect_equal(as.numeric(enabled$num.feval.fast[1L]), 1)
  expect_equal(as.numeric(disabled$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(reenabled$num.feval.fast[1L]), 1)
})

test_that("np.largeh toggles unconditional distribution bandwidth fast counts", {
  old.options <- options(np.messages = FALSE, np.tree = FALSE,
                         np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npudistbw <- getFromNamespace("npudistbw", "np")

  set.seed(104)
  n <- 80L
  dat <- data.frame(y = runif(n))

  bws <- npudistbw(
    dat = dat,
    bwmethod = "cv.cdf",
    ckertype = "epanechnikov",
    bandwidth.compute = FALSE,
    ngrid = 20L
  )
  bws$bw[] <- 1e8
  bws$bandwidth[] <- 1e8
  bws$sfactor[] <- 1e8

  options(np.largeh = TRUE)
  enabled <- npudistbw(
    dat = dat,
    bws = bws,
    bandwidth.compute = TRUE,
    nmulti = 1L,
    itmax = 1L,
    ngrid = 20L
  )
  options(np.largeh = FALSE)
  disabled <- npudistbw(
    dat = dat,
    bws = bws,
    bandwidth.compute = TRUE,
    nmulti = 1L,
    itmax = 1L,
    ngrid = 20L
  )

  expect_equal(enabled$fval, disabled$fval, tolerance = 1e-10)
  expect_gt(as.numeric(enabled$num.feval.fast[1L]), 0)
  expect_gt(as.numeric(enabled$num.feval.fast[1L]),
            as.numeric(disabled$num.feval.fast[1L]))
})
