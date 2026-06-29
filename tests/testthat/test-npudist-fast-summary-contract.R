test_that("npudist exposes unconditional fast counts in bandwidth summary", {
  library(np)
  set.seed(42)
  n <- 100L
  y <- runif(n)

  fit <- npudist(
    ~ y,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    bwsolver = "mads+powell",
    nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 20L),
    ngrid = 30L
  )

  expect_true(is.finite(as.numeric(fit$bws$num.feval[1L])))
  expect_true(is.finite(as.numeric(fit$bws$num.feval.fast[1L])))
  expect_gt(as.numeric(fit$bws$num.feval.fast[1L]), 0)
  expect_lte(as.numeric(fit$bws$num.feval.fast[1L]), as.numeric(fit$bws$num.feval[1L]))

  txt <- paste(capture.output(summary(fit$bws)), collapse = "\n")
  expect_true(grepl("Number of Function Evaluations:", txt, fixed = TRUE))
  fast.line <- sprintf("Fast CV route:\\s+%s of", format(fit$bws$num.feval.fast[1L]))
  cache.line <- sprintf(
    "Powell cache: %s repeated objective lookups avoided",
    format(fit$bws$num.feval.fast[1L])
  )
  expect_true(grepl(fast.line, txt) || grepl(cache.line, txt, fixed = TRUE))
})

test_that("npudist counts forced large-h unconditional evaluations as fast", {
  set.seed(42)
  y <- runif(100L)
  dat <- data.frame(y = y)

  bw0 <- npudistbw(
    dat = dat,
    bwmethod = "cv.cdf",
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    bandwidth.compute = FALSE
  )

  bw_big <- bw0
  bw_big$bw[] <- 1e10
  bw_big$bandwidth[] <- 1e10
  bw_big$sfactor[] <- 1e10

  forced <- npudistbw(
    dat = dat,
    bws = bw_big,
    bandwidth.compute = TRUE,
    nmulti = 1L,
    itmax = 1L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  expect_true(is.finite(as.numeric(forced$num.feval[1L])))
  expect_true(is.finite(as.numeric(forced$num.feval.fast[1L])))
  expect_gt(as.numeric(forced$num.feval.fast[1L]), 0)
  expect_lte(as.numeric(forced$num.feval.fast[1L]), as.numeric(forced$num.feval[1L]))
})
