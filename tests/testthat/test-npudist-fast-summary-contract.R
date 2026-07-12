test_that("npudist exposes unconditional fast counts in bandwidth summary", {
  library(np)
  set.seed(42)
  n <- 100L
  y <- runif(n)

  fit <- npudist(
    ~ y,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  expect_true(is.finite(as.numeric(fit$bws$num.feval[1L])))
  expect_true(is.finite(as.numeric(fit$bws$num.feval.fast[1L])))
  expect_gt(as.numeric(fit$bws$num.feval.fast[1L]), 0)

  txt <- paste(capture.output(summary(fit$bws)), collapse = "\n")
  expect_true(grepl("Number of Function Evaluations:", txt, fixed = TRUE))

  fast <- as.numeric(fit$bws$num.feval.fast[1L])
  cache.hits <- if (!is.null(fit$bws$nn.cache) &&
                    "objective.hits" %in% names(fit$bws$nn.cache)) {
    as.numeric(fit$bws$nn.cache[["objective.hits"]])
  } else {
    0
  }

  if (is.finite(cache.hits) && cache.hits > 0)
    expect_true(grepl("Evaluation cache (Powell):", txt, fixed = TRUE))

  fast.extra <- fast - if (is.finite(cache.hits)) cache.hits else 0
  expect_gte(fast.extra, 0)
  expect_lte(fast.extra, as.numeric(fit$bws$num.feval[1L]))
  if (is.finite(fast.extra) && fast.extra > 0) {
    expect_true(grepl(
      sprintf("Fast CV route: %s of", format(round(fast.extra), big.mark = ",")),
      txt,
      fixed = TRUE
    ))
  } else {
    expect_false(grepl("Fast CV route:", txt, fixed = TRUE))
  }
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
  cache.hits <- if (!is.null(forced$nn.cache) &&
                    "objective.hits" %in% names(forced$nn.cache)) {
    as.numeric(forced$nn.cache[["objective.hits"]])
  } else {
    0
  }
  fast.computations <- as.numeric(forced$num.feval.fast[1L]) - cache.hits
  expect_gte(fast.computations, 0)
  expect_lte(fast.computations, as.numeric(forced$num.feval[1L]))
})
