test_that("npcdensbw direct nomad payload preserves CV metadata", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260320)
  dat <- data.frame(x = runif(24))
  dat$y <- rbeta(nrow(dat), 1, 1)

  bw <- np::npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L,
    cxkerbound = "range",
    cykerbound = "range"
  )

  expect_identical(bw$method, "cv.ls")
  expect_identical(bw$pmethod, "Least Squares Cross-Validation")
  expect_true(length(bw$ifval) == 1L && is.na(bw$ifval))
  expect_true(length(bw$fval.history) == 1L && is.na(bw$fval.history))
  expect_true(length(bw$eval.history) == 1L && is.na(bw$eval.history))
  expect_true(length(bw$invalid.history) == 1L && is.na(bw$invalid.history))
  expect_equal(
    bw$fval,
    np:::.npcdensbw_eval_only(data.frame(x = dat$x), data.frame(y = dat$y), bw)$objective,
    tolerance = 1e-12
  )

  printed <- paste(capture.output(print(bw)), collapse = "\n")
  expect_false(grepl("Manual", printed, fixed = TRUE))
  expect_match(printed, "achieved on multistart 1", fixed = TRUE)
})

test_that("npcdens direct nomad payload retains summary bandwidth metadata", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)
  suppressPackageStartupMessages(library(np))

  set.seed(42)
  n <- 500L
  x <- runif(n)
  y <- rbeta(n, 1 + x, 2 - x)

  fit <- np::npcdens(
    y ~ x,
    nomad = TRUE,
    search.engine = "nomad",
    cykerbound = "range"
  )

  expect_false(identical(fit$bws$sumNum, NA))
  expect_false(identical(fit$bws$bandwidth, NA))
  expect_false(identical(fit$bws$sfactor, NA))
  expect_gt(as.numeric(fit$bws$num.feval.fast[1L]), 0)

  summary_text <- expect_warning(capture.output(summary(fit$bws)), NA)
  summary_text <- paste(summary_text, collapse = "\n")
  expect_match(summary_text, "Exp\\. Var\\. Name:", perl = TRUE)
  expect_match(summary_text, "Dep\\. Var\\. Name:", perl = TRUE)
})

test_that("npcdens NOMAD accounting is owner-level on fast-path fits", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)
  suppressPackageStartupMessages(library(np))

  set.seed(42)
  n <- 250L
  x <- runif(n)
  y <- runif(n)

  fit <- np::npcdens(
    y ~ x,
    nomad = TRUE,
    search.engine = "nomad",
    cykerbound = "range"
  )

  expect_gt(as.numeric(fit$bws$num.feval.fast[1L]), 0)
  expect_lte(
    as.numeric(fit$bws$num.feval.fast[1L]),
    as.numeric(fit$bws$num.feval[1L])
  )

  ev <- np:::.npcdensbw_eval_only(
    data.frame(x = x),
    data.frame(y = y),
    fit$bws,
    invalid.penalty = "baseline",
    penalty.multiplier = 10
  )

  expect_equal(as.numeric(ev$num.feval), 1)
  expect_equal(as.numeric(ev$num.feval.fast), 1)
})
