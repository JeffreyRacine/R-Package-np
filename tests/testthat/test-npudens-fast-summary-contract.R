test_that("npudens exposes unconditional fast counts in bandwidth summary", {
  library(np)
  set.seed(42)
  n <- 100L
  y <- runif(n)

  fit <- npudens(
    ~ y,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1,
    bwsolver = "mads",
    nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 40L)
  )

  expect_true(is.finite(as.numeric(fit$bws$num.feval[1L])))
  expect_true(is.finite(as.numeric(fit$bws$num.feval.fast[1L])))
  expect_gt(as.numeric(fit$bws$num.feval.fast[1L]), 0)
  expect_lte(as.numeric(fit$bws$num.feval.fast[1L]), as.numeric(fit$bws$num.feval[1L]))

  txt <- paste(capture.output(summary(fit$bws)), collapse = "\n")
  expect_true(grepl("Number of Function Evaluations:", txt, fixed = TRUE))
  expect_match(txt, sprintf("Fast CV route:\\s+%s of", format(fit$bws$num.feval.fast[1L])))
})
