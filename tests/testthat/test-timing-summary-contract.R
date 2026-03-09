library(np)

expect_bw_timing_contract <- function(obj, summary_label = "Estimation Time:") {
  expect_true(is.finite(obj$total.time))
  txt <- paste(capture.output(summary(obj)), collapse = "\n")
  expect_match(txt, summary_label, fixed = TRUE)
}

expect_fit_timing_contract <- function(obj, summary_label = "Estimation Time:") {
  expect_true(is.finite(obj$total.time))
  expect_true(is.finite(obj$fit.time))
  expect_true(is.finite(obj$optim.time))
  expect_equal(obj$total.time, obj$optim.time + obj$fit.time, tolerance = 1e-8)

  txt <- paste(capture.output(summary(obj)), collapse = "\n")
  expect_match(txt, summary_label, fixed = TRUE)
  expect_match(txt, "optim ", fixed = TRUE)
  expect_match(txt, "fit ", fixed = TRUE)
}

test_that("direct density and distribution families report elapsed optimization time", {
  set.seed(42)
  n <- 120L
  x <- runif(n)
  y <- x^2 + rnorm(n, sd = 0.1)

  bw_udens <- npudensbw(~x, nmulti = 1)
  fit_udens <- npudens(bws = bw_udens)

  bw_udist <- npudistbw(~x, nmulti = 1)
  fit_udist <- npudist(bws = bw_udist)

  bw_cdens <- npcdensbw(y ~ x, regtype = "ll", bwtype = "adaptive_nn", nmulti = 1)
  fit_cdens <- npcdens(bws = bw_cdens, txdat = data.frame(x = x), tydat = data.frame(y = y))

  bw_cdist <- npcdistbw(y ~ x, regtype = "ll", bwtype = "adaptive_nn", nmulti = 1)
  fit_cdist <- npcdist(bws = bw_cdist, txdat = data.frame(x = x), tydat = data.frame(y = y))

  for (obj in list(bw_udens, bw_udist, bw_cdens, bw_cdist))
    expect_bw_timing_contract(obj)

  for (obj in list(fit_udens, fit_udist, fit_cdens, fit_cdist))
    expect_fit_timing_contract(obj)

  expect_gt(fit_cdens$optim.time, 0)
  expect_gt(fit_cdist$optim.time, 0)
})

test_that("regression and quantile regression remain timing-contract compliant", {
  set.seed(99)
  n <- 100L
  x <- runif(n)
  y <- x^2 + rnorm(n, sd = 0.15)

  bw_reg <- npregbw(xdat = data.frame(x = x), ydat = y, regtype = "ll", bwtype = "adaptive_nn", nmulti = 1)
  fit_reg <- npreg(bws = bw_reg, txdat = data.frame(x = x), tydat = y)

  bw_qreg <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    regtype = "ll",
    bwtype = "adaptive_nn",
    nmulti = 1
  )
  fit_qreg <- npqreg(
    bws = bw_qreg,
    tau = 0.5
  )

  expect_bw_timing_contract(bw_reg)
  expect_bw_timing_contract(bw_qreg)
  expect_fit_timing_contract(fit_reg)
  expect_fit_timing_contract(fit_qreg)
})
