library(np)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

expect_bw_timing_contract <- function(obj, summary_label = "Estimation Time:") {
  expect_true(is.finite(obj$total.time))
  txt <- paste(capture.output(summary(obj)), collapse = "\n")
  expect_match(txt, summary_label, fixed = TRUE)
  if (!is.null(obj$nomad.time) && is.finite(obj$nomad.time))
    expect_match(txt, "NOMAD ", fixed = TRUE)
  if (!is.null(obj$powell.time) && is.finite(obj$powell.time))
    expect_match(txt, "Powell ", fixed = TRUE)
}

expect_fit_timing_contract <- function(obj, summary_label = "Estimation Time:") {
  expect_true(is.finite(obj$total.time))
  expect_true(is.finite(obj$fit.time))
  expect_true(is.finite(obj$optim.time))
  expect_equal(obj$total.time, obj$optim.time + obj$fit.time, tolerance = 1e-8)

  txt <- paste(capture.output(summary(obj)), collapse = "\n")
  expect_match(txt, summary_label, fixed = TRUE)
  if (!is.null(obj$nomad.time) && is.finite(obj$nomad.time)) {
    expect_match(txt, "NOMAD ", fixed = TRUE)
    if (!is.null(obj$powell.time) && is.finite(obj$powell.time))
      expect_match(txt, "Powell ", fixed = TRUE)
    expect_match(txt, "fit ", fixed = TRUE)
  } else {
    expect_match(txt, "optim ", fixed = TRUE)
    expect_match(txt, "fit ", fixed = TRUE)
  }
}

test_that("direct density and distribution families report elapsed optimization time", {
  set.seed(42)
  n <- 120L
  x <- runif(n)
  y <- x^2 + rnorm(n, sd = 0.1)

  bw_udens <- quiet_eval(npudensbw(~x, nmulti = 1))
  fit_udens <- quiet_eval(npudens(bws = bw_udens))

  bw_udist <- quiet_eval(npudistbw(~x, nmulti = 1))
  fit_udist <- quiet_eval(npudist(bws = bw_udist))

  bw_cdens <- quiet_eval(npcdensbw(y ~ x, regtype = "ll", bwtype = "adaptive_nn", nmulti = 1))
  fit_cdens <- quiet_eval(npcdens(bws = bw_cdens, txdat = data.frame(x = x), tydat = data.frame(y = y)))

  bw_cdist <- quiet_eval(npcdistbw(y ~ x, regtype = "ll", bwtype = "adaptive_nn", nmulti = 1))
  fit_cdist <- quiet_eval(npcdist(bws = bw_cdist, txdat = data.frame(x = x), tydat = data.frame(y = y)))

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

  bw_reg <- quiet_eval(npregbw(xdat = data.frame(x = x), ydat = y, regtype = "ll", bwtype = "adaptive_nn", nmulti = 1))
  fit_reg <- quiet_eval(npreg(bws = bw_reg, txdat = data.frame(x = x), tydat = y))

  bw_qreg <- quiet_eval(npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    regtype = "ll",
    bwtype = "adaptive_nn",
    nmulti = 1
  ))
  fit_qreg <- quiet_eval(npqreg(
    bws = bw_qreg,
    tau = 0.5
  ))

  expect_bw_timing_contract(bw_reg)
  expect_bw_timing_contract(bw_qreg)
  expect_fit_timing_contract(fit_reg)
  expect_fit_timing_contract(fit_qreg)
})

test_that("composite core families remain timing-contract compliant", {
  set.seed(123)

  n <- 50L
  x <- runif(n)
  z <- runif(n)
  y <- x^2 + 2 * z + rnorm(n, sd = 0.1)
  bw_pl <- quiet_eval(npplregbw(xdat = z, zdat = x, ydat = y, regtype = "ll", nmulti = 1))
  fit_pl <- quiet_eval(npplreg(bws = bw_pl))

  x1 <- runif(40L)
  x2 <- runif(40L)
  y_si <- (x1 + x2)^2 + rnorm(40L, sd = 0.1)
  bw_si <- quiet_eval(npindexbw(xdat = data.frame(x1 = x1, x2 = x2), ydat = y_si, method = "ichimura", nmulti = 1))
  fit_si <- quiet_eval(npindex(bws = bw_si))

  xs <- runif(n)
  zs <- runif(n)
  y_sc <- (0.5 + xs^2) * zs + rnorm(n, sd = 0.1)
  bw_sc <- quiet_eval(npscoefbw(xdat = xs, zdat = zs, ydat = y_sc, regtype = "ll", nmulti = 1))
  fit_sc <- quiet_eval(npscoef(bws = bw_sc))

  for (obj in list(bw_pl, bw_si, bw_sc))
    expect_bw_timing_contract(obj)

  for (obj in list(fit_pl, fit_si, fit_sc))
    expect_fit_timing_contract(obj)
})

test_that("NOMAD plus Powell timing is carried through bandwidth and fit summaries", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(24)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  bw <- quiet_eval(npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  ))
  fit <- quiet_eval(npreg(bws = bw, txdat = data.frame(x = dat$x), tydat = dat$y))

  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(is.finite(bw$nomad.time))
  expect_true(is.finite(bw$powell.time))
  expect_true(is.finite(bw$total.time))
  expect_gte(bw$total.time + 1e-8, bw$nomad.time + bw$powell.time)
  expect_equal(fit$optim.time, bw$total.time, tolerance = 1e-8)
  expect_equal(fit$nomad.time, bw$nomad.time, tolerance = 1e-8)
  expect_equal(fit$powell.time, bw$powell.time, tolerance = 1e-8)

  bw_txt <- paste(capture.output(summary(bw)), collapse = "\n")
  fit_txt <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(bw_txt, "NOMAD ", fixed = TRUE)
  expect_match(bw_txt, "Powell ", fixed = TRUE)
  expect_match(fit_txt, "NOMAD ", fixed = TRUE)
  expect_match(fit_txt, "Powell ", fixed = TRUE)
  expect_match(fit_txt, "fit ", fixed = TRUE)
})
