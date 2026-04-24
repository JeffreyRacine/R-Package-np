library(np)

quiet_nomad_eval <- function(expr) {
  out <- NULL
  capture.output(out <- force(expr))
  out
}

expect_nomad_bw_timing <- function(bw) {
  expect_true(is.finite(bw$nomad.time))
  expect_true(is.finite(bw$powell.time))
  expect_true(is.finite(bw$total.time))
  expect_equal(as.double(bw$total.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
  expect_equal(as.double(bw$degree.search$optim.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
}

expect_nomad_fit_timing <- function(fit) {
  expect_true(is.finite(fit$bws$nomad.time))
  expect_true(is.finite(fit$bws$powell.time))
  expect_true(is.finite(fit$bws$total.time))
  expect_equal(as.double(fit$bws$total.time),
               as.double(fit$bws$nomad.time + fit$bws$powell.time),
               tolerance = 1e-8)

  if (!is.null(fit$optim.time) && length(fit$optim.time))
    expect_equal(as.double(fit$optim.time),
                 as.double(fit$bws$total.time),
                 tolerance = 1e-8)

  if (!is.null(fit$total.time) && length(fit$total.time) &&
      !is.null(fit$fit.time) && length(fit$fit.time) &&
      !is.null(fit$optim.time) && length(fit$optim.time))
    expect_equal(as.double(fit$total.time),
                 as.double(fit$optim.time + fit$fit.time),
                 tolerance = 1e-8)
}

test_that("NOMAD timing is consistent across supported serial families", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260401)
  d1 <- data.frame(x = sort(runif(12)))
  d1$y <- d1$x + 0.35 * d1$x^2 + rnorm(nrow(d1), sd = 0.12)

  d2 <- data.frame(x = runif(12), z = runif(12, -1, 1))
  d2$y <- 1 + d2$x + sin(pi * d2$z) + rnorm(nrow(d2), sd = 0.18)

  d3 <- data.frame(x1 = runif(12, -1, 1), x2 = runif(12, -1, 1))
  idx <- d3$x1 + 0.75 * d3$x2
  d3$y <- sin(idx) + 0.25 * idx^2 + rnorm(nrow(d3), sd = 0.08)

  bw_reg <- quiet_nomad_eval(npregbw(
    y ~ x,
    data = d1,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  ))
  fit_reg <- quiet_nomad_eval(npreg(y ~ x, data = d1, nomad = TRUE, degree.max = 1L, nmulti = 1L))
  expect_nomad_bw_timing(bw_reg)
  expect_nomad_fit_timing(fit_reg)

  bw_cdens <- quiet_nomad_eval(npcdensbw(
    y ~ x,
    data = d1,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  ))
  fit_cdens <- quiet_nomad_eval(npcdens(y ~ x, data = d1, nomad = TRUE, degree.max = 1L, nmulti = 1L))
  expect_nomad_bw_timing(bw_cdens)
  expect_nomad_fit_timing(fit_cdens)

  bw_cdist <- quiet_nomad_eval(npcdistbw(
    y ~ x,
    data = d1,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  ))
  fit_cdist <- quiet_nomad_eval(npcdist(y ~ x, data = d1, nomad = TRUE, degree.max = 1L, nmulti = 1L, ngrid = 7L))
  expect_nomad_bw_timing(bw_cdist)
  expect_nomad_fit_timing(fit_cdist)

  bw_pl <- quiet_nomad_eval(npplregbw(
    y ~ x | z,
    data = d2,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  ))
  fit_pl <- quiet_nomad_eval(npplreg(y ~ x | z, data = d2, nomad = TRUE, degree.max = 1L, nmulti = 1L))
  expect_nomad_bw_timing(bw_pl)
  expect_nomad_fit_timing(fit_pl)

  bw_sc <- quiet_nomad_eval(npscoefbw(
    y ~ x | z,
    data = d2,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  ))
  fit_sc <- quiet_nomad_eval(npscoef(y ~ x | z, data = d2, nomad = TRUE, degree.max = 1L, nmulti = 1L))
  expect_nomad_bw_timing(bw_sc)
  expect_nomad_fit_timing(fit_sc)

  bw_si <- quiet_nomad_eval(npindexbw(
    y ~ x1 + x2,
    data = d3,
    method = "ichimura",
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    nmulti = 1L
  ))
  fit_si <- quiet_nomad_eval(npindex(y ~ x1 + x2, data = d3, method = "ichimura", nomad = TRUE, degree.max = 1L, nmulti = 1L))
  expect_nomad_bw_timing(bw_si)
  expect_nomad_fit_timing(fit_si)
})
