.beta_search_half_bounds <- function(x) {
  u <- sort(unique(as.double(x)))
  c(
    lower = u[[1L]] - (u[[2L]] - u[[1L]]) / 2,
    upper = u[[length(u)]] +
      (u[[length(u)]] - u[[length(u) - 1L]]) / 2
  )
}

test_that("beta range density CVLS certifies the asymptotic candidate", {
  set.seed(4010007)
  x <- sort(runif(250))
  bounds <- .beta_search_half_bounds(x)

  ordinary <- npudensbw(
    ~x, ckertype = "beta", ckerbound = "fixed",
    ckerlb = bounds[["lower"]], ckerub = bounds[["upper"]],
    bwmethod = "cv.ls", nmulti = 1
  )
  certified <- npudensbw(
    ~x, ckertype = "beta", ckerbound = "range",
    bwmethod = "cv.ls", nmulti = 1
  )

  expect_null(ordinary$beta.range.certification)
  expect_true(certified$fval > ordinary$fval)
  expect_true(certified$beta.range.certification$triggered)
  expect_identical(certified$beta.range.certification$selected, "asymptotic")
  expect_equal(certified$beta.range.certification$concentration,
               .Machine$double.eps, tolerance = 0)
  expect_equal(certified$beta.range.certification$tail.metric,
               unname(diff(bounds)) * 2^26, tolerance = 0)
  expect_true(certified$num.feval.certified > 0)

  evaluated <- getS3method("npudensbw", "bandwidth")(
    dat = data.frame(x = x), bws = certified,
    bandwidth.compute = TRUE, eval.only = TRUE, nmulti = 1
  )
  expect_equal(certified$fval, evaluated$fval, tolerance = 2e-12)
  expect_match(paste(capture.output(summary(certified)), collapse = "\n"),
               "Beta-range certification evaluations")
})

test_that("beta range distribution CDF uses the support-width bridge", {
  expected <- list(
    `5010001` = c(h = 4.45267153258119, objective = 0.165866148831325),
    `5010002` = c(h = 1.71948612308104, objective = 0.166450382167972)
  )

  for (seed in as.integer(names(expected))) {
    set.seed(seed)
    x <- sort(runif(250))
    certified <- npudistbw(
      ~x, ckertype = "beta", ckerbound = "range",
      bwmethod = "cv.cdf", nmulti = 1
    )
    target <- expected[[as.character(seed)]]
    expect_true(certified$beta.range.certification$triggered)
    expect_identical(certified$beta.range.certification$selected,
                     "refinement")
    expect_equal(certified$bw, target[["h"]], tolerance = 2e-7)
    expect_equal(certified$fval, target[["objective"]], tolerance = 2e-12)
    expect_true(certified$num.feval.certified > 0)

    evaluated <- getS3method("npudistbw", "dbandwidth")(
      dat = data.frame(x = x), bws = certified,
      bandwidth.compute = TRUE, eval.only = TRUE, nmulti = 1
    )
    expect_equal(certified$fval, evaluated$fval, tolerance = 2e-12)
  }
})

test_that("beta range certification retains an ordinary winner", {
  set.seed(4010001)
  x <- sort(runif(250))
  certified <- npudensbw(
    ~x, ckertype = "beta", ckerbound = "range",
    bwmethod = "cv.ls", nmulti = 1
  )
  expect_false(certified$beta.range.certification$triggered)
  expect_identical(certified$beta.range.certification$selected, "ordinary")
  expect_true(certified$num.feval.certified > 0)
})

test_that("bandwidth-object restarts receive exactly one certification", {
  set.seed(4010007)
  x <- sort(runif(250))
  density.start <- npudensbw(
    ~x, bws = 0.25, bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "range", bwmethod = "cv.ls"
  )
  density <- npudensbw(
    dat = data.frame(x = x), bws = density.start,
    bandwidth.compute = TRUE, nmulti = 1
  )
  expect_true(density$beta.range.certification$triggered)
  expect_identical(density$beta.range.certification$selected, "asymptotic")
  expect_equal(
    density$bw, density$beta.range.certification$tail.metric,
    tolerance = 2e-8
  )

  set.seed(5010001)
  x <- sort(runif(250))
  distribution.start <- npudistbw(
    ~x, bws = 0.25, bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "range", bwmethod = "cv.cdf"
  )
  distribution <- npudistbw(
    dat = data.frame(x = x), bws = distribution.start,
    bandwidth.compute = TRUE, nmulti = 1
  )
  expect_true(distribution$beta.range.certification$triggered)
  expect_identical(
    distribution$beta.range.certification$selected, "refinement"
  )
  expect_equal(distribution$bw, 4.45267153258119, tolerance = 2e-7)
})

test_that("beta range certification is isolated from adjacent routes", {
  set.seed(4010007)
  x <- sort(runif(250))
  bounds <- .beta_search_half_bounds(x)

  fixed <- npudensbw(
    ~x, ckertype = "beta", ckerbound = "fixed",
    ckerlb = bounds[["lower"]], ckerub = bounds[["upper"]],
    bwmethod = "cv.ls", nmulti = 1
  )
  cvml <- npudensbw(
    ~x, ckertype = "beta", ckerbound = "range",
    bwmethod = "cv.ml", nmulti = 1
  )
  higher <- npudensbw(
    ~x, ckertype = "beta", ckerorder = 4, ckerbound = "range",
    bwmethod = "cv.ls", nmulti = 1
  )
  gaussian <- npudensbw(
    ~x, ckertype = "gaussian", ckerbound = "range",
    bwmethod = "cv.ls", nmulti = 1
  )
  manual <- npudensbw(
    ~x, bws = 0.2, bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "range", bwmethod = "cv.ls"
  )

  expect_null(fixed$beta.range.certification)
  expect_null(cvml$beta.range.certification)
  expect_null(higher$beta.range.certification)
  expect_null(gaussian$beta.range.certification)
  expect_null(manual$beta.range.certification)
})

test_that("asymptotic candidate conversion respects bandwidth scaling", {
  x <- data.frame(x = seq(0.05, 0.95, length.out = 40))
  dati <- np:::untangle(x)
  bws <- np:::bandwidth(
    bw = 0.5, bwmethod = "cv.ls", bwscaling = TRUE,
    ckertype = "beta", ckerbound = "range", xdati = dati,
    xnames = "x", nobs = nrow(x), sdev = sd(x$x),
    nconfac = nrow(x)^(-1 / 5), bandwidth.compute = TRUE
  )
  converted <- np:::npBetaRangeCertificationBandwidths(
    bws, where = "test"
  )
  width <- bws$ckerub[bws$icon] - bws$ckerlb[bws$icon]
  factor <- sd(x$x) * nrow(x)^(-1 / 5)
  expect_equal(converted$bridge[bws$icon] * factor, width,
               tolerance = 2e-15)
  expect_equal(converted$tail[bws$icon] * factor, width * 2^26,
               tolerance = 2e-8)
})

test_that("automatic certification respects public scale-factor storage", {
  set.seed(4010007)
  x <- sort(runif(250))
  certified <- npudensbw(
    ~x, ckertype = "beta", ckerbound = "range",
    bwmethod = "cv.ls", bwscaling = TRUE, nmulti = 1
  )
  metric <- certified$bw * certified$sdev * certified$nconfac

  expect_true(certified$beta.range.certification$triggered)
  expect_identical(certified$beta.range.certification$selected, "asymptotic")
  expect_equal(
    metric, certified$beta.range.certification$tail.metric,
    tolerance = 2e-8
  )
})
