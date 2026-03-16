test_that("npudens default constructs bandwidth call in caller frame", {
  run_case <- function() {
    set.seed(101)
    x <- runif(40)
    npudens(tdat = x, bwmethod = "normal-reference")
  }

  fit <- run_case()
  expect_s3_class(fit, "npdensity")
  expect_true(any(grepl("bandwidth$", class(fit$bws))))
  expect_equal(NROW(fit$dens), NROW(fit$eval))
})

test_that("npudist default constructs bandwidth call in caller frame", {
  run_case <- function() {
    set.seed(202)
    x <- runif(45)
    npudist(tdat = x, bwmethod = "normal-reference")
  }

  fit <- run_case()
  expect_s3_class(fit, "npdistribution")
  expect_true(any(grepl("bandwidth$", class(fit$bws))))
  expect_equal(length(fit$dist), NROW(fit$eval))
})

test_that("npreg default constructs bandwidth call in caller frame", {
  run_case <- function() {
    set.seed(303)
    x <- runif(55)
    y <- sin(2 * pi * x) + rnorm(55, sd = 0.2)
    npreg(txdat = x, tydat = y, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  }

  fit <- run_case()
  expect_s3_class(fit, "npregression")
  expect_true(any(grepl("bandwidth$", class(fit$bws))))
  expect_equal(length(fit$mean), NROW(fit$eval))
})

test_that("npreg default bandwidth selection remains finite after bounded conditional selectors", {
  set.seed(20260316)
  n <- 48
  x <- runif(n)
  y <- runif(n)

  bw.cd <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cxkerorder = 2L,
    cykerorder = 2L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )

  expect_true(all(is.finite(as.numeric(bw.cd$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cd$ybw))))

  fit.reg <- npreg(
    txdat = x,
    tydat = sin(2 * pi * x) + rnorm(n, sd = 0.15),
    regtype = "lc",
    bwmethod = "cv.ls",
    nmulti = 1
  )

  expect_s3_class(fit.reg, "npregression")
  expect_true(any(grepl("bandwidth$", class(fit.reg$bws))))
  expect_equal(length(fit.reg$mean), NROW(fit.reg$eval))
})
