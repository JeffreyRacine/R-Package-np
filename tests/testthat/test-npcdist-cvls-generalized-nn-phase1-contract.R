library(np)

phase1_npcdist_cvls_gnn_fixture <- function() {
  set.seed(20260309)
  n <- 24L
  x <- data.frame(
    x1 = runif(n),
    x2 = runif(n)
  )
  y <- data.frame(
    y1 = x$x1^2 - 0.35 * x$x2 + 0.2 * sin(2 * pi * x$x1) + rnorm(n, sd = 0.08)
  )
  list(x = x, y = y)
}

phase1_npcdist_cvls_gnn_cases <- local({
  cache <- NULL

  function() {
    if (!is.null(cache))
      return(cache)

    dat <- phase1_npcdist_cvls_gnn_fixture()
    degree1 <- rep.int(1L, ncol(dat$x))
    degree2 <- rep.int(2L, ncol(dat$x))

    cache <<- list(
      dat = dat,
      degree1 = degree1,
      degree2 = degree2,
      bw.lc = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "lc",
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      ),
      bw.ll = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "ll",
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      ),
      bw.lp = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "lp",
        basis = "glp",
        degree = degree1,
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      ),
      bw.d2 = npcdistbw(
        xdat = dat$x,
        ydat = dat$y,
        regtype = "lp",
        basis = "glp",
        degree = degree2,
        bwtype = "generalized_nn",
        bwmethod = "cv.ls",
        nmulti = 1L,
        itmax = 1L
      )
    )

    cache
  }
})

test_that("phase1 npcdistbw cv.ls generalized-nn lc matches the frozen public baseline", {
  bw.lc <- phase1_npcdist_cvls_gnn_cases()$bw.lc

  expect_true(is.finite(bw.lc$fval))
  expect_equal(bw.lc$fval, 0.1096745754682811, tolerance = 1e-10)
})

test_that("phase1 npcdistbw cv.ls generalized-nn keeps ll on canonical lp degree-1 glp", {
  cases <- phase1_npcdist_cvls_gnn_cases()
  degree <- cases$degree1
  bw.ll <- cases$bw.ll
  bw.lp <- cases$bw.lp

  expect_identical(bw.ll$regtype.engine, "lp")
  expect_identical(bw.ll$basis.engine, "glp")
  expect_identical(as.integer(bw.ll$degree.engine), degree)
  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, 0.10210193919818145, tolerance = 1e-10)
  expect_equal(bw.lp$fval, 0.10210193919818145, tolerance = 1e-10)
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-10)
})

test_that("phase1 npcdistbw cv.ls generalized-nn lp degree-2 succeeds on a higher-order fixture", {
  cases <- phase1_npcdist_cvls_gnn_cases()
  bw.d1 <- cases$bw.lp
  bw.d2 <- cases$bw.d2

  expect_identical(as.integer(bw.d2$degree.engine), cases$degree2)
  expect_true(is.finite(bw.d2$fval))
  expect_equal(bw.d2$fval, 0.093788759039826752, tolerance = 1e-10)
  expect_gt(abs(bw.d2$fval - bw.d1$fval), 1e-6)
})

test_that("phase1 npcdistbw cv.ls generalized-nn avoids search-boundary collapse", {
  set.seed(20260309)
  n <- 30L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = rbeta(n, shape1 = 1 + 2 * x$x1, shape2 = 2))

  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1L,
    itmax = 1L
  )
  fit <- npcdist(
    bws = bw,
    exdat = data.frame(x1 = rep(0.5, 4L)),
    eydat = data.frame(y1 = c(0, 0.02, 0.98, 1))
  )
  pred <- fitted(fit)

  expect_true(is.finite(bw$fval))
  expect_true(all(is.finite(pred)))
  expect_true(all(pred >= -1e-6))
  expect_true(all(pred <= 1 + 1e-6))
  expect_gt(pred[3] - pred[1], 0.5)
  expect_true(all(bw$xbw > 1))
  expect_true(all(bw$ybw > 1))
  expect_true(all(bw$xbw < n))
  expect_true(all(bw$ybw < n))
})
