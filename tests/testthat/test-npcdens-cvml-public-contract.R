library(np)

public_shadow_empty <- function(n) {
  matrix(numeric(0), nrow = n, ncol = 0)
}

public_shadow_cker <- function(kernel) {
  switch(kernel,
    gaussian = 0L,
    epanechnikov = 4L,
    uniform = 8L,
    truncated = 9L,
    stop("unsupported continuous kernel")
  )
}

public_shadow_uker <- function(kernel) {
  switch(kernel,
    aitchisonaitken = 0L,
    liracine = 1L,
    stop("unsupported unordered kernel")
  )
}

public_shadow_oker <- function(kernel) {
  switch(kernel,
    wangvanryzin = 0L,
    liracine = 2L,
    racineliyan = 3L,
    stop("unsupported ordered kernel")
  )
}

public_shadow_rbw <- function(bw) {
  c(
    bw$xbw[bw$ixcon],
    bw$ybw[bw$iycon],
    bw$ybw[bw$iyuno],
    bw$ybw[bw$iyord],
    bw$xbw[bw$ixuno],
    bw$xbw[bw$ixord]
  )
}

public_shadow_bwtype <- function(bw) {
  switch(bw$type,
    fixed = 0L,
    generalized_nn = 1L,
    adaptive_nn = 2L,
    stop("unsupported bandwidth type")
  )
}

call_public_cvml_shadow <- function(bw, x, y, tree = FALSE) {
  n <- nrow(x)
  .Call(
    "C_np_shadow_cv_density_conditional",
    public_shadow_empty(n), public_shadow_empty(n), as.matrix(y),
    public_shadow_empty(n), public_shadow_empty(n), as.matrix(x),
    as.double(public_shadow_rbw(bw)),
    public_shadow_bwtype(bw),
    public_shadow_cker(bw$cykertype),
    public_shadow_uker(bw$uykertype),
    public_shadow_oker(bw$oykertype),
    public_shadow_cker(bw$cxkertype),
    public_shadow_uker(bw$uxkertype),
    public_shadow_oker(bw$oxkertype),
    FALSE,
    0L,
    0L,
    integer(0),
    tree,
    0L,
    TRUE,
    PACKAGE = "np"
  )
}

test_that("public npcdensbw cv.ml keeps lc on the legacy objective", {
  set.seed(202)
  n <- 32L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))

  bw.lc <- npcdensbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ml", nmulti = 0)
  shadow <- call_public_cvml_shadow(bw.lc, x, y)

  expect_equal(-bw.lc$fval, shadow$old, tolerance = 1e-10)
})

test_that("public npcdensbw cv.ml enforces ll == lp(glp, degree=1)", {
  set.seed(101)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 - x$x2 + rnorm(n, sd = 0.2))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(xdat = x, ydat = y, regtype = "ll", bwmethod = "cv.ml", nmulti = 0)
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ml",
    nmulti = 0
  )

  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-12)
  expect_equal(bw.ll$xbw, bw.lp$xbw, tolerance = 1e-12)
  expect_equal(bw.ll$ybw, bw.lp$ybw, tolerance = 1e-12)
})

test_that("public npcdensbw cv.ml preserves tree parity at fixed LP bandwidths", {
  old <- options(np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(111)
  n <- 34L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  options(np.tree = FALSE)
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ml",
    nmulti = 0
  )
  shadow.nt <- call_public_cvml_shadow(bw.lp, x, y, FALSE)
  shadow.tr <- call_public_cvml_shadow(bw.lp, x, y, TRUE)

  expect_equal(shadow.nt$old, shadow.tr$old, tolerance = 1e-10)
  expect_equal(shadow.nt$new, shadow.tr$new, tolerance = 1e-10)
})

test_that("public npcdensbw cv.ml generalized-nn is finite and preserves ll canonicalization", {
  set.seed(121)
  n <- 40L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.12))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 0
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ml",
    nmulti = 0
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-10)
  expect_equal(bw.ll$xbw, bw.lp$xbw, tolerance = 1e-10)
  expect_equal(bw.ll$ybw, bw.lp$ybw, tolerance = 1e-10)
})

test_that("public npcdens cv.ml LP route does not collapse gaussian bandwidths", {
  set.seed(42)
  n <- 100L
  x <- data.frame(x = rnorm(n))
  y <- data.frame(y = x$x + rnorm(n))

  bw.lc <- npcdensbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ml", nmulti = 0)
  bw.ll <- npcdensbw(xdat = x, ydat = y, regtype = "ll", bwmethod = "cv.ml", nmulti = 0)
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = 3L,
    bwmethod = "cv.ml",
    nmulti = 0
  )

  fit.ll <- npcdens(bws = bw.ll)
  fit.lp <- npcdens(bws = bw.lp)

  expect_true(all(is.finite(fitted(fit.ll))))
  expect_true(all(is.finite(fitted(fit.lp))))
  expect_gte(min(fitted(fit.ll)), 0)
  expect_gte(min(fitted(fit.lp)), 0)
  expect_lt(max(fitted(fit.ll)), 10)
  expect_lt(max(fitted(fit.lp)), 10)
  expect_gt(unname(bw.ll$bandwidth$y), 1e-4)
  expect_gt(unname(bw.lp$bandwidth$y), 1e-4)
  expect_lt(max(fitted(fit.ll)) / max(fitted(npcdens(bws = bw.lc))), 20)
  expect_lt(max(fitted(fit.lp)) / max(fitted(npcdens(bws = bw.lc))), 20)
})
