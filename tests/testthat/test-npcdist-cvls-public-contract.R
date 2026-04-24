library(np)

cdist_shadow_empty <- function(n) {
  matrix(numeric(0), nrow = n, ncol = 0)
}

cdist_shadow_cker <- function(kernel) {
  switch(kernel,
    gaussian = 0L,
    epanechnikov = 4L,
    uniform = 8L,
    truncated = 9L,
    stop("unsupported continuous kernel")
  )
}

cdist_shadow_uker <- function(kernel) {
  switch(kernel,
    aitchisonaitken = 0L,
    liracine = 1L,
    stop("unsupported unordered kernel")
  )
}

cdist_shadow_oker <- function(kernel) {
  switch(kernel,
    wangvanryzin = 0L,
    liracine = 2L,
    racineliyan = 3L,
    stop("unsupported ordered kernel")
  )
}

cdist_shadow_rbw <- function(bw) {
  c(
    bw$xbw[bw$ixcon],
    bw$ybw[bw$iycon],
    bw$ybw[bw$iyuno],
    bw$ybw[bw$iyord],
    bw$xbw[bw$ixuno],
    bw$xbw[bw$ixord]
  )
}

cdist_shadow_bwtype <- function(bw) {
  switch(bw$type,
    fixed = 0L,
    generalized_nn = 1L,
    adaptive_nn = 2L,
    stop("unsupported bandwidth type")
  )
}

cdist_shadow_safe_call <- function(name, ...) {
  on.exit(
    tryCatch(.Call("C_np_shadow_reset_state", PACKAGE = "np"),
             error = function(e) NULL),
    add = TRUE
  )
  .Call(name, ..., PACKAGE = "np")
}

call_public_cdist_cvls_shadow <- function(bw, x, ytrain, yeval = ytrain, cdfontrain = FALSE) {
  n <- nrow(x)
  ne <- nrow(yeval)
  cdist_shadow_safe_call(
    "C_np_shadow_cv_distribution_conditional",
    cdist_shadow_empty(n), cdist_shadow_empty(n), as.matrix(ytrain),
    cdist_shadow_empty(ne), cdist_shadow_empty(ne), as.matrix(yeval),
    cdist_shadow_empty(n), cdist_shadow_empty(n), as.matrix(x),
    as.double(cdist_shadow_rbw(bw)),
    cdist_shadow_bwtype(bw),
    cdist_shadow_cker(bw$cykertype),
    cdist_shadow_uker(bw$uykertype),
    cdist_shadow_oker(bw$oykertype),
    cdist_shadow_cker(bw$cxkertype),
    cdist_shadow_uker(bw$uxkertype),
    cdist_shadow_oker(bw$oxkertype),
    FALSE,
    0L,
    integer(0),
    FALSE,
    0L,
    cdfontrain,
    TRUE
  )
}

public_cdist_eval_grid <- function(ydat, ngrid = 100L) {
  probs <- seq(0, 1, length.out = ngrid)
  evy <- ydat[seq_len(ngrid), , drop = FALSE]
  for (i in seq_len(ncol(evy))) {
    evy[, i] <- uocquantile(ydat[, i], probs)
  }
  evy
}

test_that("public npcdistbw cv.ls keeps lc on the legacy objective", {
  set.seed(303)
  n <- 32L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))

  bw.lc <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  shadow <- call_public_cdist_cvls_shadow(bw.lc, x, y, yeval = public_cdist_eval_grid(y))

  expect_equal(bw.lc$fval, shadow$old, tolerance = 1e-10)
})

test_that("public npcdistbw cv.ls fixed LP/LL route activates with ll == lp parity", {
  set.seed(20260311)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )

  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdistbw cv.ls generalized-nn LP route activates with ll == lp parity", {
  set.seed(20260312)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdistbw cv.ls adaptive-nn LP route activates with ll == lp parity", {
  set.seed(20260313)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.lp <- npcdistbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})
