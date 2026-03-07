library(np)

cvls_shadow_empty <- function(n) {
  matrix(numeric(0), nrow = n, ncol = 0)
}

cvls_shadow_cker <- function(kernel) {
  switch(kernel,
    gaussian = 0L,
    epanechnikov = 4L,
    uniform = 8L,
    truncated = 9L,
    stop("unsupported continuous kernel")
  )
}

cvls_shadow_uker <- function(kernel) {
  switch(kernel,
    aitchisonaitken = 0L,
    liracine = 1L,
    stop("unsupported unordered kernel")
  )
}

cvls_shadow_oker <- function(kernel) {
  switch(kernel,
    wangvanryzin = 0L,
    liracine = 2L,
    racineliyan = 3L,
    stop("unsupported ordered kernel")
  )
}

cvls_shadow_rbw <- function(bw) {
  c(
    bw$xbw[bw$ixcon],
    bw$ybw[bw$iycon],
    bw$ybw[bw$iyuno],
    bw$ybw[bw$iyord],
    bw$xbw[bw$ixuno],
    bw$xbw[bw$ixord]
  )
}

cvls_shadow_bwtype <- function(bw) {
  switch(bw$type,
    fixed = 0L,
    generalized_nn = 1L,
    adaptive_nn = 2L,
    stop("unsupported bandwidth type")
  )
}

cvls_shadow_regtype <- function(bw) {
  if (identical(bw$regtype.engine, "lp")) 2L else 0L
}

cvls_shadow_degree <- function(bw) {
  if (identical(bw$regtype.engine, "lp")) as.integer(bw$degree.engine) else integer(0)
}

cvls_shadow_basis <- function(basis_engine, regtype_engine) {
  if (!identical(regtype_engine, "lp")) {
    return(0L)
  }
  switch(basis_engine,
    additive = 0L,
    glp = 1L,
    tensor = 2L,
    stop("unsupported LP basis")
  )
}

call_public_cvls_shadow <- function(bw, x, y, compare_old = identical(bw$regtype.engine, "lc")) {
  n <- nrow(x)
  .Call(
    "C_np_shadow_cv_density_conditional",
    cvls_shadow_empty(n), cvls_shadow_empty(n), as.matrix(y),
    cvls_shadow_empty(n), cvls_shadow_empty(n), as.matrix(x),
    as.double(cvls_shadow_rbw(bw)),
    cvls_shadow_bwtype(bw),
    cvls_shadow_cker(bw$cykertype),
    cvls_shadow_uker(bw$uykertype),
    cvls_shadow_oker(bw$oykertype),
    cvls_shadow_cker(bw$cxkertype),
    cvls_shadow_uker(bw$uxkertype),
    cvls_shadow_oker(bw$oxkertype),
    FALSE,
    1L,
    cvls_shadow_regtype(bw),
    cvls_shadow_degree(bw),
    isTRUE(bw$bernstein.basis.engine),
    cvls_shadow_basis(bw$basis.engine, bw$regtype.engine),
    compare_old,
    PACKAGE = "np"
  )
}

test_that("public npcdensbw cv.ls keeps lc on the legacy objective", {
  set.seed(222)
  n <- 32L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))

  bw.lc <- npcdensbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ls", nmulti = 0)
  shadow <- call_public_cvls_shadow(bw.lc, x, y)

  expect_equal(bw.lc$fval, -shadow$old, tolerance = 1e-10)
})

test_that("public npcdensbw cv.ls fixed LP/LL route activates with ll == lp parity", {
  set.seed(141)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    nmulti = 0
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 0
  )

  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ls fixed LP route preserves tree parity", {
  set.seed(142)
  n <- 34L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.12))
  degree <- rep.int(1L, ncol(x))

  bw.serial <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 0
  )

  old_opt <- getOption("np.tree")
  on.exit(options(np.tree = old_opt), add = TRUE)
  options(np.tree = TRUE)

  bw.tree <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 0
  )

  expect_equal(bw.tree$fval, bw.serial$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ls generalized-nn LP route activates with ll == lp parity", {
  set.seed(143)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 0
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 0
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ls adaptive-nn LP route activates with ll == lp parity", {
  set.seed(144)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 0
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 0
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})
