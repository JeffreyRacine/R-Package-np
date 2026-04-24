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

cvls_shadow_safe_call <- function(name, ...) {
  on.exit(
    tryCatch(.Call("C_np_shadow_reset_state", PACKAGE = "np"),
             error = function(e) NULL),
    add = TRUE
  )
  .Call(name, ..., PACKAGE = "np")
}

call_public_cvls_shadow <- function(bw, x, y, compare_old = identical(bw$regtype.engine, "lc")) {
  n <- nrow(x)
  cvls_shadow_safe_call(
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
    compare_old
  )
}

test_that("public npcdensbw cv.ls lc matches the production fixed-point objective", {
  set.seed(222)
  n <- 32L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))

  bw.lc <- npcdensbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  shadow <- call_public_cvls_shadow(bw.lc, x, y)

  expect_equal(bw.lc$fval, -shadow$prod, tolerance = 1e-10)
})

test_that("provided fixed lc cv.ls eval_only remains finite", {
  set.seed(42)
  x <- data.frame(x = runif(80L))
  y <- data.frame(y = rbeta(80L, 1, 1))

  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bwtype = "fixed",
    bwmethod = "cv.ls",
    bws = c(0.4, 0.4),
    bandwidth.compute = FALSE
  )

  out <- np:::.npcdensbw_eval_only(x, y, bw)

  expect_true(is.list(out))
  expect_true(is.finite(out$objective))
  expect_equal(out$num.feval, 1)
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
    nmulti = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 1
  )

  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ls fixed LP tree and serial evaluators agree at fixed points", {
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
    nmulti = 1
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
    nmulti = 1
  )

  options(np.tree = FALSE)
  serial.at.serial <- np:::.npcdensbw_eval_only(x, y$y1, bw.serial)$objective
  serial.at.tree <- np:::.npcdensbw_eval_only(x, y$y1, bw.tree)$objective

  options(np.tree = TRUE)
  tree.at.serial <- np:::.npcdensbw_eval_only(x, y$y1, bw.serial)$objective
  tree.at.tree <- np:::.npcdensbw_eval_only(x, y$y1, bw.tree)$objective

  expect_equal(tree.at.serial, serial.at.serial, tolerance = 2e-2)
  expect_equal(tree.at.tree, serial.at.tree, tolerance = 2e-2)
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
    nmulti = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1
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
    nmulti = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})
