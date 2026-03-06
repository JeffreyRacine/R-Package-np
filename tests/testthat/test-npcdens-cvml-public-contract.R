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

test_that("public npcdensbw cv.ml LP/LL routes fail fast during containment", {
  set.seed(101)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 - x$x2 + rnorm(n, sd = 0.2))
  degree <- rep.int(1L, ncol(x))

  expect_error(
    npcdensbw(xdat = x, ydat = y, regtype = "ll", bwmethod = "cv.ml", nmulti = 0),
    "temporarily disabled pending low-memory shadow CV remediation"
  )
  expect_error(
    npcdensbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      basis = "glp",
      degree = degree,
      bwmethod = "cv.ml",
      nmulti = 0
    ),
    "temporarily disabled pending low-memory shadow CV remediation"
  )
  expect_error(
    npcdensbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      basis = "glp",
      degree = degree,
      bwtype = "generalized_nn",
      bwmethod = "cv.ml",
      nmulti = 0
    ),
    "temporarily disabled pending low-memory shadow CV remediation"
  )
})
