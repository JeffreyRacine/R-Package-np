library(np)

adaptive_shadow_empty <- function(n) {
  matrix(numeric(0), nrow = n, ncol = 0)
}

adaptive_shadow_cker <- function(kernel) {
  switch(kernel,
    gaussian = 0L,
    epanechnikov = 4L,
    uniform = 8L,
    truncated = 9L,
    stop("unsupported continuous kernel")
  )
}

adaptive_shadow_uker <- function(kernel) {
  switch(kernel,
    aitchisonaitken = 0L,
    liracine = 1L,
    stop("unsupported unordered kernel")
  )
}

adaptive_shadow_oker <- function(kernel) {
  switch(kernel,
    wangvanryzin = 0L,
    liracine = 2L,
    racineliyan = 3L,
    stop("unsupported ordered kernel")
  )
}

adaptive_shadow_rbw <- function(bw) {
  c(
    bw$xbw[bw$ixcon],
    bw$ybw[bw$iycon],
    bw$ybw[bw$iyuno],
    bw$ybw[bw$iyord],
    bw$xbw[bw$ixuno],
    bw$xbw[bw$ixord]
  )
}

adaptive_shadow_bwtype <- function(bw) {
  switch(bw$type,
    fixed = 0L,
    generalized_nn = 1L,
    adaptive_nn = 2L,
    stop("unsupported bandwidth type")
  )
}

adaptive_shadow_regtype <- function(bw) {
  if (identical(bw$regtype.engine, "lp")) 2L else 0L
}

adaptive_shadow_degree <- function(bw) {
  if (identical(bw$regtype.engine, "lp")) as.integer(bw$degree.engine) else integer(0)
}

adaptive_shadow_basis <- function(basis_engine, regtype_engine) {
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

adaptive_shadow_safe_call <- function(name, ...) {
  on.exit(
    tryCatch(.Call("C_np_shadow_reset_state", PACKAGE = "np"),
             error = function(e) NULL),
    add = TRUE
  )
  .Call(name, ..., PACKAGE = "np")
}

call_adaptive_density_shadow <- function(bw, x, y, criterion = c("cv.ml", "cv.ls")) {
  criterion <- match.arg(criterion)
  n <- nrow(x)
  adaptive_shadow_safe_call(
    "C_np_shadow_cv_density_conditional",
    adaptive_shadow_empty(n), adaptive_shadow_empty(n), as.matrix(y),
    adaptive_shadow_empty(n), adaptive_shadow_empty(n), as.matrix(x),
    as.double(adaptive_shadow_rbw(bw)),
    adaptive_shadow_bwtype(bw),
    adaptive_shadow_cker(bw$cykertype),
    adaptive_shadow_uker(bw$uykertype),
    adaptive_shadow_oker(bw$oykertype),
    adaptive_shadow_cker(bw$cxkertype),
    adaptive_shadow_uker(bw$uxkertype),
    adaptive_shadow_oker(bw$oxkertype),
    FALSE,
    if (identical(criterion, "cv.ml")) 0L else 1L,
    adaptive_shadow_regtype(bw),
    adaptive_shadow_degree(bw),
    isTRUE(bw$bernstein.basis.engine),
    adaptive_shadow_basis(bw$basis.engine, bw$regtype.engine),
    identical(bw$regtype.engine, "lc")
  )
}

call_adaptive_distribution_shadow <- function(bw, x, ytrain, yeval = ytrain, cdfontrain = FALSE) {
  n <- nrow(x)
  ne <- nrow(yeval)
  adaptive_shadow_safe_call(
    "C_np_shadow_cv_distribution_conditional",
    adaptive_shadow_empty(n), adaptive_shadow_empty(n), as.matrix(ytrain),
    adaptive_shadow_empty(ne), adaptive_shadow_empty(ne), as.matrix(yeval),
    adaptive_shadow_empty(n), adaptive_shadow_empty(n), as.matrix(x),
    as.double(adaptive_shadow_rbw(bw)),
    adaptive_shadow_bwtype(bw),
    adaptive_shadow_cker(bw$cykertype),
    adaptive_shadow_uker(bw$uykertype),
    adaptive_shadow_oker(bw$oykertype),
    adaptive_shadow_cker(bw$cxkertype),
    adaptive_shadow_uker(bw$uxkertype),
    adaptive_shadow_oker(bw$oxkertype),
    FALSE,
    adaptive_shadow_regtype(bw),
    adaptive_shadow_degree(bw),
    isTRUE(bw$bernstein.basis.engine),
    adaptive_shadow_basis(bw$basis.engine, bw$regtype.engine),
    cdfontrain,
    identical(bw$regtype.engine, "lc")
  )
}

adaptive_cdist_eval_grid <- function(ydat, ngrid = 48L) {
  probs <- seq(0, 1, length.out = ngrid)
  evy <- ydat[rep(1L, ngrid), , drop = FALSE]
  for (i in seq_len(ncol(evy))) {
    evy[, i] <- uocquantile(ydat[, i], probs)
  }
  evy
}

adaptive_cv_fixture <- function(seed) {
  set.seed(seed)
  n <- 36L
  x <- data.frame(
    x1 = runif(n),
    x2 = runif(n)
  )
  y <- data.frame(
    y1 = sin(2 * pi * x$x1) + 0.5 * x$x2^2 + rnorm(n, sd = 0.08)
  )
  list(
    x = x,
    y = y,
    degree1 = rep.int(1L, ncol(x)),
    degree2 = rep.int(2L, ncol(x)),
    bws = c(5, 7, 6)
  )
}

make_adaptive_npcdens_bw <- function(x, y, regtype, bws, degree = NULL, bandwidth.compute = FALSE) {
  args <- list(
    xdat = x,
    ydat = y,
    bwtype = "adaptive_nn",
    bws = bws,
    bandwidth.compute = bandwidth.compute,
    regtype = regtype
  )
  if (identical(regtype, "lp")) {
    args$basis <- "glp"
    args$degree <- degree
  }
  do.call(npcdensbw, args)
}

make_adaptive_npcdist_bw <- function(x, y, regtype, bws, degree = NULL, bandwidth.compute = FALSE) {
  args <- list(
    xdat = x,
    ydat = y,
    bwtype = "adaptive_nn",
    bws = bws,
    bandwidth.compute = bandwidth.compute,
    regtype = regtype,
    itmax = 1L
  )
  if (identical(regtype, "lp")) {
    args$basis <- "glp"
    args$degree <- degree
  }
  do.call(npcdistbw, args)
}

test_that("adaptive conditional density shadow CV LS activates LP production route", {
  fixture <- adaptive_cv_fixture(4101)
  bw.lc <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lc", fixture$bws)
  bw.ll <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "ll", fixture$bws)
  bw.lp1 <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree1)
  bw.lp2 <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree2)

  shadow.lc <- call_adaptive_density_shadow(bw.lc, fixture$x, fixture$y, criterion = "cv.ls")
  shadow.ll <- call_adaptive_density_shadow(bw.ll, fixture$x, fixture$y, criterion = "cv.ls")
  shadow.lp1 <- call_adaptive_density_shadow(bw.lp1, fixture$x, fixture$y, criterion = "cv.ls")
  shadow.lp2 <- call_adaptive_density_shadow(bw.lp2, fixture$x, fixture$y, criterion = "cv.ls")

  expect_true(is.finite(shadow.ll$new))
  expect_true(is.finite(shadow.lp1$new))
  expect_true(is.finite(shadow.lp2$new))
  expect_true(is.finite(shadow.lp1$prod))
  expect_equal(shadow.ll$new, shadow.lp1$new, tolerance = 1e-10)
  expect_equal(shadow.ll$prod, shadow.lp1$prod, tolerance = 1e-10)
  expect_gt(abs(shadow.lc$new - shadow.lp1$new), 1e-6)
  expect_gt(abs(shadow.lp2$new - shadow.lp1$new), 1e-6)
})

test_that("adaptive conditional density shadow CV ML activates LP production route", {
  fixture <- adaptive_cv_fixture(4102)
  bw.lc <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lc", fixture$bws)
  bw.ll <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "ll", fixture$bws)
  bw.lp1 <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree1)
  bw.lp2 <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree2)

  shadow.lc <- call_adaptive_density_shadow(bw.lc, fixture$x, fixture$y, criterion = "cv.ml")
  shadow.ll <- call_adaptive_density_shadow(bw.ll, fixture$x, fixture$y, criterion = "cv.ml")
  shadow.lp1 <- call_adaptive_density_shadow(bw.lp1, fixture$x, fixture$y, criterion = "cv.ml")
  shadow.lp2 <- call_adaptive_density_shadow(bw.lp2, fixture$x, fixture$y, criterion = "cv.ml")

  expect_true(is.finite(shadow.ll$new))
  expect_true(is.finite(shadow.lp1$new))
  expect_true(is.finite(shadow.lp2$new))
  expect_true(is.finite(shadow.lp1$prod))
  expect_equal(shadow.ll$new, shadow.lp1$new, tolerance = 1e-10)
  expect_equal(shadow.ll$prod, shadow.lp1$prod, tolerance = 1e-10)
  expect_gt(abs(shadow.lc$new - shadow.lp1$new), 1e-6)
  expect_gt(abs(shadow.lp2$new - shadow.lp1$new), 1e-6)
})

test_that("adaptive conditional distribution shadow CV LS activates LP production route", {
  fixture <- adaptive_cv_fixture(4103)
  yeval <- adaptive_cdist_eval_grid(fixture$y)
  bw.lc <- make_adaptive_npcdist_bw(fixture$x, fixture$y, "lc", fixture$bws)
  bw.ll <- make_adaptive_npcdist_bw(fixture$x, fixture$y, "ll", fixture$bws)
  bw.lp1 <- make_adaptive_npcdist_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree1)
  bw.lp2 <- make_adaptive_npcdist_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree2)

  shadow.lc <- call_adaptive_distribution_shadow(bw.lc, fixture$x, fixture$y, yeval = yeval)
  shadow.ll <- call_adaptive_distribution_shadow(bw.ll, fixture$x, fixture$y, yeval = yeval)
  shadow.lp1 <- call_adaptive_distribution_shadow(bw.lp1, fixture$x, fixture$y, yeval = yeval)
  shadow.lp2 <- call_adaptive_distribution_shadow(bw.lp2, fixture$x, fixture$y, yeval = yeval)

  expect_true(is.finite(shadow.ll$new))
  expect_true(is.finite(shadow.lp1$new))
  expect_true(is.finite(shadow.lp2$new))
  expect_true(is.finite(shadow.lp1$prod))
  expect_equal(shadow.ll$new, shadow.lp1$new, tolerance = 1e-10)
  expect_equal(shadow.lp1$new, shadow.lp1$prod, tolerance = 1e-10)
  expect_gt(abs(shadow.lc$new - shadow.lp1$new), 1e-6)
  expect_gt(abs(shadow.lp2$new - shadow.lp1$new), 1e-6)
})

test_that("adaptive public conditional density CV LS separates lc from LP while preserving ll canonicalization", {
  fixture <- adaptive_cv_fixture(4201)
  bw.lc <- npcdensbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lc",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1
  )
  bw.ll <- npcdensbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1
  )
  bw.lp1 <- npcdensbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lp",
    basis = "glp",
    degree = fixture$degree1,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1
  )
  bw.lp2 <- npcdensbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lp",
    basis = "glp",
    degree = fixture$degree2,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1
  )

  expect_equal(bw.ll$fval, bw.lp1$fval, tolerance = 1e-8)
  expect_gt(abs(bw.lc$fval - bw.lp1$fval), 1e-6)
  expect_gt(abs(bw.lp2$fval - bw.lp1$fval), 1e-6)
})

test_that("adaptive public conditional density CV ML separates lc from LP while preserving ll canonicalization", {
  fixture <- adaptive_cv_fixture(4202)
  bw.lc <- npcdensbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lc",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1
  )
  bw.ll <- npcdensbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1
  )
  bw.lp1 <- npcdensbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lp",
    basis = "glp",
    degree = fixture$degree1,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1
  )
  bw.lp2 <- npcdensbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lp",
    basis = "glp",
    degree = fixture$degree2,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ml",
    nmulti = 1
  )

  expect_equal(bw.ll$fval, bw.lp1$fval, tolerance = 1e-8)
  expect_gt(abs(bw.lc$fval - bw.lp1$fval), 1e-6)
  expect_gt(abs(bw.lp2$fval - bw.lp1$fval), 1e-6)
})

test_that("adaptive public conditional distribution CV LS separates lc from LP while preserving ll canonicalization", {
  fixture <- adaptive_cv_fixture(4203)
  bw.lc <- npcdistbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lc",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.ll <- npcdistbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.lp1 <- npcdistbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lp",
    basis = "glp",
    degree = fixture$degree1,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )
  bw.lp2 <- npcdistbw(
    xdat = fixture$x,
    ydat = fixture$y,
    regtype = "lp",
    basis = "glp",
    degree = fixture$degree2,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1L
  )

  expect_equal(bw.ll$fval, bw.lp1$fval, tolerance = 1e-8)
  expect_gt(abs(bw.lc$fval - bw.lp1$fval), 1e-6)
  expect_gt(abs(bw.lp2$fval - bw.lp1$fval), 1e-6)
})
