library(np)

adaptive_density_objective <- function(bw, x, y) {
  np:::.npcdensbw_eval_only(x, y, bw)$objective
}

adaptive_distribution_objective <- function(bw, x, y) {
  np:::.npcdistbw_eval_only(x, y, bws = bw)$objective
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
    # The degree-(2,2) generalized basis has nine terms. Keep every adaptive
    # neighborhood comfortably identified; the retired QR path could return a
    # finite value for the old 5--7-neighbor, underdetermined fixture.
    bws = c(18, 20, 19)
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

test_that("adaptive conditional density CV LS activates the canonical LP route", {
  fixture <- adaptive_cv_fixture(4101)
  bw.lc <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lc", fixture$bws)
  bw.ll <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "ll", fixture$bws)
  bw.lp1 <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree1)
  bw.lp2 <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree2)

  objective <- vapply(
    list(bw.lc, bw.ll, bw.lp1, bw.lp2),
    adaptive_density_objective,
    numeric(1L),
    x = fixture$x,
    y = fixture$y
  )

  expect_true(all(is.finite(objective)))
  expect_equal(objective[[2L]], objective[[3L]], tolerance = 1e-10)
  expect_gt(abs(objective[[1L]] - objective[[3L]]), 1e-6)
  expect_gt(abs(objective[[4L]] - objective[[3L]]), 1e-6)
})

test_that("adaptive conditional density CV ML activates the canonical LP route", {
  fixture <- adaptive_cv_fixture(4102)
  bw.lc <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lc", fixture$bws)
  bw.ll <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "ll", fixture$bws)
  bw.lp1 <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree1)
  bw.lp2 <- make_adaptive_npcdens_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree2)

  objective <- vapply(
    list(bw.lc, bw.ll, bw.lp1, bw.lp2),
    adaptive_density_objective,
    numeric(1L),
    x = fixture$x,
    y = fixture$y
  )

  expect_true(all(is.finite(objective)))
  expect_equal(objective[[2L]], objective[[3L]], tolerance = 1e-10)
  expect_gt(abs(objective[[1L]] - objective[[3L]]), 1e-6)
  expect_gt(abs(objective[[4L]] - objective[[3L]]), 1e-6)
})

test_that("adaptive conditional distribution CV LS activates the canonical LP route", {
  fixture <- adaptive_cv_fixture(4103)
  bw.lc <- make_adaptive_npcdist_bw(fixture$x, fixture$y, "lc", fixture$bws)
  bw.ll <- make_adaptive_npcdist_bw(fixture$x, fixture$y, "ll", fixture$bws)
  bw.lp1 <- make_adaptive_npcdist_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree1)
  bw.lp2 <- make_adaptive_npcdist_bw(fixture$x, fixture$y, "lp", fixture$bws, degree = fixture$degree2)

  objective <- vapply(
    list(bw.lc, bw.ll, bw.lp1, bw.lp2),
    adaptive_distribution_objective,
    numeric(1L),
    x = fixture$x,
    y = fixture$y
  )

  expect_true(all(is.finite(objective)))
  expect_equal(objective[[2L]], objective[[3L]], tolerance = 1e-10)
  expect_gt(abs(objective[[1L]] - objective[[3L]]), 1e-6)
  expect_gt(abs(objective[[4L]] - objective[[3L]]), 1e-6)
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

test_that("adaptive conditional density CV ML degree zero is canonical LC", {
  fixture <- adaptive_cv_fixture(4104)
  eval_only <- getFromNamespace(".npcdensbw_eval_only", "np")
  make_bw <- function(regtype, bernstein = FALSE) {
    args <- list(
      xdat = fixture$x,
      ydat = fixture$y,
      regtype = regtype,
      bwtype = "adaptive_nn",
      bwmethod = "cv.ml",
      bws = fixture$bws,
      bandwidth.compute = FALSE
    )
    if (identical(regtype, "lp")) {
      args$basis <- "glp"
      args$degree <- c(0L, 0L)
      args$bernstein.basis <- bernstein
    }
    do.call(npcdensbw, args)
  }
  objective <- function(bw) {
    eval_only(fixture$x, fixture$y, bw)$objective
  }

  lc <- objective(make_bw("lc"))
  raw <- objective(make_bw("lp", FALSE))
  bernstein <- objective(make_bw("lp", TRUE))

  expect_identical(raw, lc)
  expect_identical(bernstein, lc)
})
