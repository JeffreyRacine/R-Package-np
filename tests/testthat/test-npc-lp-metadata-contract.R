library(np)

set.seed(101)

make_xy <- function(n = 24L) {
  x <- data.frame(x = stats::runif(n))
  y <- data.frame(y = x$x + stats::rnorm(n, sd = 0.1))
  list(x = x, y = y)
}

make_cat_xy <- function(n = 80L) {
  x <- data.frame(
    f = factor(stats::rbinom(n, 1L, 0.45)),
    g = ordered(sample(letters[1:3], n, replace = TRUE))
  )
  y <- data.frame(
    y = 0.4 * as.integer(x$f) + 0.2 * as.integer(x$g) +
      stats::rnorm(n, sd = 0.3)
  )
  list(x = x, y = y)
}

bw_num <- function(bw) as.numeric(unlist(bw$bw, use.names = FALSE))

test_that("npcdensbw stores canonical ll/lp metadata", {
  d <- make_xy()

  bw.lc <- npcdensbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE
  )
  expect_identical(bw.lc$regtype, "lc")
  expect_identical(bw.lc$regtype.engine, "lc")

  bw.ll <- npcdensbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  expect_identical(bw.ll$regtype, "ll")
  expect_identical(bw.ll$regtype.engine, "lp")
  expect_identical(bw.ll$basis.engine, "glp")
  expect_identical(as.integer(bw.ll$degree.engine), 1L)
  expect_false(isTRUE(bw.ll$bernstein.basis.engine))

  bw.lp <- npcdensbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "tensor",
    degree = 2L,
    bernstein.basis = TRUE
  )
  expect_identical(bw.lp$regtype, "lp")
  expect_identical(bw.lp$regtype.engine, "lp")
  expect_identical(bw.lp$basis.engine, "tensor")
  expect_identical(as.integer(bw.lp$degree.engine), 2L)
  expect_true(isTRUE(bw.lp$bernstein.basis.engine))
})

test_that("npcdistbw stores canonical ll/lp metadata", {
  d <- make_xy()

  bw.ll <- npcdistbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  expect_identical(bw.ll$regtype, "ll")
  expect_identical(bw.ll$regtype.engine, "lp")
  expect_identical(bw.ll$basis.engine, "glp")
  expect_identical(as.integer(bw.ll$degree.engine), 1L)

  bw.lp <- npcdistbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "additive",
    degree = 3L
  )
  expect_identical(bw.lp$regtype.engine, "lp")
  expect_identical(bw.lp$basis.engine, "additive")
  expect_identical(as.integer(bw.lp$degree.engine), 3L)
})

test_that("npc* conditional regtype argument contracts fail fast", {
  d <- make_xy()

  expect_error(
    npcdensbw(
      xdat = d$x,
      ydat = d$y,
      bws = c(0.35, 0.35),
      bandwidth.compute = FALSE,
      regtype = "ll",
      degree = 2L
    ),
    "canonical LP\\(degree=1"
  )

  expect_error(
    npcdistbw(
      xdat = d$x,
      ydat = d$y,
      bws = c(0.35, 0.35),
      bandwidth.compute = FALSE,
      regtype = "lc",
      basis = "glp"
    ),
    "regtype='lc' does not accept basis/degree/bernstein.basis"
  )
})

test_that("npcdens ll matches lp(degree=1, basis='glp')", {
  d <- make_xy()

  bw.ll <- npcdensbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdensbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 1L
  )

  fit.ll <- npcdens(bws = bw.ll, txdat = d$x, tydat = d$y, gradients = TRUE)
  fit.lp <- npcdens(bws = bw.lp, txdat = d$x, tydat = d$y, gradients = TRUE)

  expect_equal(fitted(fit.ll), fitted(fit.lp), tolerance = 1e-10)
  expect_equal(fit.ll$congrad, fit.lp$congrad, tolerance = 1e-10)
})

test_that("npcdist ll matches lp(degree=1, basis='glp')", {
  d <- make_xy()

  bw.ll <- npcdistbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  bw.lp <- npcdistbw(
    xdat = d$x,
    ydat = d$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 1L
  )

  fit.ll <- npcdist(bws = bw.ll, txdat = d$x, tydat = d$y, gradients = TRUE)
  fit.lp <- npcdist(bws = bw.lp, txdat = d$x, tydat = d$y, gradients = TRUE)

  expect_equal(fitted(fit.ll), fitted(fit.lp), tolerance = 1e-10)
  expect_equal(fit.ll$congrad, fit.lp$congrad, tolerance = 1e-10)
})

test_that("npcdens categorical-only predictors reject impossible degree structure", {
  d <- make_cat_xy()

  set.seed(90210)
  bw.lc <- npcdensbw(
    xdat = d$x,
    ydat = d$y,
    regtype = "lc",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  expect_error(
    npcdensbw(xdat = d$x, ydat = d$y, regtype = "ll",
              bwmethod = "cv.ls", nmulti = 1L),
    "requires at least one continuous predictor"
  )
  expect_error(
    npcdensbw(xdat = d$x, ydat = d$y, regtype = "lp", degree = 1L,
              bwmethod = "cv.ls", nmulti = 1L),
    "degree must be 0"
  )
  expect_error(
    npcdensbw(xdat = d$x, ydat = d$y, regtype = "lp",
              bwmethod = "cv.ls", nmulti = 1L),
    "degree must be 0"
  )
  expect_error(
    npcdensbw(xdat = d$x, ydat = d$y, regtype = "lp",
              degree.select = "coordinate", search.engine = "cell",
              bwmethod = "cv.ls", nmulti = 1L),
    "automatic degree search requires at least one continuous"
  )
  set.seed(90210)
  bw.lp <- npcdensbw(
    xdat = d$x,
    ydat = d$y,
    regtype = "lp",
    degree = 0L,
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(bw.lp$regtype, "lp")
  expect_identical(bw.lp$regtype.engine, "lc")
  expect_equal(as.numeric(bw.lc$fval), as.numeric(bw.lp$fval), tolerance = 1e-12)
  expect_equal(bw_num(bw.lc), bw_num(bw.lp), tolerance = 1e-12)
  expect_error(
    npcdensbw(xdat = d$x, ydat = d$y, nomad = TRUE),
    "nomad=TRUE requires at least one continuous predictor"
  )
})

test_that("npcdist categorical-only predictors reject impossible degree structure", {
  d <- make_cat_xy()

  set.seed(90210)
  bw.lc <- npcdistbw(
    xdat = d$x,
    ydat = d$y,
    regtype = "lc",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  expect_error(
    npcdistbw(xdat = d$x, ydat = d$y, regtype = "ll",
              bwmethod = "cv.ls", nmulti = 1L),
    "requires at least one continuous predictor"
  )
  expect_error(
    npcdistbw(xdat = d$x, ydat = d$y, regtype = "lp", degree = 1L,
              bwmethod = "cv.ls", nmulti = 1L),
    "degree must be 0"
  )
  expect_error(
    npcdistbw(xdat = d$x, ydat = d$y, regtype = "lp",
              bwmethod = "cv.ls", nmulti = 1L),
    "degree must be 0"
  )
  expect_error(
    npcdistbw(xdat = d$x, ydat = d$y, regtype = "lp",
              degree.select = "coordinate", search.engine = "cell",
              bwmethod = "cv.ls", nmulti = 1L),
    "automatic degree search requires at least one continuous"
  )
  set.seed(90210)
  bw.lp <- npcdistbw(
    xdat = d$x,
    ydat = d$y,
    regtype = "lp",
    degree = 0L,
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(bw.lp$regtype, "lp")
  expect_identical(bw.lp$regtype.engine, "lc")
  expect_equal(as.numeric(bw.lc$fval), as.numeric(bw.lp$fval), tolerance = 1e-12)
  expect_equal(bw_num(bw.lc), bw_num(bw.lp), tolerance = 1e-12)
  expect_error(
    npcdistbw(xdat = d$x, ydat = d$y, nomad = TRUE),
    "nomad=TRUE requires at least one continuous predictor"
  )
})
