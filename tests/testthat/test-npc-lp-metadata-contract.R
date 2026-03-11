library(npRmpi)

set.seed(101)

make_xy <- function(n = 24L) {
  x <- data.frame(x = stats::runif(n))
  y <- data.frame(y = x$x + stats::rnorm(n, sd = 0.1))
  list(x = x, y = y)
}

test_that("npcdensbw stores canonical ll/lp metadata", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

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
