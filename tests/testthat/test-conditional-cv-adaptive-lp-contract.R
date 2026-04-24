library(npRmpi)

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
    degree2 = rep.int(2L, ncol(x))
  )
}

test_that("adaptive conditional CV LP source contract is widened in jksum", {
  src <- testthat::test_path("..", "..", "src", "jksum.c")
  skip_if_not(file.exists(src))

  lines <- readLines(src, warn = FALSE)

  expect_true(any(grepl("BW_GEN_NN\\) \\|\\| \\(BANDWIDTH_den == BW_ADAP_NN\\)", lines)))
  expect_true(any(grepl("BANDWIDTH_den_extern == BW_ADAP_NN", lines)))
})

test_that("adaptive conditional distribution CVLS uses local compiled routing contract in npRmpi", {
  src <- testthat::test_path("..", "..", "R", "np.condistribution.bw.R")
  skip_if_not(file.exists(src))

  lines <- readLines(src, warn = FALSE)

  expect_true(any(grepl("\\.npRmpi_with_local_cdist_eval", lines)))
  expect_true(any(grepl("identical\\(bws\\$type, \\\"adaptive_nn\\\"\\)", lines)))
})

test_that("adaptive public conditional density CV LS separates lc from LP while preserving ll canonicalization", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- adaptive_cv_fixture(5201)

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- adaptive_cv_fixture(5202)

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- adaptive_cv_fixture(5203)

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
