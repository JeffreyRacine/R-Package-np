library(npRmpi)

test_that("unconditional bounded cv.ls quadrature uses source-of-truth fixed grid", {
  src <- testthat::test_path("..", "..", "src", "jksum.c")
  skip_if_not(file.exists(src), "source file unavailable")

  lines <- readLines(src, warn = FALSE)
  start <- grep("np_density_cvls_bounded_i1_quadrature_general", lines, fixed = TRUE)[1L]
  skip_if(is.na(start), "bounded density quadrature helper unavailable")
  end <- start + grep("const int block_size", lines[start:length(lines)], fixed = TRUE)[1L] - 1L
  helper_header <- paste(lines[start:end], collapse = "\n")

  expect_match(helper_header, "const int q = np_bounded_cvls_grid_points\\(ncon\\);")
  expect_false(grepl("np_bounded_cvls_conditional_grid_points", helper_header, fixed = TRUE))
})

chisq_support_fixture <- function(n, seed) {
  set.seed(seed)
  x <- runif(n, 0, 1)
  y <- rchisq(n, df = 2 + 4 * (x - 0.5)^2)
  list(x = data.frame(x = x), y = data.frame(y = y))
}

test_that("npcdensbw stores the cv.ls quadrature grid mode on conditional bandwidth objects", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- chisq_support_fixture(n = 40L, seed = 20260423L)

  bw_default <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed"
  )
  bw_hybrid <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cvls.quadrature.grid = "hybrid"
  )
  bw_uniform <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cvls.quadrature.grid = "uniform"
  )

  expect_false("cvls.i1.rescue" %in% names(formals(getS3method("npcdensbw", "default"))))
  expect_false("cvls.quadrature.adaptive" %in% names(formals(getS3method("npcdensbw", "default"))))
  expect_false("cvls.i1.rescue" %in% names(bw_default))
  expect_false("cvls.quadrature.adaptive" %in% names(bw_default))
  expect_equal(
    eval(formals(getS3method("npcdensbw", "default"))$cvls.quadrature.ratios),
    c(0.20, 0.55, 0.25)
  )
  expect_error(
    npcdensbw(
      xdat = dat$x,
      ydat = dat$y,
      bws = c(0.35, 0.35),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      cvls.i1.rescue = FALSE
    ),
    "cvls.i1.rescue has been removed"
  )
  expect_identical(bw_default$cvls.quadrature.grid, "hybrid")
  expect_equal(bw_default$cvls.quadrature.ratios, c(0.20, 0.55, 0.25))
  expect_identical(bw_hybrid$cvls.quadrature.grid, "hybrid")
  expect_identical(bw_uniform$cvls.quadrature.grid, "uniform")
})

test_that("cv.ls hybrid quadrature ratios validate and persist", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- chisq_support_fixture(n = 40L, seed = 20260426L)

  bw <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cvls.quadrature.ratios = c(0.2, 0.5, 0.3)
  )

  expect_equal(bw$cvls.quadrature.ratios, c(0.2, 0.5, 0.3))
  expect_error(
    npcdensbw(
      xdat = dat$x,
      ydat = dat$y,
      bws = c(0.35, 0.35),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      cvls.quadrature.ratios = c(0.5, 0.5)
    ),
    "three-element"
  )
  expect_error(
    npcdensbw(
      xdat = dat$x,
      ydat = dat$y,
      bws = c(0.35, 0.35),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      cvls.quadrature.ratios = c(0.5, -0.1, 0.6)
    ),
    "non-negative"
  )
  expect_error(
    npcdensbw(
      xdat = dat$x,
      ydat = dat$y,
      bws = c(0.35, 0.35),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      cvls.quadrature.ratios = c(0.5, 0.25, 0.2)
    ),
    "summing to one"
  )
})

test_that("cv.ls quadrature grid modes are stable finite objective controls", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- chisq_support_fixture(n = 80L, seed = 20260423L)

  bw_uniform <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lp",
    degree = 0,
    cxkerbound = "fixed",
    cxkerlb = 0,
    cxkerub = 1,
    cykerbound = "fixed",
    cykerlb = 0,
    cykerub = Inf,
    cvls.quadrature.grid = "uniform",
    cvls.quadrature.points = c(81L, 31L)
  )
  bw_hybrid <- bw_uniform
  bw_hybrid$cvls.quadrature.grid <- "hybrid"
  bw_sample <- bw_uniform
  bw_sample$cvls.quadrature.grid <- "sample"

  obj_uniform <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_uniform)$objective
  obj_hybrid <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_hybrid)$objective
  obj_sample <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_sample)$objective

  expect_true(is.finite(obj_uniform))
  expect_true(is.finite(obj_hybrid))
  expect_true(is.finite(obj_sample))
  expect_gt(abs(obj_uniform - obj_hybrid), 1e-10)
  expect_gt(abs(obj_uniform - obj_sample), 1e-10)
})

test_that("hybrid cv.ls grid is honored for scalar continuous plus discrete responses", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260425)
  n <- 42L
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(
    y = 0.25 + rchisq(n, df = 3),
    z = factor(sample(letters[1:3], n, replace = TRUE))
  )

  bw_uniform <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.36, 0.24, 0.31),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lp",
    degree = 0,
    cxkerbound = "fixed",
    cxkerlb = 0,
    cxkerub = 1,
    cykerbound = "fixed",
    cykerlb = 0,
    cykerub = Inf,
    cvls.quadrature.grid = "uniform",
    cvls.quadrature.points = c(41L, 17L)
  )
  bw_hybrid <- bw_uniform
  bw_hybrid$cvls.quadrature.grid <- "hybrid"
  bw_sample <- bw_uniform
  bw_sample$cvls.quadrature.grid <- "sample"

  obj_uniform <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw_uniform)$objective
  obj_hybrid <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw_hybrid)$objective
  obj_sample <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw_sample)$objective

  expect_true(is.finite(obj_uniform))
  expect_true(is.finite(obj_hybrid))
  expect_true(is.finite(obj_sample))
  expect_gt(abs(obj_uniform - obj_hybrid), 1e-10)
  expect_gt(abs(obj_uniform - obj_sample), 1e-10)
})

test_that("hybrid cv.ls grid improves the known bad one-sided tiny-hy candidate", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- chisq_support_fixture(n = 400L, seed = 600007L)

  bw_uniform <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(1.94042638343838e-05, 2455873.66968089),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lp",
    degree = 0,
    cxkerbound = "fixed",
    cxkerlb = 0,
    cxkerub = 1,
    cykerbound = "fixed",
    cykerlb = 0,
    cykerub = Inf,
    cvls.quadrature.grid = "uniform",
    cvls.quadrature.extend.factor = 2,
    cvls.quadrature.points = c(81L, 31L)
  )
  bw_hybrid <- bw_uniform
  bw_hybrid$cvls.quadrature.grid <- "hybrid"

  obj_uniform <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_uniform)$objective
  obj_hybrid <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_hybrid)$objective

  expect_true(is.finite(obj_uniform))
  expect_true(is.finite(obj_hybrid))
  expect_gt(obj_uniform, 1)
  expect_lt(obj_hybrid, obj_uniform)
})
