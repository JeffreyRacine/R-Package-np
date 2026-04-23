library(npRmpi)

chisq_support_fixture <- function(n, seed) {
  set.seed(seed)
  x <- runif(n, 0, 1)
  y <- rchisq(n, df = 2 + 4 * (x - 0.5)^2)
  list(x = data.frame(x = x), y = data.frame(y = y))
}

test_that("npcdensbw stores the cv.ls I1 rescue toggle on conditional bandwidth objects", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- chisq_support_fixture(n = 40L, seed = 20260423L)

  bw_off <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed"
  )
  bw_on <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cvls.i1.rescue = TRUE
  )

  expect_false(isTRUE(bw_off$cvls.i1.rescue))
  expect_true(isTRUE(bw_on$cvls.i1.rescue))
})

test_that("cv.ls I1 rescue leaves resolved bounded objectives unchanged", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- chisq_support_fixture(n = 80L, seed = 20260423L)

  bw_off <- npcdensbw(
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
    cvls.i1.rescue = FALSE
  )
  bw_on <- bw_off
  bw_on$cvls.i1.rescue <- TRUE

  obj_off <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_off)$objective
  obj_on <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_on)$objective

  expect_true(is.finite(obj_off))
  expect_true(is.finite(obj_on))
  expect_equal(obj_on, obj_off, tolerance = 1e-12)
})

test_that("cv.ls I1 rescue penalizes the known bad one-sided tiny-hy candidate", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- chisq_support_fixture(n = 400L, seed = 600007L)

  bw_off <- npcdensbw(
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
    cvls.i1.rescue = FALSE
  )
  bw_on <- bw_off
  bw_on$cvls.i1.rescue <- TRUE

  obj_off <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_off)$objective
  obj_on <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_on)$objective

  expect_true(is.finite(obj_off))
  expect_true(is.finite(obj_on))
  expect_gt(obj_off, 1)
  expect_lt(obj_on, obj_off)
})
