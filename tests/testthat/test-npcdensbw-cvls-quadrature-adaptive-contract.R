library(npRmpi)

chisq_support_fixture <- function(n, seed) {
  set.seed(seed)
  x <- runif(n, 0, 1)
  y <- rchisq(n, df = 2 + 4 * (x - 0.5)^2)
  list(x = data.frame(x = x), y = data.frame(y = y))
}

test_that("npcdensbw stores the cv.ls adaptive quadrature toggle on conditional bandwidth objects", {
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
  bw_on <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cvls.quadrature.adaptive = TRUE
  )
  bw_off <- npcdensbw(
    xdat = dat$x,
    ydat = dat$y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cvls.quadrature.adaptive = FALSE
  )

  expect_false("cvls.i1.rescue" %in% names(formals(getS3method("npcdensbw", "default"))))
  expect_false("cvls.i1.rescue" %in% names(bw_default))
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
  expect_true(isTRUE(bw_default$cvls.quadrature.adaptive))
  expect_true(isTRUE(bw_on$cvls.quadrature.adaptive))
  expect_false(isTRUE(bw_off$cvls.quadrature.adaptive))
})

test_that("cv.ls adaptive quadrature toggle disables all adaptive triggers", {
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
    cvls.quadrature.adaptive = FALSE
  )
  bw_off_triggered <- bw_off
  bw_off_triggered$cvls.quadrature.adaptive.tol <- 1e6
  bw_off_triggered$cvls.quadrature.adaptive.grid.hy.ratio <- 0
  bw_off_triggered$cvls.quadrature.adaptive.floor.tol <- 1e6

  obj_off <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_off)$objective
  obj_off_triggered <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_off_triggered)$objective

  expect_true(is.finite(obj_off))
  expect_true(is.finite(obj_off_triggered))
  expect_equal(obj_off_triggered, obj_off, tolerance = 1e-12)
})

test_that("cv.ls adaptive quadrature penalizes the known bad one-sided tiny-hy candidate", {
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
    cvls.quadrature.adaptive = FALSE,
    cvls.quadrature.extend.factor = 2,
    cvls.quadrature.points = c(81L, 31L)
  )
  bw_on <- bw_off
  bw_on$cvls.quadrature.adaptive <- TRUE

  obj_off <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_off)$objective
  obj_on <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_on)$objective

  expect_true(is.finite(obj_off))
  expect_true(is.finite(obj_on))
  expect_gt(obj_off, 1)
  expect_lt(obj_on, obj_off)
})
