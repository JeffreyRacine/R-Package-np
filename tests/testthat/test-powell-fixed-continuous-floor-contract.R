suppressPackageStartupMessages(library(npRmpi))

test_that("Powell hot-start respects fixed continuous lower bound", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  floor_value <- function(y, n, p = 2L, q = 1L) {
    0.1 * npRmpi:::EssDee(y) * n^(-1 / (2 * p + q))
  }

  set.seed(600007)
  n <- 400
  x <- runif(n)
  y <- rchisq(n, df = 2 + 4 * (x - 0.5)^2)
  floor_y <- floor_value(y, n)

  nom <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bwmethod = "cv.ls",
    regtype = "lp",
    bwtype = "fixed",
    search.engine = "nomad",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 10L,
    bernstein.basis = TRUE,
    nmulti = 2,
    nomad = TRUE,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykerbound = "fixed", cykerlb = 0, cykerub = Inf
  )

  expect_gte(nom$ybw[1], floor_y)

  nom$degree <- 0L
  nom$degree.engine <- 0L
  nom$regtype <- "lp"
  nom$pregtype <- "Local-Polynomial"

  hot <- npRmpi:::npcdensbw.conbandwidth(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = nom,
    bandwidth.compute = TRUE,
    nmulti = 1
  )

  expect_gte(hot$ybw[1], floor_y)
})

test_that("non-search explicit bandwidths are not clamped", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(77)
  n <- 80
  x <- runif(n)
  y <- rexp(n)

  base <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bwmethod = "cv.ls",
    regtype = "lp",
    degree = 0L,
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykerbound = "fixed", cykerlb = 0, cykerub = Inf
  )

  base$xbw[1] <- 0.12345
  base$ybw[1] <- 1e-06

  out <- npRmpi:::npcdensbw.conbandwidth(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = base,
    bandwidth.compute = FALSE
  )

  expect_equal(out$xbw[1], 0.12345, tolerance = 0)
  expect_equal(out$ybw[1], 1e-06, tolerance = 0)
})

test_that("public eval-only bounded cv.ls objective still accepts tiny fixed bandwidths", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 350L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rbeta(n, 1, 1))

  bw.bad <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.0610566549, 0.0009802178),
    bwtype = "fixed",
    bwmethod = "cv.ls",
    regtype = "lp",
    degree = 3,
    basis = "glp",
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )

  obj.bad <- npRmpi:::.npcdensbw_eval_only(x, y, bw.bad)$objective
  expect_true(is.finite(obj.bad))
})
