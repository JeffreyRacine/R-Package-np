suppressPackageStartupMessages(library(npRmpi))

chisq_support_fixture <- function(n, seed) {
  set.seed(seed)
  x <- runif(n, 0, 1)
  y <- rchisq(n, df = 2 + 4 * (x - 0.5)^2)
  list(x = data.frame(x = x), y = data.frame(y = y))
}

make_bad_seed_bandwidth <- function(scale.factor.search.lower = NULL) {
  dat <- chisq_support_fixture(n = 400L, seed = 600007L)
  args <- list(
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
    cykerub = Inf
  )
  if (!is.null(scale.factor.search.lower)) {
    args$scale.factor.search.lower <- scale.factor.search.lower
  }
  list(
    data = dat,
    bw = do.call(npcdensbw, args)
  )
}

test_that("omitted scale-factor floor matches explicit default 0.1", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  default_case <- make_bad_seed_bandwidth()
  strict_case <- make_bad_seed_bandwidth(scale.factor.search.lower = 0.1)
  legacy_case <- make_bad_seed_bandwidth(scale.factor.search.lower = 0.01)

  obj_default <- npRmpi:::.npcdensbw_eval_only(default_case$data$x, default_case$data$y, default_case$bw)$objective
  obj_strict <- npRmpi:::.npcdensbw_eval_only(strict_case$data$x, strict_case$data$y, strict_case$bw)$objective
  obj_legacy <- npRmpi:::.npcdensbw_eval_only(legacy_case$data$x, legacy_case$data$y, legacy_case$bw)$objective

  expect_equal(default_case$bw$scale.factor.search.lower, 0.1, tolerance = 0)
  expect_equal(strict_case$bw$scale.factor.search.lower, 0.1, tolerance = 0)
  expect_equal(legacy_case$bw$scale.factor.search.lower, 0.01, tolerance = 0)
  expect_equal(obj_default, obj_strict, tolerance = 1e-12)
  expect_equal(obj_default, obj_legacy, tolerance = 1e-12)
})

test_that("explicit bandwidth storage is unchanged when bandwidth.compute is FALSE", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  strict_case <- make_bad_seed_bandwidth(scale.factor.search.lower = 0.1)

  expect_equal(strict_case$bw$scale.factor.search.lower, 0.1, tolerance = 0)
  expect_equal(strict_case$bw$ybw[1L], 1.94042638343838e-05, tolerance = 0)
  expect_equal(strict_case$bw$xbw[1L], 2455873.66968089, tolerance = 0)
})

test_that("explicit 0.1 floor does not reinterpret eval-only objective values", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  default_case <- make_bad_seed_bandwidth()
  strict_case <- make_bad_seed_bandwidth(scale.factor.search.lower = 0.1)

  obj_default <- npRmpi:::.npcdensbw_eval_only(default_case$data$x, default_case$data$y, default_case$bw)$objective
  obj_strict <- npRmpi:::.npcdensbw_eval_only(strict_case$data$x, strict_case$data$y, strict_case$bw)$objective

  expect_equal(obj_strict, obj_default, tolerance = 1e-12)
})

test_that("explicit 0.1 floor keeps nomad and nomad+powell aligned on bad seed", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- chisq_support_fixture(n = 400L, seed = 600007L)

  common_args <- list(
    xdat = dat$x,
    ydat = dat$y,
    bwmethod = "cv.ls",
    regtype = "lp",
    bwtype = "fixed",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 10L,
    bernstein.basis = TRUE,
    nmulti = 2,
    scale.factor.search.lower = 0.1,
    cvls.quadrature.points = c(81L, 31L),
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykerbound = "fixed", cykerlb = 0, cykerub = Inf
  )

  nomad <- do.call(npcdensbw, c(common_args, list(search.engine = "nomad")))
  hot <- do.call(npcdensbw, c(common_args, list(search.engine = "nomad+powell")))

  expect_equal(hot$degree, nomad$degree, tolerance = 0)
  expect_gte(hot$sfactor$y[1L], hot$scale.factor.search.lower)
  expect_true(is.finite(hot$fval[1L]))
})

test_that("explicit high floor is enforced during conditional-density Powell search", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 40L
  dat <- data.frame(x = rnorm(n), y = rnorm(n))

  out <- npcdensbw(
    y ~ x,
    data = dat,
    bwmethod = "cv.ls",
    scale.factor.search.lower = 1,
    nmulti = 1L,
    itmax = 200L
  )

  expect_true(is.finite(out$fval[1L]))
  expect_gte(out$sfactor$x[1L], out$scale.factor.search.lower)
  expect_gte(out$sfactor$y[1L], out$scale.factor.search.lower)
})
