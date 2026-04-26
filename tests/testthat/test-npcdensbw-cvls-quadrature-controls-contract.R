library(npRmpi)

quadrature_control_fixture <- function(n = 48L, seed = 20260424L) {
  set.seed(seed)
  x <- data.frame(x = runif(n))
  y <- data.frame(y = 0.25 + x$x + rchisq(n, df = 3))
  list(x = x, y = y)
}

quadrature_control_bw <- function(dat,
                                  cykerlb = 0,
                                  cykerub = Inf,
                                  ...) {
  npcdensbw(
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
    cykerlb = cykerlb,
    cykerub = cykerub,
    ...
  )
}

test_that("npcdensbw validates cv.ls quadrature controls", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- quadrature_control_fixture(n = 12L)

  expect_no_error(quadrature_control_bw(dat, cvls.quadrature.extend.factor = 0.5))
  expect_no_error(quadrature_control_bw(dat, cvls.quadrature.extend.factor = 1))
  expect_no_error(quadrature_control_bw(dat, cvls.quadrature.extend.factor = 2))
  expect_no_error(quadrature_control_bw(dat, cvls.quadrature.points = c(41L, 17L)))
  expect_no_error(quadrature_control_bw(dat, cvls.quadrature.grid = "hybrid"))
  expect_no_error(quadrature_control_bw(dat, cvls.quadrature.grid = "uniform"))
  expect_no_error(quadrature_control_bw(dat, cvls.quadrature.grid = "sample"))

  removed_args <- list(
    cvls.i1.rescue = FALSE,
    cvls.quadrature.adaptive = TRUE,
    cvls.quadrature.adaptive.tol = 0,
    cvls.quadrature.adaptive.grid.hy.ratio = 0,
    cvls.quadrature.adaptive.floor.tol = 0
  )
  for (nm in names(removed_args)) {
    args <- list(dat)
    args[[nm]] <- removed_args[[nm]]
    expect_error(
      do.call(quadrature_control_bw, args),
      sprintf("%s has been removed", nm),
      fixed = TRUE
    )
  }

  bad_extend <- list(0, -1, NA_real_, NaN, Inf, "2", c(1, 2))
  for (value in bad_extend) {
    expect_error(
      quadrature_control_bw(dat, cvls.quadrature.extend.factor = value),
      "cvls.quadrature.extend.factor"
    )
  }

  bad_points <- list(1, c(81), c(81, 1), c(81, NA), c(81, Inf), c(81.5, 31), "81", c(81, 31, 21))
  for (value in bad_points) {
    expect_error(
      quadrature_control_bw(dat, cvls.quadrature.points = value),
      "cvls.quadrature.points"
    )
  }

  bad_grid <- list(NA, c("hybrid", "uniform"), TRUE, "adaptive")
  for (value in bad_grid) {
    expect_error(
      quadrature_control_bw(dat, cvls.quadrature.grid = value),
      "cvls.quadrature.grid"
    )
  }
})

test_that("npcdensbw stores cv.ls quadrature controls and old objects use defaults", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- quadrature_control_fixture(n = 16L)

  bw_default <- quadrature_control_bw(dat)
  bw_explicit <- quadrature_control_bw(
    dat,
    cvls.quadrature.grid = "sample",
    cvls.quadrature.extend.factor = 1.5,
    cvls.quadrature.points = c(43L, 19L)
  )

  expect_identical(bw_default$cvls.quadrature.grid, "hybrid")
  expect_equal(bw_default$cvls.quadrature.extend.factor, 1)
  expect_identical(unname(bw_default$cvls.quadrature.points), c(100L, 50L))
  expect_identical(bw_explicit$cvls.quadrature.grid, "sample")
  expect_equal(bw_explicit$cvls.quadrature.extend.factor, 1.5)
  expect_identical(unname(bw_explicit$cvls.quadrature.points), c(43L, 19L))

  bw_old <- bw_default
  bw_old$cvls.quadrature.grid <- NULL
  bw_old$cvls.quadrature.extend.factor <- NULL
  bw_old$cvls.quadrature.points <- NULL
  expect_true(is.finite(npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_old)$objective))
})

test_that("resident npcdens NOMAD shadow accepts cv.ls quadrature controls", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- quadrature_control_fixture(n = 18L)
  bw <- quadrature_control_bw(
    dat,
    cvls.quadrature.grid = "uniform",
    cvls.quadrature.points = c(41L, 17L),
    cvls.quadrature.extend.factor = 1.5
  )
  prep <- npRmpi:::.npcdensbw_nomad_shadow_prepare_args(
    xdat = dat$x,
    ydat = dat$y,
    bws = bw,
    invalid.penalty = "baseline"
  )

  expect_length(prep$myopti, 29L)
  expect_length(prep$myoptd, 21L)
  expect_equal(prep$myopti[[28L]], 0L)
  expect_equal(prep$myopti[[29L]], 41L)
  expect_equal(prep$myoptd[[21L]], 1.5)

  prepare_body <- paste(
    deparse(body(npRmpi:::npRmpiNomadShadowPrepareConditionalDensity)),
    collapse = "\n"
  )
  expect_match(prepare_body, "length\\(myoptd\\) <= 20L")
  expect_false(grepl("length(myoptd) <= 23L", prepare_body, fixed = TRUE))
})

test_that("explicit infinite response bounds warn when quadrature points are implicit", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- quadrature_control_fixture(n = 12L)

  expect_warning(
    npRmpi:::.npcdensbw_warn_infinite_response_quadrature(
      kerlb = 0,
      kerub = Inf,
      kerbound = "fixed",
      points.supplied = FALSE,
      where = "npcdensbw()"
    ),
    "fixed infinite response bounds",
    fixed = TRUE
  )
  expect_silent(
    npRmpi:::.npcdensbw_warn_infinite_response_quadrature(
      kerlb = 0,
      kerub = Inf,
      kerbound = "fixed",
      points.supplied = TRUE,
      where = "npcdensbw()"
    )
  )
  expect_silent(
    npRmpi:::.npcdensbw_warn_infinite_response_quadrature(
      kerlb = 0,
      kerub = max(dat$y$y),
      kerbound = "fixed",
      points.supplied = FALSE,
      where = "npcdensbw()"
    )
  )
})

test_that("finite response bounds are invariant to cv.ls quadrature extend factor", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- quadrature_control_fixture(n = 44L)
  finite_ub <- max(dat$y$y) + 0.5

  bw_factor1 <- quadrature_control_bw(
    dat,
    cykerub = finite_ub,
    cvls.quadrature.extend.factor = 1
  )
  bw_factor2 <- bw_factor1
  bw_factor2$cvls.quadrature.extend.factor <- 2

  obj_factor1 <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_factor1)$objective
  obj_factor2 <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_factor2)$objective

  expect_true(is.finite(obj_factor1))
  expect_equal(obj_factor2, obj_factor1, tolerance = 1e-12)
})

test_that("cv.ls quadrature point vector controls one- and two-dimensional grids", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- quadrature_control_fixture(n = 40L)

  bw_1d_default <- quadrature_control_bw(dat, cvls.quadrature.points = c(100L, 50L))
  bw_1d_coarse <- quadrature_control_bw(dat, cvls.quadrature.points = c(41L, 31L))
  obj_1d_default <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_1d_default)$objective
  obj_1d_coarse <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_1d_coarse)$objective

  expect_true(is.finite(obj_1d_default))
  expect_true(is.finite(obj_1d_coarse))
  expect_gt(abs(obj_1d_default - obj_1d_coarse), 1e-10)

  set.seed(20260424)
  n2 <- 16L
  x2 <- data.frame(x = runif(n2))
  y2 <- data.frame(y1 = rbeta(n2, 2, 4), y2 = rbeta(n2, 3, 3))
  bw_2d_default <- npcdensbw(
    xdat = x2,
    ydat = y2,
    bws = c(0.16, 0.18, 0.22),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lc",
    cxkerbound = "range",
    cykerbound = "range",
    cvls.quadrature.points = c(100L, 50L)
  )
  bw_2d_coarse <- bw_2d_default
  bw_2d_coarse$cvls.quadrature.points <- c(81L, 17L)

  obj_2d_default <- npRmpi:::.npcdensbw_eval_only(x2, y2, bw_2d_default)$objective
  obj_2d_coarse <- npRmpi:::.npcdensbw_eval_only(x2, y2, bw_2d_coarse)$objective

  expect_true(is.finite(obj_2d_default))
  expect_true(is.finite(obj_2d_coarse))
  expect_gt(abs(obj_2d_default - obj_2d_coarse), 1e-10)
  expect_identical(bw_2d_default$cvls.quadrature.grid, "uniform")
  expect_error(
    npcdensbw(
      xdat = x2,
      ydat = y2,
      bws = c(0.16, 0.18, 0.22),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      regtype = "lc",
      cxkerbound = "range",
      cykerbound = "range",
      cvls.quadrature.grid = "hybrid"
    ),
    "scalar continuous responses"
  )
})
