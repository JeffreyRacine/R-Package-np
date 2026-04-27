library(npRmpi)

compact_window_route_fixture <- function() {
  set.seed(20260423)
  n <- 24L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = 0.4 + 0.8 * x$x + rchisq(n, df = 4))
  list(x = x, y = y)
}

compact_window_span_bounds <- function(vals, lb, ub, factor = 2) {
  rng <- range(vals)
  span <- diff(rng)
  if (!is.finite(span) || span <= 0)
    span <- 1
  c(
    lb = if (is.finite(lb)) lb else rng[1L] - factor * span,
    ub = if (is.finite(ub)) ub else rng[2L] + factor * span
  )
}

compact_window_clone_bw <- function(bw, bound, lb = NULL, ub = NULL) {
  if (bound == "range") {
    yvals <- as.numeric(bw$ydati)
    lb <- min(yvals)
    ub <- max(yvals)
  }
  if (bound == "none") {
    lb <- -Inf
    ub <- Inf
  }

  out <- npRmpi:::conbandwidth(
    xbw = bw$xbw,
    ybw = bw$ybw,
    bwmethod = bw$method,
    bwscaling = bw$scaling,
    bwtype = bw$type,
    cxkertype = bw$cxkertype,
    cxkerorder = bw$cxkerorder,
    cxkerbound = bw$cxkerbound,
    cxkerlb = bw$cxkerlb,
    cxkerub = bw$cxkerub,
    uxkertype = bw$uxkertype,
    oxkertype = bw$oxkertype,
    cykertype = bw$cykertype,
    cykerorder = bw$cykerorder,
    cykerbound = bound,
    cykerlb = lb,
    cykerub = ub,
    uykertype = bw$uykertype,
    oykertype = bw$oykertype,
    fval = bw$fval,
    ifval = bw$ifval,
    num.feval = bw$num.feval,
    num.feval.fast = bw$num.feval.fast,
    fval.history = bw$fval.history,
    eval.history = bw$eval.history,
    invalid.history = bw$invalid.history,
    nobs = bw$nobs,
    xdati = bw$xdati,
    ydati = bw$ydati,
    xnames = bw$xnames,
    ynames = bw$ynames,
    sfactor = bw$sfactor,
    bandwidth = bw$bandwidth,
    rows.omit = bw$rows.omit,
    nconfac = bw$nconfac,
    ncatfac = bw$ncatfac,
    sdev = bw$sdev,
    bandwidth.compute = TRUE,
    timing = bw$timing,
    total.time = bw$total.time,
    regtype = bw$regtype,
    pregtype = bw$pregtype,
    basis = bw$basis,
    degree = bw$degree,
    bernstein.basis = bw$bernstein.basis,
    regtype.engine = bw$regtype.engine,
    basis.engine = bw$basis.engine,
    degree.engine = bw$degree.engine,
    bernstein.basis.engine = bw$bernstein.basis.engine
  )
  out$scale.factor.search.lower <- bw$scale.factor.search.lower
  out$cvls.quadrature.grid <- bw$cvls.quadrature.grid
  out$cvls.quadrature.extend.factor <- bw$cvls.quadrature.extend.factor
  out$cvls.quadrature.points <- bw$cvls.quadrature.points
  out
}

test_that("one-sided fixed infinite bounds use the configured span surrogate", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- compact_window_route_fixture()

  bw_upper <- npcdensbw(
    ydat = dat$y,
    xdat = dat$x,
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
    cvls.quadrature.points = c(31L, 17L)
  )

  bw_lower <- npcdensbw(
    ydat = dat$y,
    xdat = dat$x,
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
    cykerlb = -Inf,
    cykerub = max(dat$y$y) + 0.25,
    cvls.quadrature.grid = "uniform",
    cvls.quadrature.points = c(31L, 17L)
  )
  bw_upper_span2 <- bw_upper
  bw_upper_span2$cvls.quadrature.extend.factor <- 2
  bw_lower_span2 <- bw_lower
  bw_lower_span2$cvls.quadrature.extend.factor <- 2

  upper_span1 <- compact_window_span_bounds(dat$y$y, 0, Inf, factor = 1)
  upper_span2 <- compact_window_span_bounds(dat$y$y, 0, Inf, factor = 2)
  lower_span1 <- compact_window_span_bounds(dat$y$y, -Inf, max(dat$y$y) + 0.25, factor = 1)
  lower_span2 <- compact_window_span_bounds(dat$y$y, -Inf, max(dat$y$y) + 0.25, factor = 2)

  upper_obj <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_upper)$objective
  lower_obj <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_lower)$objective
  upper_obj_factor2 <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_upper_span2)$objective
  lower_obj_factor2 <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_lower_span2)$objective

  upper_obj_span1 <- npRmpi:::.npcdensbw_eval_only(
    dat$x, dat$y, compact_window_clone_bw(bw_upper, "fixed", upper_span1[["lb"]], upper_span1[["ub"]])
  )$objective
  upper_obj_span2 <- npRmpi:::.npcdensbw_eval_only(
    dat$x, dat$y, compact_window_clone_bw(bw_upper, "fixed", upper_span2[["lb"]], upper_span2[["ub"]])
  )$objective

  lower_obj_span1 <- npRmpi:::.npcdensbw_eval_only(
    dat$x, dat$y, compact_window_clone_bw(bw_lower, "fixed", lower_span1[["lb"]], lower_span1[["ub"]])
  )$objective
  lower_obj_span2 <- npRmpi:::.npcdensbw_eval_only(
    dat$x, dat$y, compact_window_clone_bw(bw_lower, "fixed", lower_span2[["lb"]], lower_span2[["ub"]])
  )$objective

  expect_equal(upper_obj, upper_obj_span1, tolerance = 1e-12)
  expect_equal(upper_obj_factor2, upper_obj_span2, tolerance = 1e-12)
  expect_gt(abs(upper_obj - upper_obj_factor2), 1e-8)

  expect_equal(lower_obj, lower_obj_span1, tolerance = 1e-12)
  expect_equal(lower_obj_factor2, lower_obj_span2, tolerance = 1e-12)
  expect_gt(abs(lower_obj - lower_obj_factor2), 1e-8)
})

test_that("explicit fixed [-Inf, Inf] survives and uses the configured span surrogate", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- compact_window_route_fixture()

  bw_two_inf <- npcdensbw(
    ydat = dat$y,
    xdat = dat$x,
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
    cykerlb = -Inf,
    cykerub = Inf,
    cvls.quadrature.grid = "uniform",
    cvls.quadrature.points = c(31L, 17L)
  )
  bw_two_inf_span2 <- bw_two_inf
  bw_two_inf_span2$cvls.quadrature.extend.factor <- 2

  span1 <- compact_window_span_bounds(dat$y$y, -Inf, Inf, factor = 1)
  span2 <- compact_window_span_bounds(dat$y$y, -Inf, Inf, factor = 2)
  obj_two_inf <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_two_inf)$objective
  obj_two_inf_factor2 <- npRmpi:::.npcdensbw_eval_only(dat$x, dat$y, bw_two_inf_span2)$objective
  obj_span1 <- npRmpi:::.npcdensbw_eval_only(
    dat$x,
    dat$y,
    compact_window_clone_bw(bw_two_inf, "fixed", span1[["lb"]], span1[["ub"]])
  )$objective
  obj_span2 <- npRmpi:::.npcdensbw_eval_only(
    dat$x,
    dat$y,
    compact_window_clone_bw(bw_two_inf, "fixed", span2[["lb"]], span2[["ub"]])
  )$objective

  expect_identical(as.character(bw_two_inf$cykerbound), "fixed")
  expect_true(is.infinite(bw_two_inf$cykerlb[which(bw_two_inf$iycon)][1L]))
  expect_true(is.infinite(bw_two_inf$cykerub[which(bw_two_inf$iycon)][1L]))
  expect_equal(obj_two_inf, obj_span1, tolerance = 1e-12)
  expect_equal(obj_two_inf_factor2, obj_span2, tolerance = 1e-12)
  expect_gt(abs(obj_two_inf - obj_two_inf_factor2), 1e-8)
})
