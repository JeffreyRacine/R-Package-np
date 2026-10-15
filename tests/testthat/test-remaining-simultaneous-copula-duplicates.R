test_that("copula display quotient preserves original rows and rejects conflicts", {
  grid.eval <- getFromNamespace(".npcopula_grid_eval", "npRmpi")
  grid.matrix <- getFromNamespace(".npcopula_grid_matrix", "npRmpi")
  u <- expand.grid(u1 = c(.5, .2, .5), u2 = c(.25, .5, .75))
  value <- u$u1 + 2 * u$u2
  dat <- data.frame(copula = value, u, x = u$u1 - 1, y = u$u2 + 1)
  fit <- structure(list(eval = dat, copula = value, xnames = c("x", "y"),
                         evaluation = "grid", grid.dim = c(3L, 3L), u.grid = u),
                    class = "npcopula")
  original <- serialize(fit, NULL)
  grid <- grid.eval(fit)
  expect_identical(grid$xgrid, dat[c("x", "y")])
  expect_length(grid$row.order, 9L)
  expect_identical(dim(grid$z), c(2L, 3L))
  expect_identical(grid$z, outer(c(.2, .5), c(.25, .5, .75), function(x,y) x + 2*y))
  expect_identical(serialize(fit, NULL), original)
  lower <- value - .1
  missing <- u$u1 == .5 & u$u2 == .5
  lower[missing] <- NA_real_
  expect_identical(which(is.na(grid.matrix(lower, grid))), 4L)
  conflict <- lower
  conflict[which(duplicated(grid$display.cell))[1L]] <- 17
  expect_error(grid.matrix(conflict, grid), "conflicting surface values", fixed = TRUE)
  conflict <- lower
  conflict[which(missing)[1L]] <- 0
  expect_error(grid.matrix(conflict, grid), "conflicting surface values", fixed = TRUE)
  incomplete <- fit
  incomplete$u.grid[1L, 1L] <- .9
  expect_error(grid.eval(incomplete), "not rectangular", fixed = TRUE)
})

test_that("repeated copula simultaneous targets are computed once and expanded exactly", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  set.seed(20261015)
  n <- 80L
  # Preserve the discovery fixture's RNG stream (including its other responses).
  xx <- sort(runif(n, -1, 1)); zz <- runif(n, -1, 1)
  yy <- .4 + sin(1.5 * xx) + rnorm(n, sd = .35)
  d <- data.frame(a = xx, b = .3 * xx + rnorm(n, sd = .7))
  pkg <- unname(getNamespaceName(environment(npcopula)))
  grid.eval <- getFromNamespace(".npcopula_grid_eval", pkg)
  grid.matrix <- getFromNamespace(".npcopula_grid_matrix", pkg)
  axis <- c(.3, .5, .7)
  unique.u <- data.frame(u1 = axis, u2 = axis)
  repeated.u <- expand.grid(u1 = axis[c(3, 1, 2)], u2 = axis[c(2, 3, 1)])
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  b <- npudensbw(dat = d, bws = c(32, 32), bwtype = "generalized_nn",
    bandwidth.compute = FALSE)
  fit <- npcopula(b, data = d, u = repeated.u, n.quasi.inv = 31L, se = TRUE)
  one <- npcopula(b, data = d, u = unique.u, n.quasi.inv = 31L, se = TRUE)
  before <- serialize(fit, NULL)
  grid <- grid.eval(fit); map <- grid$display.cell
  run <- function(x, output = "data", band = "all") {
    args <- list(x = x, errors = "asymptotic", band = band,
      view = "fixed", output = output, legend = FALSE)
    if (band %in% c("all", "simultaneous")) {
      args$gaussian.draws <- 9999L
      args$random.seed <- 42L
    }
    do.call(plot, args)
  }
  dat <- run(fit); reference <- run(one)
  expect_identical(nrow(dat), 81L)
  expect_identical(dat[names(as.data.frame(fit))], as.data.frame(fit))
  for (nm in c("simultaneous.lower", "simultaneous.upper")) {
    expect_true(all(is.finite(dat[[nm]])))
    # A separate unique-grid public request is the oracle, not a renderer
    # tolerance or an average of inconsistent duplicate endpoints.
    expect_identical(dat[[nm]], reference[[nm]][map])
    expect_identical(grid.matrix(dat[[nm]], grid), matrix(reference[[nm]], 3L, 3L))
  }
  cal <- attr(dat, "band.calibration")
  expect_identical(cal$grid.size, 9L)
  expect_identical(cal$evaluation.rows, 81L)
  expect_identical(cal$evaluation.index, map)
  expect_identical(run(fit, "plot-data"), dat)
  expect_silent(run(fit, "plot"))
  sim <- run(fit, band = "simultaneous")
  expect_identical(sim$lower, dat$simultaneous.lower)
  expect_identical(sim$upper, dat$simultaneous.upper)
  expect_identical(run(fit, "plot-data", "simultaneous"), sim)
  expect_identical(serialize(fit, NULL), before)
  # Legacy marginal intervals retain their original multiplicity, including
  # Bonferroni's original row count; the repair never recalibrates them.
  for (band in c("pointwise", "bonferroni")) {
    marginal <- run(fit, band = band)
    expect_identical(marginal$lower, dat[[paste0(band, ".lower")]])
    expect_identical(marginal$upper, dat[[paste0(band, ".upper")]])
  }
})

test_that("copula quotient rejects conflicting physical targets and exact endpoint conflicts", {
  pkg <- unname(getNamespaceName(environment(npcopula)))
  interval <- getFromNamespace(".np_plot_copula_interval", pkg)
  matrix.values <- getFromNamespace(".npcopula_grid_matrix", pkg)
  grid <- list(display.rows = c(1L, 2L), display.cell = c(1L, 2L, 1L),
    row.order = c(1L, 3L, 2L), u1 = c(.2, .8), u2 = .5)
  ex <- data.frame(x = c(0, 1, 0))
  object <- list(copula = c(.2, .8, .2))
  invoke <- function(x = ex, obj = object, se = rep(.1, 3L))
    interval(NULL, obj, data.frame(x = 1:8), x, se, .05, "all", 9999L, grid)
  bad <- ex; bad$x[3L] <- .Machine$double.eps
  expect_error(invoke(x = bad), "conflicting surface values", fixed = TRUE)
  bad <- object; bad$copula[3L] <- NA_real_
  expect_error(invoke(obj = bad), "conflicting surface values", fixed = TRUE)
  expect_error(invoke(se = c(.1, .1, .1 + .Machine$double.eps)),
    "conflicting surface values", fixed = TRUE)
  # Even a one-ULP contradiction remains an error at the unchanged renderer.
  expect_error(matrix.values(c(1, 2, 1 + .Machine$double.eps), grid),
    "conflicting surface values", fixed = TRUE)
})

test_that("copula quotient preserves unavailable payloads and avoids duplicate calibration work", {
  pkg <- unname(getNamespaceName(environment(npcopula)))
  adapter <- getFromNamespace(".np_plot_copula_interval", pkg)
  # The synthetic owner fixture sits below npRmpi's unchanged master-local
  # wrapper; the public matrix above exercises the real active-pool path.
  if (identical(pkg, "npRmpi"))
    adapter <- getFromNamespace(".np_plot_copula_interval_local", pkg)
  scope <- new.env(parent = environment(adapter)); environment(adapter) <- scope
  map <- rep(c(2L, 1L), 300L)
  grid <- list(display.rows = c(2L, 1L), display.cell = map,
    row.order = order(map), u1 = c(.2, .8), u2 = .5)
  ex <- data.frame(x = c(.1, .9)[map])
  object <- list(copula = c(.2, .8)[map])
  se <- c(.1, .2)[map]
  original <- list(err = matrix(NA_real_, 600L, 3L),
    all.err = list(pointwise = matrix(.1, 600L, 2L)))
  scope$.np_plot_copula_influence <- function(object, data, ex, se) {
    expect_identical(object$copula, c(.2, .8))
    expect_identical(ex$x, c(.1, .9))
    expect_identical(se, c(.1, .2))
    list(A = rbind(c(.1, 0), c(0, .2)), exact.diagonal = TRUE)
  }
  set.seed(18)
  got <- adapter(original, object, NULL, ex, se, .05, "all", 9999L, grid)
  expect_identical(got$calibration$grid.size, 2L)
  expect_identical(got$calibration$evaluation.rows, 600L)
  expect_identical(got$calibration$simultaneous.se, se)
  expect_identical(got$all.err$pointwise, original$all.err$pointwise)
  scope$.np_plot_copula_influence <- function(...) {
    getFromNamespace(".np_plot_joint_unavailable", pkg)("declared test limitation")
  }
  expect_warning(unavailable <- adapter(original, object, NULL, ex, se,
    .05, "all", 9999L, grid), "declared test limitation", fixed = TRUE)
  expect_identical(unavailable, original)
})
