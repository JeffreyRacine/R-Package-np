test_that("copula render matrices follow requested coordinates without reordering fits", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(-1, 1, length.out = 24), y = sin(seq_len(24)))
  u <- data.frame(u1 = c(.8, .2, .5), u2 = c(.7, .3, .5))
  b <- npudistbw(dat = d, bws = c(.3, .4), bandwidth.compute = FALSE)
  fit <- npcopula(b, data = d, u = u, n.quasi.inv = 40L)
  original <- as.data.frame(fit)
  sorted <- npcopula(b, data = d, u = as.data.frame(lapply(u, sort)), n.quasi.inv = 40L)
  grid <- npRmpi:::.npcopula_grid_eval(fit)
  expect_identical(grid$u1, sort(u$u1))
  expect_identical(grid$u2, sort(u$u2))
  expect_identical(grid$z, matrix(fitted(sorted), 3L, 3L))
  expect_identical(grid$xgrid, npRmpi:::.npcopula_eval_xgrid(fit))
  expect_identical(as.data.frame(fit), original)
  old.fit <- fit; old.fit$u.grid <- NULL
  expect_identical(npRmpi:::.npcopula_grid_eval(old.fit)$z, grid$z)
  repeated <- npcopula(b, data = d, u = data.frame(u1 = c(.2, .2), u2 = c(.3, .7)),
    n.quasi.inv = 40L)
  expect_error(npRmpi:::.npcopula_grid_eval(repeated), "not rectangular")
  malformed <- fit; malformed$u.grid[1L, ] <- malformed$u.grid[2L, ]
  expect_error(npRmpi:::.npcopula_grid_eval(malformed), "not rectangular")
  endpoints <- npcopula(b, data = d,
    u = data.frame(u1 = c(1, 0, .5), u2 = c(1, 0, .5)), n.quasi.inv = 40L)
  endpoint.grid <- npRmpi:::.npcopula_grid_eval(endpoints)
  expect_identical(endpoint.grid$u1, c(0, .5, 1))
  expect_identical(endpoint.grid$u2, c(0, .5, 1))
  expect_identical(endpoint.grid$z, matrix(fitted(endpoints)[
    order(endpoints$u.grid[[2L]], endpoints$u.grid[[1L]])], 3L, 3L))
})

test_that("all copula renderers and uncertainty layers share the coordinate map", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(-1, 1, length.out = 24), y = sin(seq_len(24)))
  u <- data.frame(u1 = c(.8, .2, .5), u2 = c(.7, .3, .5))
  b <- npudistbw(dat = d, bws = c(.3, .4), bandwidth.compute = FALSE)
  fit <- npcopula(b, data = d, u = u, n.quasi.inv = 40L)
  original <- as.data.frame(fit)
  permutation <- order(fit$u.grid[[2L]], fit$u.grid[[1L]])
  expected <- matrix(fitted(fit)[permutation], 3L, 3L)
  captured <- new.env(parent = emptyenv())
  captured$frames <- list(); captured$errors <- NULL
  capture.frame <- function(x, y, z, ...) {
    captured$frames[[length(captured$frames) + 1L]] <- list(x = x, y = y, z = z)
    invisible(NULL)
  }
  capture.errors <- function(lerr, herr, lerr.all, herr.all, ...) {
    captured$errors <- list(lower = lerr, upper = herr,
      lower.all = lerr.all, upper.all = herr.all)
  }
  local_mocked_bindings(contour = capture.frame, image = capture.frame,
    .package = "graphics")
  local_mocked_bindings(
    .np_plot_render_surface_base_frame = function(persp.args, ...) {
      do.call(capture.frame, persp.args); diag(4L)
    },
    .np_plot_draw_error_wireframes_persp = capture.errors,
    .np_plot_render_surface_rgl = function(x, y, z, draw.extras, ...) {
      capture.frame(x, y, z); draw.extras(); invisible(NULL)
    },
    .np_plot_error_surfaces_rgl = capture.errors,
    .package = "npRmpi")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  for (view in c("contour", "image")) {
    captured$frames <- list()
    out <- plot(fit, view = view, output = "plot-data")
    expect_identical(out, original)
    expect_identical(captured$frames[[1L]]$z, expected)
    expect_identical(captured$frames[[1L]]$x, sort(u$u1))
    expect_identical(captured$frames[[1L]]$y, sort(u$u2))
  }
  for (renderer in c("base", "rgl")) {
    captured$frames <- list(); captured$errors <- NULL
    out <- plot(fit, view = "fixed", renderer = renderer, output = "plot-data",
      errors = "asymptotic", band = "all", legend = FALSE)
    expect_identical(captured$frames[[1L]]$z, expected)
    expect_identical(out[names(original)], original)
    expect_identical(captured$errors$lower, matrix(out$lower[permutation], 3L, 3L))
    expect_identical(captured$errors$upper, matrix(out$upper[permutation], 3L, 3L))
    for (name in names(captured$errors$lower.all)) {
      expect_identical(captured$errors$lower.all[[name]],
        matrix(out[[paste0(name, ".lower")]][permutation], 3L, 3L))
      expect_identical(captured$errors$upper.all[[name]],
        matrix(out[[paste0(name, ".upper")]][permutation], 3L, 3L))
    }
  }
  captured$frames <- list()
  out <- plot(fit, view = "all", output = "plot-data")
  expect_identical(out$copula, original)
  expect_identical(captured$frames[[1L]]$z, expected)
  expect_identical(captured$frames[[2L]]$z, expected)
  expect_identical(captured$frames[[3L]]$z, matrix(out$density$copula[permutation], 3L, 3L))
  expect_identical(as.data.frame(fit), original)
  captured$frames <- list(); captured$errors <- NULL
  set.seed(121)
  out <- plot(fit, view = "fixed", output = "plot-data", errors = "bootstrap",
    B = 3L, band = "pmzsd", legend = FALSE)
  seed <- .Random.seed
  expect_identical(captured$frames[[1L]]$z, expected)
  expect_identical(captured$errors$lower, matrix(out$lower[permutation], 3L, 3L))
  expect_identical(captured$errors$upper, matrix(out$upper[permutation], 3L, 3L))
  set.seed(121)
  same <- plot(fit, output = "data", errors = "bootstrap", B = 3L, band = "pmzsd")
  expect_identical(same, out)
  expect_identical(.Random.seed, seed)
})
