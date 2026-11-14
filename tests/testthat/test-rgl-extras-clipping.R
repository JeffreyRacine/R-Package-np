# Unique contract: all layers share the existing rgl clipregion; geometry is retained.
# Deterministic rendering only; no estimator, simulation, optimizer or MPI search.
test_that("rgl extras inherit the surface clipping owner without changing vertices", {
  skip_if_not_installed("rgl")
  old <- options(rgl.useNULL = TRUE)
  on.exit(options(old), add = TRUE)
  old.env <- Sys.getenv("RGL_USE_NULL", unset = NA_character_)
  Sys.setenv(RGL_USE_NULL = "TRUE")
  on.exit(if (is.na(old.env)) Sys.unsetenv("RGL_USE_NULL") else
    Sys.setenv(RGL_USE_NULL = old.env), add = TRUE)
  before <- as.integer(rgl::rgl.dev.list())
  axis <- c(0, .5, 1)
  z <- outer(axis, axis, "+")/2
  extra.z <- z + 2
  points <- matrix(c(-1, 2, -2, 3, -3, 4), ncol = 3L)
  render <- .np_plot_render_surface_rgl
  pkg <- getNamespaceName(environment(render))
  original.persp <- rgl::persp3d
  original.grid <- rgl::grid3d
  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(persp3d = function(...) {
    result <- original.persp(...)
    seen$surface.ids <- result
    seen$parent <- rgl::currentSubscene3d()
    result
  }, grid3d = function(...) {
    seen$grid.owner <- rgl::currentSubscene3d()
    seen$bbox <- as.numeric(rgl::par3d("bbox"))
    original.grid(...)
  }, rglwidget = function(x, ...) list(scene = x), .package = "rgl")
  testthat::local_mocked_bindings(.np_plot_rgl_hidpi_widget = function(widget, ...) widget,
                                .package = pkg)
  cases <- list(none = list(), z = list(zlim = c(0, 1)),
                xy = list(xlim = c(.1, .9), ylim = c(.2, .8)),
                all = list(xlim = c(.1, .9), ylim = c(.2, .8), zlim = c(0, 1)),
                override = list(zlim = c(0, 4), persp3d.args = list(zlim = c(0, 1))))
  for (name in names(cases)) {
    args <- c(list(x = axis, y = axis, z = z, xlab = "x", ylab = "y", zlab = "z",
                   main = "", theta = 0, phi = 20, draw.extras = function() {
                     seen$extra.owner <- rgl::currentSubscene3d()
                     id <- rgl::points3d(points)
                     seen$points <- rgl::rgl.attrib(id, "vertices")
                     id <- rgl::surface3d(axis, axis, extra.z)
                     seen$surface <- rgl::rgl.attrib(id, "vertices")
                   }), cases[[name]])
    capture.output(do.call(render, args))
    clip <- unname(seen$surface.ids[names(seen$surface.ids) == "clipregion"])
    expect_identical(seen$grid.owner, seen$parent)
    expect_equal(unname(seen$points), points, tolerance = 0)
    expect_equal(sort(seen$surface[, 3L]), sort(as.vector(extra.z)), tolerance = 0)
    expect_identical(as.integer(rgl::rgl.dev.list()), before)
    if (name == "none") {
      expect_length(clip, 0L)
      expect_identical(seen$extra.owner, seen$parent)
      expect_equal(seen$bbox, c(-1, 2, -2, 3, -3, 4), tolerance = 0)
    } else {
      expect_length(clip, 1L)
      expect_equal(as.integer(seen$extra.owner), as.integer(clip))
      if (name %in% c("z", "all", "override")) expect_equal(seen$bbox[5:6], c(0, 1), tolerance = 1e-6)
      if (name %in% c("xy", "all")) expect_equal(seen$bbox[1:4], c(.1, .9, .2, .8), tolerance = 1e-6)
    }
  }
})


test_that("rgl all-band legends retain the view background owner", {
  skip_if_not_installed("rgl")
  withr::local_options(rgl.useNULL = TRUE)
  withr::local_envvar(RGL_USE_NULL = "TRUE")
  before <- as.integer(rgl::rgl.dev.list())
  axis <- c(0, .5, 1)
  z <- outer(axis, axis, "+") / 2
  render <- .np_plot_render_surface_rgl
  pkg <- getNamespaceName(environment(render))
  seen <- new.env(parent = emptyenv())
  surfaces <- function(scene) {
    ids <- rgl::ids3d(type = "shapes", subscene = scene)
    ids$id[ids$type == "surface"]
  }
  original.legend <- rgl::legend3d
  testthat::local_mocked_bindings(legend3d = function(...) {
    seen$legend.owner <- rgl::currentSubscene3d()
    original.legend(...)
  }, rglwidget = function(x, ...) list(scene = x), .package = "rgl")
  testthat::local_mocked_bindings(.np_plot_rgl_hidpi_widget = function(widget, ...) widget,
                                .package = pkg)
  for (limit in list(NULL, c(0, 1))) {
    capture.output(render(x = axis, y = axis, z = z, zlim = limit,
      xlab = "x", ylab = "y", zlab = "z", main = "", theta = 0, phi = 20,
      draw.extras = function() {
        current <- rgl::currentSubscene3d()
        root <- rgl::subsceneInfo("root")$id
        before.child <- surfaces(current)
        before.root <- surfaces(root)
        lower <- list(pointwise = z - 2, simultaneous = z - 3, bonferroni = z - 4)
        upper <- list(pointwise = z + 2, simultaneous = z + 3, bonferroni = z + 4)
        .np_plot_error_surfaces_rgl(axis, axis, "all", "asymptotic",
          lerr.all = lower, herr.all = upper)
        seen$band.ids <- setdiff(surfaces(current),
                                before.child)
        seen$root.band.ids <- setdiff(surfaces(root),
                                     before.root)
        seen$band.z <- sort(unlist(lapply(seen$band.ids, function(id)
          rgl::rgl.attrib(id, "vertices")[, 3L]), use.names = FALSE))
        seen$expected.z <- sort(unlist(c(lower, upper), use.names = FALSE))
        seen$root <- root
        seen$current <- current
        seen$after <- rgl::currentSubscene3d()
        seen$root.backgrounds <- nrow(rgl::ids3d(type = "background", subscene = root))
        seen$child.backgrounds <- if (!identical(current, root))
          nrow(rgl::ids3d(type = "background", subscene = current)) else 0L
      }))
    expect_identical(seen$legend.owner, seen$root)
    expect_identical(seen$after, seen$current)
    expect_equal(seen$root.backgrounds, 1L)
    expect_equal(seen$child.backgrounds, 0L)
    expect_length(seen$band.ids, 6L)
    expect_equal(seen$band.z, seen$expected.z, tolerance = 0)
    if (!is.null(limit)) expect_length(seen$root.band.ids, 0L)
    expect_identical(as.integer(rgl::rgl.dev.list()), before)
  }
})


test_that("rgl legend errors restore the clipped geometry owner", {
  skip_if_not_installed("rgl")
  withr::local_options(rgl.useNULL = TRUE)
  withr::local_envvar(RGL_USE_NULL = "TRUE")
  before <- as.integer(rgl::rgl.dev.list())
  axis <- c(0, .5, 1)
  z <- outer(axis, axis, "+") / 2
  render <- .np_plot_render_surface_rgl
  pkg <- getNamespaceName(environment(render))
  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(legend3d = function(...) {
    stop("controlled legend failure", call. = FALSE)
  }, rglwidget = function(x, ...) list(scene = x), .package = "rgl")
  testthat::local_mocked_bindings(.np_plot_rgl_hidpi_widget = function(widget, ...) widget,
                                .package = pkg)
  capture.output(render(x = axis, y = axis, z = z, zlim = c(0, 1),
    xlab = "x", ylab = "y", zlab = "z", main = "", theta = 0, phi = 20,
    draw.extras = function() {
      seen$incoming <- rgl::currentSubscene3d()
      seen$root <- rgl::subsceneInfo("root")$id
      seen$error <- tryCatch({
        .np_plot_error_surfaces_rgl(axis, axis, "all", "asymptotic",
          lerr.all = list(pointwise = z - 2), herr.all = list(pointwise = z + 2))
        NULL
      }, error = conditionMessage)
      seen$after <- rgl::currentSubscene3d()
      point <- rgl::points3d(.5, .5, .5)
      seen$point.in.child <- point %in% rgl::ids3d(type = "shapes", subscene = seen$incoming)$id
      seen$point.in.root <- point %in% rgl::ids3d(type = "shapes", subscene = seen$root)$id
    }))
  expect_identical(seen$error, "controlled legend failure")
  expect_false(identical(seen$incoming, seen$root))
  expect_identical(seen$after, seen$incoming)
  expect_true(seen$point.in.child)
  expect_false(seen$point.in.root)
  expect_identical(as.integer(rgl::rgl.dev.list()), before)
})
