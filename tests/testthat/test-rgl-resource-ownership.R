# Unique contract: exact renderer-owned device cleanup and caller restoration.
# Existing RGL module/HiDPI tests cover scene controls and widget presentation;
# these bounded 2x2 surfaces add error/resource coverage without MPI or searches.
.rgl_resource_package <- "np"

.rgl_resource_scope <- function(run) {
  skip_if_not(suppressWarnings(requireNamespace("rgl", quietly = TRUE)))
  scope <- new.env(parent = emptyenv())
  scope$opened <- integer()
  scope$initial.current <- as.integer(rgl::cur3d())
  open <- rgl::open3d
  on.exit({
    for (device in rev(scope$opened)) {
      if (device %in% rgl::rgl.dev.list())
        try(rgl::close3d(dev = device, silent = TRUE), silent = TRUE)
    }
    if (scope$initial.current != 0L &&
        scope$initial.current %in% rgl::rgl.dev.list())
      try(rgl::set3d(scope$initial.current, silent = TRUE), silent = TRUE)
  }, add = TRUE)
  testthat::with_mocked_bindings(
    run(scope),
    open3d = function(...) {
      device <- open(...)
      scope$opened <- c(scope$opened, as.integer(device))
      device
    },
    .package = "rgl"
  )
}

.rgl_resource_render <- function(...) {
  render <- getFromNamespace(".np_plot_render_surface_rgl", .rgl_resource_package)
  render(x = c(-1, 2), y = c(0, 4),
         z = matrix(c(.2, .6, .3, .9), 2L, 2L),
         xlab = "x", ylab = "y", zlab = "z", main = "",
         theta = 0, phi = 20, ...)
}

.rgl_resource_ids <- function() as.integer(rgl::rgl.dev.list())

.rgl_resource_state <- function() {
  list(options = options()[c("rgl.useNULL", "rgl.printRglwidget")],
       environment = Sys.getenv("RGL_USE_NULL", unset = NA_character_))
}

test_that("RGL01 ordinary rendering errors release the exact opened device", {
  .rgl_resource_scope(function(scope) {
    before <- .rgl_resource_ids()
    current <- as.integer(rgl::cur3d())
    state <- .rgl_resource_state()
    expect_error(.rgl_resource_render(col = "not-a-color"),
                 "renderer='rgl' failed to draw the surface (invalid color name 'not-a-color')",
                 fixed = TRUE)
    expect_length(scope$opened, 1L)
    expect_false(scope$opened[[1L]] %in% .rgl_resource_ids())
    expect_identical(.rgl_resource_ids(), before)
    expect_identical(as.integer(rgl::cur3d()), current)
    expect_identical(.rgl_resource_state(), state)
  })
})

test_that("RGL02 success closes before widget creation and restores prior current device", {
  .rgl_resource_scope(function(scope) {
    older <- as.integer(rgl::open3d(useNULL = TRUE, silent = TRUE))
    newer <- as.integer(rgl::open3d(useNULL = TRUE, silent = TRUE))
    rgl::set3d(older, silent = TRUE)
    before <- .rgl_resource_ids()
    state <- .rgl_resource_state()
    widget <- list(marker = "unchanged widget")
    rendered <- NULL
    testthat::with_mocked_bindings(
      testthat::with_mocked_bindings({
        capture.output(rendered <- .rgl_resource_render(draw.extras = function() {
          scope$plot.device <- as.integer(rgl::cur3d())
        }))
      }, rglwidget = function(x, ...) {
        scope$scene <- x
        scope$ids.at.widget <- .rgl_resource_ids()
        widget
      }, .package = "rgl"),
      .np_plot_rgl_hidpi_widget = function(widget, ...) widget,
      .package = .rgl_resource_package
    )
    expect_identical(scope$plot.device, tail(scope$opened, 1L))
    expect_false(scope$plot.device %in% scope$ids.at.widget)
    expect_true(is.list(scope$scene))
    expect_identical(rendered, widget)
    expect_identical(.rgl_resource_ids(), before)
    expect_true(all(c(older, newer) %in% .rgl_resource_ids()))
    expect_identical(as.integer(rgl::cur3d()), older)
    expect_identical(.rgl_resource_state(), state)
  })
})

test_that("RGL03 callback devices survive cleanup after a callback error", {
  .rgl_resource_scope(function(scope) {
    older <- as.integer(rgl::open3d(useNULL = TRUE, silent = TRUE))
    newer <- as.integer(rgl::open3d(useNULL = TRUE, silent = TRUE))
    rgl::set3d(older, silent = TRUE)
    before <- .rgl_resource_ids()
    state <- .rgl_resource_state()
    expect_error(.rgl_resource_render(draw.extras = function() {
      scope$plot.device <- as.integer(rgl::cur3d())
      scope$callback.device <- as.integer(rgl::open3d(useNULL = TRUE, silent = TRUE))
      stop("resource callback failure")
    }), "renderer='rgl' failed to draw the surface (resource callback failure)",
    fixed = TRUE)
    expect_false(scope$plot.device %in% .rgl_resource_ids())
    expect_true(scope$callback.device %in% .rgl_resource_ids())
    expect_setequal(.rgl_resource_ids(), c(before, scope$callback.device))
    expect_true(all(c(older, newer) %in% .rgl_resource_ids()))
    expect_identical(as.integer(rgl::cur3d()), older)
    expect_identical(.rgl_resource_state(), state)
  })
})

test_that("RGL04 late widget errors preserve caller devices after early closure", {
  .rgl_resource_scope(function(scope) {
    older <- as.integer(rgl::open3d(useNULL = TRUE, silent = TRUE))
    newer <- as.integer(rgl::open3d(useNULL = TRUE, silent = TRUE))
    rgl::set3d(older, silent = TRUE)
    before <- .rgl_resource_ids()
    state <- .rgl_resource_state()
    testthat::with_mocked_bindings(
      expect_error(.rgl_resource_render(draw.extras = function() {
        scope$plot.device <- as.integer(rgl::cur3d())
      }), "renderer='rgl' failed to draw the surface (late widget failure)",
      fixed = TRUE),
      rglwidget = function(...) {
        scope$ids.at.widget <- .rgl_resource_ids()
        stop("late widget failure")
      }, .package = "rgl"
    )
    expect_false(scope$plot.device %in% scope$ids.at.widget)
    expect_identical(.rgl_resource_ids(), before)
    expect_true(all(c(older, newer) %in% .rgl_resource_ids()))
    expect_identical(as.integer(rgl::cur3d()), older)
    expect_identical(.rgl_resource_state(), state)
  })
})
