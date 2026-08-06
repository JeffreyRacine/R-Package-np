test_that("rgl presentation defaults are centralized and immutable", {
  defaults <- getFromNamespace(".np_plot_rgl_defaults", "np")

  first <- defaults()
  second <- defaults()

  expect_identical(first$fov, 55)
  expect_identical(first$max.pixel.ratio, 2)
  expect_identical(first$legend.magnify, 2)
  expect_identical(first$windowRect, c(900, 100, 1540, 740))

  first$fov <- 1
  expect_identical(second$fov, 55)
})

test_that("rgl high-DPI adapter preserves the serialized scene and delegates upstream", {
  skip_if_not(suppressWarnings(requireNamespace("rgl", quietly = TRUE)))
  skip_if_not(suppressWarnings(requireNamespace("htmlwidgets", quietly = TRUE)))

  old.options <- options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
  on.exit(options(old.options), add = TRUE)
  old.environment <- Sys.getenv("RGL_USE_NULL", unset = NA_character_)
  Sys.setenv(RGL_USE_NULL = "TRUE")
  on.exit({
    if (is.na(old.environment)) {
      Sys.unsetenv("RGL_USE_NULL")
    } else {
      Sys.setenv(RGL_USE_NULL = old.environment)
    }
  }, add = TRUE)

  rgl::open3d(useNULL = TRUE, silent = TRUE)
  on.exit(try(rgl::close3d(silent = TRUE), silent = TRUE), add = TRUE)
  rgl::points3d(0, 0, 0)
  source.widget <- rgl::rglwidget(x = rgl::scene3d())

  adapt <- getFromNamespace(".np_plot_rgl_hidpi_widget", "np")
  result <- adapt(source.widget)

  expect_identical(result$x, source.widget$x)
  expect_length(result$jsHooks$render, 1L)
  hook <- result$jsHooks$render[[1L]]$code
  expect_match(hook, "window.devicePixelRatio", fixed = TRUE)
  expect_match(hook, "upstreamResize(host)", fixed = TRUE)
  expect_match(hook, "upstreamRelMouseCoords(event)", fixed = TRUE)
  expect_match(hook, "upstreamRestartCanvas()", fixed = TRUE)
  expect_match(hook, "npHiDpiMode = 'native'", fixed = TRUE)
  expect_match(hook, "incompatible partial rgl high-DPI", fixed = TRUE)
})

test_that("rgl legend magnification retains apparent geometry", {
  scale.legend <- getFromNamespace(".np_plot_rgl_scale_legend_geometry", "np")

  result <- scale.legend(list(
    magnify = 2,
    cex = 0.8,
    pt.cex = c(0.5, 1),
    lwd = 1.25,
    box.lwd = 0.75,
    seg.len = 2
  ))

  expect_identical(result$magnify, 2)
  expect_equal(result$cex, 1.6)
  expect_equal(result$pt.cex, c(1, 2))
  expect_equal(result$lwd, 2.5)
  expect_equal(result$box.lwd, 1.5)
  expect_identical(result$seg.len, 2)

  expect_error(scale.legend(list(magnify = 0)), "finite positive")
  expect_error(scale.legend(list(magnify = Inf)), "finite positive")
})

test_that("rgl DPI and legend controls retain user precedence", {
  user.args <- getFromNamespace(".np_plot_user_args", "np")(
    list(magnify = 1.5, unsupported = TRUE),
    type = "rgl.legend3d"
  )
  expect_identical(user.args, list(magnify = 1.5))

  make.hook <- getFromNamespace(".np_plot_rgl_hidpi_hook", "np")
  expect_error(make.hook(0.5), "no smaller than one")
  expect_match(make.hook(1.5), "var maxRatio = 1.5", fixed = TRUE)
})

test_that("rgl uses the base-matched camera unless the user overrides it", {
  skip_if_not(suppressWarnings(requireNamespace("rgl", quietly = TRUE)))
  skip_if_not(suppressWarnings(requireNamespace("htmlwidgets", quietly = TRUE)))

  view.angles <- getFromNamespace(".np_plot_rgl_view_angles", "np")
  expect_identical(view.angles(0, 20), list(theta = 0, phi = -70))
  expect_identical(view.angles(25, 15), list(theta = 25, phi = 15))

  old.options <- options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
  on.exit(options(old.options), add = TRUE)
  old.environment <- Sys.getenv("RGL_USE_NULL", unset = NA_character_)
  Sys.setenv(RGL_USE_NULL = "TRUE")
  on.exit({
    if (is.na(old.environment)) {
      Sys.unsetenv("RGL_USE_NULL")
    } else {
      Sys.setenv(RGL_USE_NULL = old.environment)
    }
  }, add = TRUE)

  render <- getFromNamespace(".np_plot_render_surface_rgl", "np")
  mapped <- view.angles(0, 20)
  default.widget <- suppressWarnings(render(
    0:1,
    0:1,
    matrix(0, 2, 2),
    "x",
    "y",
    "z",
    "",
    mapped$theta,
    mapped$phi
  ))
  override.widget <- suppressWarnings(render(
    0:1,
    0:1,
    matrix(0, 2, 2),
    "x",
    "y",
    "z",
    "",
    mapped$theta,
    mapped$phi,
    view3d.args = list(fov = 45, zoom = 0.9)
  ))

  subscene <- function(widget) {
    Filter(
      function(object) identical(object$type, "subscene"),
      widget$x$objects
    )[[1L]]
  }

  expect_identical(subscene(default.widget)$par3d$FOV, 55)
  expect_identical(subscene(override.widget)$par3d$FOV, 45)
  expect_equal(subscene(override.widget)$par3d$zoom, 0.9, tolerance = 1e-7)
})
