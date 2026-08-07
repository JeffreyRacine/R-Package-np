surface_ns <- asNamespace(if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np")
surface_internal <- function(name) get(name, envir = surface_ns, inherits = FALSE)

test_that("surface palette mapping has one exact shared contract", {
  surface_palette <- surface_internal(".np_plot_surface_palette")
  surface_values <- surface_internal(".np_plot_surface_palette_values")
  rgl_colors <- surface_internal(".np_plot_rgl_surface_colors")
  base_colors <- surface_internal(".np_plot_persp_surface_colors")
  z <- matrix(c(-1, 0, 1, 2, 3, 5), nrow = 2L)
  palette <- surface_palette(7L)
  z.range <- range(z, finite = TRUE)

  mapped <- surface_values(z, z.range, 7L)
  scaled <- 1L + floor(6 * (z - z.range[1L]) / diff(z.range))

  expect_identical(mapped, palette[scaled])
  expect_identical(
    rgl_colors(z, num.colors = 7L),
    matrix(mapped, nrow = nrow(z), ncol = ncol(z))
  )
  expect_identical(
    base_colors(z, col = c("red", "blue")),
    c("red", "blue")
  )
})

test_that("surface grid and opacity defaults are centrally owned", {
  plot_color <- surface_internal(".np_plot_color")
  plot_lwd <- surface_internal(".np_plot_lwd")
  surface_alpha <- surface_internal(".np_plot_surface_alpha")
  merge_args <- surface_internal(".np_plot_merge_override_args")
  expect_identical(plot_color("surface_grid"), "gray")
  expect_identical(plot_lwd("surface_grid", 1), 1)
  expect_identical(surface_alpha("base"), 0.5)
  expect_identical(surface_alpha("rgl"), 0.6)
  expect_error(plot_color("support_grid"), "unknown plot color role")

  merged <- merge_args(
    list(side = c("x", "y+", "z"), col = plot_color("surface_grid"), lwd = 1),
    list(col = "red", lwd = 2)
  )
  expect_identical(merged$col, "red")
  expect_identical(merged$lwd, 2)
})

test_that("canonical base frame preserves the persp transform", {
  draw_grid <- surface_internal(".np_plot_draw_box_grid_persp")
  render_frame <- surface_internal(".np_plot_render_surface_base_frame")
  x <- seq(0, 1, length.out = 5L)
  y <- seq(-1, 1, length.out = 6L)
  z <- outer(x, y, function(a, b) a + b^2)
  zlim <- range(z, finite = TRUE)
  args <- list(x = x, y = y, z = z, zlim = zlim,
               theta = 15, phi = 25, ticktype = "detailed")

  file.one <- tempfile(fileext = ".pdf")
  file.two <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(file.one, file.two)), add = TRUE)

  grDevices::pdf(file.one)
  legacy <- do.call(graphics::persp, args)
  draw_grid(range(x), range(y), zlim, legacy)
  grDevices::dev.off()

  grDevices::pdf(file.two)
  canonical <- render_frame(args, range(x), range(y), zlim)
  grDevices::dev.off()

  expect_identical(canonical, legacy)
})
