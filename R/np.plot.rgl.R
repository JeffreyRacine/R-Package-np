.np_plot_error_surfaces_rgl <- function(x,
                                        y,
                                        plot.errors.type,
                                        lerr = NULL,
                                        herr = NULL,
                                        lerr.all = NULL,
                                        herr.all = NULL,
                                        surface3d.args = list(),
                                        legend3d.args = list()) {
  draw_one <- function(z, color) {
    if (is.null(z))
      return(invisible(FALSE))
    if (!any(is.finite(z)))
      return(invisible(FALSE))
    surf.args <- .np_plot_merge_override_args(
      list(
        x = x,
        y = y,
        z = z,
        color = color,
        alpha = 0.2,
        front = "lines",
        back = "lines",
        lit = FALSE
      ),
      surface3d.args
    )
    do.call(rgl::surface3d, surf.args)
    invisible(TRUE)
  }

  if (identical(plot.errors.type, "all") &&
      !is.null(lerr.all) &&
      !is.null(herr.all)) {
    band.cols <- .np_plot_all_band_colors()
    band.alpha <- .np_plot_all_band_alpha()
    drawn.bands <- character(0L)
    for (bn in c("pointwise", "simultaneous", "bonferroni")) {
      drawn.lower <- draw_one(
        lerr.all[[bn]],
        grDevices::adjustcolor(band.cols[[bn]], alpha.f = band.alpha[[bn]])
      )
      drawn.upper <- draw_one(
        herr.all[[bn]],
        grDevices::adjustcolor(band.cols[[bn]], alpha.f = band.alpha[[bn]])
      )
      if (isTRUE(drawn.lower) || isTRUE(drawn.upper))
        drawn.bands <- c(drawn.bands, bn)
    }
    if (!length(drawn.bands)) {
      return(invisible(FALSE))
    }
    legend3d.call <- .np_plot_merge_override_args(
      list(
        "topright",
        legend = c(pointwise = "Pointwise",
                   simultaneous = "Simultaneous",
                   bonferroni = "Bonferroni")[drawn.bands],
        col = unname(band.cols[drawn.bands]),
        lty = .np_plot_lty("solid"),
        lwd = .np_plot_lwd("band_all_surface"),
        cex = .np_plot_cex("legend"),
        bg = .np_plot_color("legend_bg"),
        bty = "n"
      ),
      legend3d.args
    )
    do.call(rgl::legend3d, legend3d.call)
    return(invisible(TRUE))
  }

  draw_one(lerr, .np_plot_color("context_wire"))
  draw_one(herr, .np_plot_color("context_wire"))
  invisible(TRUE)
}

.np_plot_rgl_view_angles <- function(theta, phi) {
  theta <- as.double(theta)[1L]
  phi <- as.double(phi)[1L]

  if (isTRUE(all.equal(theta, 0.0)) && isTRUE(all.equal(phi, 20.0))) {
    phi <- -70.0
  }

  list(theta = theta, phi = phi)
}

.np_plot_rgl_surface_colors <- function(z, col = NULL, num.colors = 1000L) {
  if (!is.null(col))
    return(col)

  z <- as.matrix(z)
  z.range <- range(z, finite = TRUE)
  palette_fun <- function(n) grDevices::hcl.colors(as.integer(n), palette = "viridis")

  if (!all(is.finite(z.range)))
    return(palette_fun(1L))

  if (isTRUE(all.equal(z.range[1L], z.range[2L]))) {
    return(matrix(palette_fun(1L), nrow = nrow(z), ncol = ncol(z)))
  }

  colorlut <- palette_fun(num.colors)
  scaled <- 1L + floor((length(colorlut) - 1L) * (z - z.range[1L]) / diff(z.range))
  scaled[!is.finite(scaled)] <- 1L
  scaled <- pmax.int(1L, pmin.int(length(colorlut), scaled))

  matrix(colorlut[scaled], nrow = nrow(z), ncol = ncol(z))
}

.np_plot_validate_renderer_request <- function(renderer,
                                               route,
                                               perspective,
                                               supported.route = TRUE,
                                               view = "fixed",
                                               gradients = FALSE,
                                               coef = FALSE,
                                               plot.errors.method = "none",
                                               plot.data.overlay = FALSE,
                                               plot.behavior = "plot",
                                               allow.plot.errors = FALSE,
                                               allow.plot.data.overlay = FALSE) {
  renderer <- .np_plot_match_renderer(renderer)
  plot.behavior <- as.character(plot.behavior)[1L]

  if (!identical(renderer, "rgl"))
    return(renderer)

  if (identical(plot.behavior, "data"))
    return(renderer)

  if (!isTRUE(perspective)) {
    stop("renderer='rgl' requires perspective=TRUE.", call. = FALSE)
  }

  if (!isTRUE(supported.route)) {
    stop(sprintf(
      "renderer='rgl' is currently implemented only for supported 2D surface routes in %s",
      route
    ), call. = FALSE)
  }

  if (isTRUE(gradients)) {
    stop("renderer='rgl' does not yet support gradients=TRUE. Use renderer='base'.",
         call. = FALSE)
  }

  if (isTRUE(coef)) {
    stop("renderer='rgl' does not yet support coefficient-mode plots. Use renderer='base'.",
         call. = FALSE)
  }

  if (!identical(as.character(view)[1L], "fixed")) {
    stop("renderer='rgl' currently supports view='fixed' only in this rollout tranche.",
         call. = FALSE)
  }

  if (!identical(as.character(plot.errors.method)[1L], "none")) {
    if (!isTRUE(allow.plot.errors)) {
      stop("renderer=\"rgl\" does not yet support errors != \"none\". Use renderer=\"base\".",
           call. = FALSE)
    }
  }

  if (isTRUE(plot.data.overlay)) {
    if (!isTRUE(allow.plot.data.overlay)) {
      stop("renderer='rgl' does not yet support plot.data.overlay=TRUE. Use renderer='base' or disable overlay.",
           call. = FALSE)
    }
  }

  if (!(plot.behavior %in% c("plot", "plot-data"))) {
    stop("renderer=\"rgl\" currently supports behavior %in% c(\"plot\", \"plot-data\", \"data\") only in this rollout tranche.",
         call. = FALSE)
  }

  if (!isTRUE(suppressWarnings(requireNamespace("rgl", quietly = TRUE)))) {
    stop("renderer='rgl' requires the suggested package 'rgl'. Please install it with install.packages('rgl').",
         call. = FALSE)
  }

  renderer
}

.np_plot_rgl_finalize <- function(rgl.out,
                                  plot.behavior = "plot",
                                  plot.data = NULL) {
  plot.behavior <- as.character(plot.behavior)[1L]

  if (identical(plot.behavior, "plot-data"))
    return(plot.data)

  if (!is.null(rgl.out))
    return(rgl.out)

  invisible(NULL)
}

.np_plot_render_surface_rgl <- function(x,
                                        y,
                                        z,
                                        xlab,
                                        ylab,
                                        zlab,
                                        main,
                                        theta,
                                        phi,
                                        col = NULL,
                                        border = .np_plot_color("surface_border"),
                                        zlim = NULL,
                                        par3d.args = list(),
                                        view3d.args = list(),
                                        persp3d.args = list(),
                                        grid3d.args = list(),
                                        widget.args = list(),
                                        draw.extras = NULL) {
  tryCatch({
    old.opts <- options(
      rgl.useNULL = TRUE,
      rgl.printRglwidget = TRUE
    )
    on.exit(options(old.opts), add = TRUE)
    old.env <- Sys.getenv("RGL_USE_NULL", unset = NA_character_)
    Sys.setenv(RGL_USE_NULL = "TRUE")
    on.exit({
      if (is.na(old.env)) {
        Sys.unsetenv("RGL_USE_NULL")
      } else {
        Sys.setenv(RGL_USE_NULL = old.env)
      }
    }, add = TRUE)
    devices.before <- try(rgl::rgl.dev.list(), silent = TRUE)
    if (inherits(devices.before, "try-error") || is.null(devices.before))
      devices.before <- integer(0L)

    rgl::open3d(useNULL = TRUE, silent = TRUE)
    par3d.call <- .np_plot_merge_override_args(
      list(windowRect = c(900, 100, 900 + 640, 100 + 640)),
      par3d.args
    )
    do.call(rgl::par3d, par3d.call)

    view3d.call <- .np_plot_merge_override_args(
      list(theta = theta, phi = phi, fov = 80),
      view3d.args
    )
    do.call(rgl::view3d, view3d.call)

    persp3d.call <- .np_plot_merge_override_args(list(
      x = x,
      y = y,
      z = z,
      zlim = zlim,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      ticktype = "detailed",
      border = border,
      color = .np_plot_rgl_surface_colors(z = z, col = col),
      alpha = 0.6,
      back = "lines",
      main = main
    ), persp3d.args)
    do.call(rgl::persp3d, persp3d.call)

    if (!is.null(draw.extras))
      draw.extras()

    grid.side <- c("x", "y+", "z")
    if (!is.null(grid3d.args$side)) {
      grid.side <- grid3d.args$side
      grid3d.args$side <- NULL
    }
    grid3d.call <- c(list(grid.side), grid3d.args)
    do.call(rgl::grid3d, grid3d.call)
    if (isTRUE(rgl::rgl.useNULL()) || isTRUE(getOption("rgl.printRglwidget"))) {
      scene <- rgl::scene3d()
      devices.after <- try(rgl::rgl.dev.list(), silent = TRUE)
      if (inherits(devices.after, "try-error") || is.null(devices.after))
        devices.after <- integer(0L)
      new.devices <- setdiff(devices.after, devices.before)
      if (length(new.devices))
        try(rgl::close3d(dev = new.devices, silent = TRUE), silent = TRUE)
      else
        try(rgl::close3d(silent = TRUE), silent = TRUE)
      widget <- do.call(rgl::rglwidget, c(list(x = scene), widget.args))
      print(widget)
      return(widget)
    }
    NULL
  }, error = function(e) {
    stop(sprintf("renderer='rgl' failed to draw the surface (%s)", conditionMessage(e)),
         call. = FALSE)
  })
}

.np_plot_collect_rgl_args <- function(dots, type, prefix) {
  direct.args <- .np_plot_user_args(dots, type)
  prefixed.args <- .np_plot_extract_prefixed_args(dots, prefix)
  .np_plot_merge_override_args(direct.args, prefixed.args)
}

.np_plot_rgl_controls <- function(dots, legend = TRUE) {
  legend3d <- .np_plot_collect_rgl_args(
    dots,
    "rgl.legend3d",
    "rgl.legend3d."
  )

  list(
    persp3d = .np_plot_collect_rgl_args(
      dots,
      "rgl.persp3d",
      "rgl.persp3d."
    ),
    view3d = .np_plot_collect_rgl_args(
      dots,
      "rgl.view3d",
      "rgl.view3d."
    ),
    par3d = .np_plot_collect_rgl_args(
      dots,
      "rgl.par3d",
      "rgl.par3d."
    ),
    grid3d = .np_plot_collect_rgl_args(
      dots,
      "rgl.grid3d",
      "rgl.grid3d."
    ),
    widget = .np_plot_collect_rgl_args(
      dots,
      "rgl.widget",
      "rgl.widget."
    ),
    legend3d = .np_plot_merge_rgl_legend_control(legend3d, legend),
    points3d = .np_plot_extract_prefixed_args(dots, "rgl.points3d."),
    surface3d = .np_plot_extract_prefixed_args(dots, "rgl.surface3d.")
  )
}

.np_plot_merge_rgl_legend_control <- function(legend3d.args, legend = TRUE) {
  legend.value <- if (is.list(legend) && any(names(legend) %in% c("tau", "bands"))) {
      if (!is.null(legend$bands)) legend$bands else TRUE
    } else {
      legend
    }
  if (is.null(legend.value))
    return(.np_plot_merge_override_args(legend3d.args, list(plot = FALSE)))
  if (is.logical(legend.value)) {
    if (length(legend.value) != 1L)
      stop("legend must be TRUE/FALSE, NULL, NA, a legend position string, or a list of graphics::legend arguments",
           call. = FALSE)
    if (is.na(legend.value) || !isTRUE(legend.value))
      return(.np_plot_merge_override_args(legend3d.args, list(plot = FALSE)))
    return(legend3d.args)
  }
  if (is.character(legend.value) && length(legend.value) == 1L && !is.na(legend.value))
    return(.np_plot_merge_override_args(list(x = legend.value), legend3d.args))
  if (is.list(legend.value)) {
    show <- legend.value$show
    if (!is.null(show)) {
      if (!is.logical(show) || length(show) != 1L)
        stop("legend$show must be TRUE or FALSE", call. = FALSE)
      legend.value$show <- NULL
      if (is.na(show) || !isTRUE(show))
        return(.np_plot_merge_override_args(legend3d.args, list(plot = FALSE)))
    }
    return(.np_plot_merge_override_args(legend.value, legend3d.args))
  }
  stop("legend must be TRUE/FALSE, NULL, NA, a legend position string, or a list of graphics::legend arguments",
       call. = FALSE)
}
