# A display-only fallback permits replay of an entirely unavailable target.
# It does not change retained values, intervals or their statuses.
.np_plot_data_limits <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) range(x) else c(0, 1)
}

.np_plot_data_surface <- function(grid, values, columns) {
  x <- grid[[columns[1L]]]; y <- grid[[columns[2L]]]
  ux <- sort(unique(x)); uy <- sort(unique(y))
  index <- match(x, ux) + (match(y, uy) - 1L) * length(ux)
  if (length(unique(index)) != length(ux) * length(uy))
    stop("plot surface grid is not rectangular", call. = FALSE)
  for (i in which(duplicated(index))) {
    first <- match(index[i], index)
    if (!identical(values[i], values[first]))
      stop("duplicate surface coordinates have conflicting values", call. = FALSE)
  }
  z <- matrix(NA_real_, length(ux), length(uy)); z[index] <- values
  list(x = ux, y = uy, z = z)
}

# Pure replay: no estimator, bandwidth-selection, bootstrap or RNG operation is
# reachable from this method. Every endpoint comes from the interval table.
plot.np_plot_data <- function(x, y = NULL, ..., panel = NULL, target = NULL,
                               band = NULL, renderer = c("base", "rgl")) {
  .np_plot_check_public_names(.np_plot_dot_names(substitute(list(...))[-1L]),
                              replay = TRUE)
  if (!is.null(y)) stop("y is not used when replaying plot data", call. = FALSE)
  .np_plot_data_validate(x)
  renderer <- match.arg(renderer)
  dots <- list(...)
  if (x$family %in% c("npcopula", "npconmode"))
    .np_plot_reject_unconstructed_grid(dots, "saved copula/conditional-mode replay",
                                      "create new plot data to change the grid")
  panels <- if (is.null(panel)) names(x$panels) else .np_plot_data_select(panel, names(x$panels), "panel")
  jobs <- list()
  for (pid in panels) {
    ids <- unique(x$panels[[pid]]$values$target_id)
    if (!is.null(target)) ids <- .np_plot_data_select(target, ids, "target")
    for (tid in ids) jobs[[length(jobs) + 1L]] <- list(panel = pid, target = tid)
  }
  if (!length(jobs)) stop("plot data contains no targets to draw", call. = FALSE)
  if (length(jobs) > 1L && renderer == "base") {
    old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
    graphics::par(mfrow = grDevices::n2mfrow(length(jobs)))
  }
  widget <- NULL
  for (job in jobs) {
    p <- x$panels[[job$panel]]
    values <- p$values[p$values$target_id == job$target, , drop = FALSE]
    intervals <- p$intervals[p$intervals$target_id == job$target, , drop = FALSE]
    if (!is.null(band)) {
      selected <- .np_plot_data_select(band, unique(intervals$band), "band")
      intervals <- intervals[intervals$band == selected, , drop = FALSE]
    }
    grid <- p$grid[values$row_id, , drop = FALSE]
    label <- x$targets$label[match(job$target, x$targets$target_id)]
    display <- p$display
    style <- utils::modifyList(.np_plot_data_field(display, "style", list()), dots)
    if (identical(display$kind, "iv")) {
      spec <- list(curve.args = c(list(x = grid[[1L]], y = values$estimate), display$curve_args),
        frame.args = display$frame_args,
        continuous = display$continuous, overlay = !is.null(display$overlay),
        overlay.x = display$overlay$x, overlay.y = display$overlay$y,
        overlay.args = display$overlay_args, rug = !is.null(display$rug), rug.x = display$rug)
      .np_iv_plot_render(.np_iv_plot_replay_spec(spec, dots))
      next
    }
    varying <- if (!is.null(display$coordinates)) match(display$coordinates, names(grid)) else
      which(vapply(grid, function(g) length(unique(g)) > 1L, logical(1L)))
    if (identical(display$kind, "coefficients")) {
      if (renderer != "base") stop("coefficient displays use the base renderer", call. = FALSE)
      positions <- do.call(graphics::barplot, .np_plot_literal_labels(utils::modifyList(list(height = values$estimate,
        names.arg = grid[[1L]], ylab = label, xlab = "Regressor"),
        style[intersect(names(style), c("col", "border", "main", "sub", "xlab", "ylab", "ylim", "xlim"))])))
      if (nrow(intervals)) draw.error.bars(as.vector(positions), intervals$lower, intervals$upper)
      next
    }
    if (identical(display$kind, "scatter")) {
      if (renderer != "base") stop("probability-pair replay uses the base renderer", call. = FALSE)
      do.call(graphics::plot, .np_plot_literal_labels(utils::modifyList(list(x = grid[[1L]], y = grid[[2L]],
        xlab = p$variables$label[1L], ylab = p$variables$label[2L], main = label), style)))
      next
    }
    if (identical(display$kind, "surface") && length(varying) == 2L) {
      g <- .np_plot_data_surface(grid, values$estimate, varying)
      if (isTRUE(display$view %in% c("contour", "image"))) {
        if (renderer != "base") stop("contour/image replay uses the base renderer", call. = FALSE)
        fun <- if (display$view == "contour") graphics::contour else graphics::image
        do.call(fun, .np_plot_literal_labels(utils::modifyList(list(x = g$x, y = g$y, z = g$z,
          xlab = p$variables$label[varying[1L]], ylab = p$variables$label[varying[2L]], main = label), style)))
        next
      }
      lower <- upper <- list()
      for (bn in unique(intervals$band)) {
        d <- intervals[intervals$band == bn, , drop = FALSE]
        order <- match(values$row_id, d$row_id)
        lower[[bn]] <- .np_plot_data_surface(grid, d$lower[order], varying)$z
        upper[[bn]] <- .np_plot_data_surface(grid, d$upper[order], varying)$z
      }
      args <- utils::modifyList(list(x = as.numeric(g$x), y = as.numeric(g$y), z = g$z,
        xlab = p$variables$label[varying[1L]], ylab = p$variables$label[varying[2L]],
        zlab = label, main = label, theta = 0, phi = 20,
        zlim = .np_plot_expand_surface_range(.np_plot_data_limits(
          c(values$estimate, intervals$lower, intervals$upper, display$overlay$value)))),
        style[!names(style) %in% c("type", "lty", "lwd")])
      if (identical(x$family, "npcopula") && !is.null(style$zlim)) {
        # Match fresh copula display clipping, leaving all saved numbers intact.
        args$zlim <- .npcopula_validate_zlim(style$zlim)
        args$z <- .npcopula_clip_z_for_display(args$z, args$zlim)
        lower <- .npcopula_clip_z_list_for_display(lower, args$zlim)
        upper <- .npcopula_clip_z_list_for_display(upper, args$zlim)
      }
      all <- length(lower) > 1L
      if (renderer == "rgl") {
        args$draw.extras <- function() {
          if (!is.null(display$overlay)) do.call(.np_plot_overlay_points_rgl,
            c(list(x1 = as.numeric(display$overlay$grid[[1L]]),
              x2 = as.numeric(display$overlay$grid[[2L]]), y = display$overlay$value), display$points_args))
          if (!is.null(display$rug)) .np_plot_draw_floor_rug_rgl(
            x1 = as.numeric(display$rug[[1L]]), x2 = as.numeric(display$rug[[2L]]), zlim = args$zlim)
          if (length(lower)) .np_plot_error_surfaces_rgl(
            x = args$x, y = args$y, plot.errors.type = if (all) "all" else names(lower),
            plot.errors.method = x$request$errors,
            lerr = lower[[1L]], herr = upper[[1L]],
            lerr.all = if (all) lower else NULL, herr.all = if (all) upper else NULL)
        }
        widget <- do.call(.np_plot_render_surface_rgl, .np_plot_literal_labels(args))
      } else {
        mat <- .np_plot_render_surface_base_frame(
          args, range(args$x), range(args$y), args$zlim)
        if (!is.null(display$overlay)) do.call(.np_plot_overlay_points_persp,
          c(list(x1 = as.numeric(display$overlay$grid[[1L]]),
            x2 = as.numeric(display$overlay$grid[[2L]]), y = display$overlay$value,
            persp.mat = mat), display$points_args))
        if (!is.null(display$rug)) .np_plot_draw_floor_rug_persp(
          x1 = as.numeric(display$rug[[1L]]), x2 = as.numeric(display$rug[[2L]]),
          zlim = args$zlim, persp.mat = mat)
        if (length(lower)) .np_plot_draw_error_wireframes_persp(
          x = args$x, y = args$y, persp.mat = mat,
          plot.errors.type = if (all) "all" else names(lower),
          lerr = lower[[1L]], herr = upper[[1L]],
          lerr.all = if (all) lower else NULL, herr.all = if (all) upper else NULL)
      }
      next
    }
    if (renderer != "base") stop("one-dimensional replay uses the base renderer", call. = FALSE)
    column <- if (length(varying)) varying[1L] else 1L
    coord <- grid[[column]]
    # A specific Z-role override wins; otherwise preserve the generic replay
    # xlim override, or the retained Z-role limit when neither is supplied.
    if (identical(p$variables$role[column], "Z") &&
        ("zlim" %in% names(dots) || !("xlim" %in% names(dots))))
      style$xlim <- style$zlim
    if (identical(p$variables$role[column], "Z") &&
        ("zlab" %in% names(dots) || !("xlab" %in% names(dots))))
      style$xlab <- style$zlab
    categorical <- is.factor(coord) || is.character(coord)
    xpos <- if (categorical) match(as.character(coord), unique(as.character(coord))) else coord
    args <- utils::modifyList(list(x = xpos, y = values$estimate,
      type = if (categorical) "p" else "l", xlab = p$variables$label[column],
      ylab = label, main = label,
      ylim = .np_plot_data_limits(c(values$estimate, values$bias_corrected, intervals$lower, intervals$upper,
        display$overlay$value))),
      style[!names(style) %in% c("theta", "phi", "zlim", "zlab")])
    if (categorical) args$xaxt <- "n"
    do.call(graphics::plot, .np_plot_literal_labels(args))
    if (categorical) graphics::axis(1, at = unique(xpos), labels = unique(as.character(coord)))
    if (!is.null(display$overlay)) {
      fun <- if (categorical) .np_plot_overlay_points_factor else .np_plot_overlay_points_1d
      do.call(fun, c(list(x = display$overlay$grid[[1L]], y = display$overlay$value), display$points_args))
    }
    if (!is.null(display$rug) && !categorical) .np_plot_draw_rug_1d(display$rug[[1L]])
    box <- display$boxplots[[job$target]]
    if (categorical && isTRUE(display$boxplot_enabled) && !is.null(box)) {
      graphics::bxp(box, add = TRUE, outline = isTRUE(display$interval_style$boxplot_outliers))
      next
    }
    if (any(is.finite(values$bias_corrected))) graphics::lines(xpos, values$bias_corrected, lty = 3)
    colors <- .np_plot_all_band_colors()
    for (bn in unique(intervals$band)) {
      d <- intervals[intervals$band == bn, , drop = FALSE]
      ord <- match(values$row_id, d$row_id)
      col <- if (length(unique(intervals$band)) > 1L && bn %in% names(colors)) colors[[bn]] else "black"
      if (categorical || identical(display$interval_style$style[1L], "bar")) {
        keep <- which(is.finite(xpos) & is.finite(d$lower[ord]) & is.finite(d$upper[ord]))
        count <- display$interval_style$bar_num
        if (!categorical && length(count) && length(keep) > count)
          keep <- keep[unique(round(seq(1, length(keep), length.out = count)))]
        if (length(keep)) graphics::arrows(xpos[keep], d$lower[ord][keep], xpos[keep], d$upper[ord][keep],
          angle = 90, code = 3, length = if (identical(display$interval_style$bar[1L], "I")) .03 else 0, col = col)
      } else {
        geometry <- .np_plot_interval_geometry(xpos, d$lower[ord], d$upper[ord])
        graphics::lines(geometry$x, geometry$lower, col = col, lty = 2)
        graphics::lines(geometry$x, geometry$upper, col = col, lty = 2)
      }
    }
  }
  if (renderer == "rgl") widget else invisible(NULL)
}
