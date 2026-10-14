# Interval geometry for the explicitly isolated density/distribution/conmode
# renderers. Legacy shared renderers remain unchanged for protected consumers.
.np_plot_interval_geometry <- function(ex, ely, ehy) {
  n <- min(length(ex), length(ely), length(ehy))
  idx <- seq_len(n)
  x <- as.double(ex[idx])
  lower <- as.double(ely[idx])
  upper <- as.double(ehy[idx])
  finite <- is.finite(x) & is.finite(lower) & is.finite(upper)
  # Keep the declared row positions: NA separates finite line runs.
  x[!finite] <- lower[!finite] <- upper[!finite] <- NA_real_
  list(x = x, lower = lower, upper = upper, finite = finite)
}

.np_plot_error_vectors_with_gaps <- function(x, estimate, err,
                                             plotOnEstimate = TRUE) {
  bounds <- .np_plot_error_bounds(estimate, err, plotOnEstimate)
  .np_plot_interval_geometry(x, bounds$lower, bounds$upper)
}

.np_plot_draw_errors_with_gaps <- function(ex, ely, ehy,
                                           plot.errors.style,
                                           plot.errors.bar,
                                           plot.errors.bar.num,
                                           lty,
                                           col = par("col"),
                                           annotation = NULL) {
  if (plot.errors.style == "bar")
    return(draw.errors(ex, ely, ehy, plot.errors.style, plot.errors.bar,
                       plot.errors.bar.num, lty, col, annotation))
  geometry <- .np_plot_interval_geometry(ex, ely, ehy)
  if (!any(geometry$finite))
    return(invisible(NULL))
  if (plot.errors.style == "band")
    draw.error.bands(ex = geometry$x, ely = geometry$lower,
                     ehy = geometry$upper, lty = lty, col = col)
  .np_plot_draw_variability_annotation(annotation)
}

.np_plot_draw_all_errors_with_gaps <- function(ex, center, all.err,
                                               plot.errors.style = "band",
                                               plot.errors.bar = "|",
                                               plot.errors.bar.num = min(length(ex), 25),
                                               lty = .np_plot_lty("interval"),
                                               add.legend = TRUE,
                                               legend.loc = "topleft",
                                               legend = TRUE,
                                               xi.factor = FALSE,
                                               plot.errors.method,
                                               annotation = NULL) {
  if (is.null(all.err)) return(invisible(NULL))
  if (xi.factor) {
    plot.errors.style <- "bar"
    plot.errors.bar <- "I"
  }
  draw_one <- function(err, col) {
    if (is.null(err)) return(invisible(FALSE))
    panel.center <- .np_plot_range_center(center, list(err))
    lower <- panel.center - err[, 1L]
    upper <- panel.center + err[, 2L]
    geometry <- .np_plot_interval_geometry(ex, lower, upper)
    if (!any(geometry$finite)) return(invisible(FALSE))
    .np_plot_draw_errors_with_gaps(
      ex = geometry$x, ely = geometry$lower, ehy = geometry$upper,
      plot.errors.style = plot.errors.style,
      plot.errors.bar = plot.errors.bar,
      plot.errors.bar.num = plot.errors.bar.num, lty = lty, col = col)
    invisible(TRUE)
  }
  band.cols <- .np_plot_all_band_colors()
  drawn <- c(
    pointwise = draw_one(all.err$pointwise, band.cols[["pointwise"]]),
    simultaneous = draw_one(all.err$simultaneous, band.cols[["simultaneous"]]),
    bonferroni = draw_one(all.err$bonferroni, band.cols[["bonferroni"]])
  )
  if (add.legend)
    .np_plot_draw_all_band_legend(
      legend = legend, x = legend.loc, lty = .np_plot_lty("interval"),
      lwd = .np_plot_lwd("band_all_1d"),
      plot.errors.method = plot.errors.method, drawn = drawn)
  if (any(drawn))
    .np_plot_draw_variability_annotation(annotation)
  invisible(NULL)
}
