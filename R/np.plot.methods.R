.np_plot_from_slot <- function(object, slot = "bws", ...) {
  bws <- object[[slot]]
  if (is.null(bws))
    stop("plot object does not contain expected bandwidth slot")
  npplot(bws = bws, ...)
}

.np_plot_npregression <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_npdensity <- function(object, ...) {
  dots <- list(...)

  direct.args <- c("plot.behavior", "plot.errors.method", "plot.errors.type",
                   "plot.errors.boot.num", "plot.errors.boot.method",
                   "plot.errors.alpha", "perspective", "gradients",
                   "xdat", "data", "neval", "xtrim", "xq")
  use.direct <- isTRUE(object$ndim == 1) &&
    isTRUE(object$trainiseval) &&
    !any(direct.args %in% names(dots))

  if (use.direct) {
    ex <- object$eval[[1]]
    if (!is.factor(ex)) {
      ord <- order(ex)
      xlab.val <- if (!is.null(dots$xlab)) dots$xlab else gen.label(object$xnames[1], "X1")
      ylab.val <- if (!is.null(dots$ylab)) dots$ylab else "Density"
      type.val <- if (!is.null(dots$type)) dots$type else "l"

      dots$xlab <- NULL
      dots$ylab <- NULL
      dots$type <- NULL

      do.call(plot, c(list(x = as.numeric(ex[ord]),
                           y = as.numeric(object$dens[ord]),
                           type = type.val,
                           xlab = xlab.val,
                           ylab = ylab.val),
                      dots))
      return(invisible(object))
    }
  }

  .np_plot_from_slot(object, "bws", ...)
}
.np_plot_condensity <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_condistribution <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_npdistribution <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_qregression <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_singleindex <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_smoothcoefficient <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_plregression <- function(object, ...) .np_plot_from_slot(object, "bw", ...)

plot.bandwidth <- function(...) npplot(...)
plot.rbandwidth <- function(...) npplot(...)
plot.dbandwidth <- function(...) npplot(...)
plot.conbandwidth <- function(...) npplot(...)
plot.condbandwidth <- function(...) npplot(...)
plot.plbandwidth <- function(...) npplot(...)
plot.sibandwidth <- function(...) npplot(...)
plot.scbandwidth <- function(...) npplot(...)

plot.npregression <- function(x, ...) .np_plot_npregression(x, ...)
plot.npdensity <- function(x, ...) .np_plot_npdensity(x, ...)
plot.condensity <- function(x, ...) .np_plot_condensity(x, ...)
plot.condistribution <- function(x, ...) .np_plot_condistribution(x, ...)
plot.npdistribution <- function(x, ...) .np_plot_npdistribution(x, ...)
plot.qregression <- function(x, ...) .np_plot_qregression(x, ...)
plot.singleindex <- function(x, ...) .np_plot_singleindex(x, ...)
plot.smoothcoefficient <- function(x, ...) .np_plot_smoothcoefficient(x, ...)
plot.plregression <- function(x, ...) .np_plot_plregression(x, ...)
