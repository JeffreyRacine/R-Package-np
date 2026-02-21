.np_plot_call_method <- function(method, bws, ..., where = "npplot()") {
  dots <- list(...)
  random.seed <- if (!is.null(dots$random.seed)) dots$random.seed else 42L
  dots$random.seed <- NULL

  # Keep backward compatibility with legacy plot argument aliases.
  if (!is.null(dots$gradient) && is.null(dots$gradients)) {
    dots$gradients <- dots$gradient
  }
  dots$gradient <- NULL

  if (!is.null(dots$persp) && is.null(dots$perspective)) {
    dots$perspective <- dots$persp
  }
  dots$persp <- NULL

  .npRmpi_require_active_slave_pool(where = where)
  .npRmpi_guard_no_auto_object_in_manual_bcast(bws, where = where)

  .np_with_seed(
    random.seed,
    do.call(method, c(list(bws = .npRmpi_autodispatch_untag(bws)), dots))
  )
}

.np_plot_from_slot <- function(object, slot = "bws", ...) {
  bws <- object[[slot]]
  if (is.null(bws))
    stop("plot object does not contain expected bandwidth slot")
  .np_plot_call_method(npplot, bws = bws, ..., where = "npplot()")
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

.np_plot_bandwidth <- function(bws, ...) {
  .np_plot_call_method(npplot.bandwidth, bws = bws, ..., where = "plot.bandwidth()")
}
.np_plot_rbandwidth <- function(bws, ...) {
  .np_plot_call_method(npplot.rbandwidth, bws = bws, ..., where = "plot.rbandwidth()")
}
.np_plot_dbandwidth <- function(bws, ...) {
  .np_plot_call_method(npplot.dbandwidth, bws = bws, ..., where = "plot.dbandwidth()")
}
.np_plot_conbandwidth <- function(bws, ...) {
  .np_plot_call_method(npplot.conbandwidth, bws = bws, ..., where = "plot.conbandwidth()")
}
.np_plot_condbandwidth <- function(bws, ...) {
  .np_plot_call_method(npplot.condbandwidth, bws = bws, ..., where = "plot.condbandwidth()")
}
.np_plot_plbandwidth <- function(bws, ...) {
  .np_plot_call_method(npplot.plbandwidth, bws = bws, ..., where = "plot.plbandwidth()")
}
.np_plot_sibandwidth <- function(bws, ...) {
  .np_plot_call_method(npplot.sibandwidth, bws = bws, ..., where = "plot.sibandwidth()")
}
.np_plot_scbandwidth <- function(bws, ...) {
  .np_plot_call_method(npplot.scbandwidth, bws = bws, ..., where = "plot.scbandwidth()")
}

plot.bandwidth <- function(x, ...) .np_plot_bandwidth(x, ...)
plot.rbandwidth <- function(x, ...) .np_plot_rbandwidth(x, ...)
plot.dbandwidth <- function(x, ...) .np_plot_dbandwidth(x, ...)
plot.conbandwidth <- function(x, ...) .np_plot_conbandwidth(x, ...)
plot.condbandwidth <- function(x, ...) .np_plot_condbandwidth(x, ...)
plot.plbandwidth <- function(x, ...) .np_plot_plbandwidth(x, ...)
plot.sibandwidth <- function(x, ...) .np_plot_sibandwidth(x, ...)
plot.scbandwidth <- function(x, ...) .np_plot_scbandwidth(x, ...)

plot.npregression <- function(x, ...) .np_plot_npregression(x, ...)
plot.npdensity <- function(x, ...) .np_plot_npdensity(x, ...)
plot.condensity <- function(x, ...) .np_plot_condensity(x, ...)
plot.condistribution <- function(x, ...) .np_plot_condistribution(x, ...)
plot.npdistribution <- function(x, ...) .np_plot_npdistribution(x, ...)
plot.qregression <- function(x, ...) .np_plot_qregression(x, ...)
plot.singleindex <- function(x, ...) .np_plot_singleindex(x, ...)
plot.smoothcoefficient <- function(x, ...) .np_plot_smoothcoefficient(x, ...)
plot.plregression <- function(x, ...) .np_plot_plregression(x, ...)
