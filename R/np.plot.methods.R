.np_plot_from_slot <- function(object, slot = "bws", ...) {
  bws <- object[[slot]]
  if (is.null(bws))
    stop("plot object does not contain expected bandwidth slot")
  npplot(bws = bws, ...)
}

.np_plot_npregression <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_condensity <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_condistribution <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_npdistribution <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_qregression <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_singleindex <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_smoothcoefficient <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_plregression <- function(object, ...) .np_plot_from_slot(object, "bw", ...)

.np_plot_bandwidth <- function(...) npplot(...)
.np_plot_rbandwidth <- function(...) npplot(...)
.np_plot_dbandwidth <- function(...) npplot(...)
.np_plot_conbandwidth <- function(...) npplot(...)
.np_plot_condbandwidth <- function(...) npplot(...)
.np_plot_plbandwidth <- function(...) npplot(...)
.np_plot_sibandwidth <- function(...) npplot(...)
.np_plot_scbandwidth <- function(...) npplot(...)
