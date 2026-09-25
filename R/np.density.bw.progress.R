.np_density_bw_progress <- function(bandwidth.compute, bwsolver, eval.only, expr) {
  if (isTRUE(bandwidth.compute) && !isTRUE(eval.only) &&
      is.character(bwsolver) && length(bwsolver) > 0L &&
      bwsolver[[1L]] %in% c("mads", "mads+powell") &&
      !.np_progress_bandwidth_active()) {
    .np_progress_select_bandwidth_enhanced(
      label = .np_progress_bandwidth_title(), expr = expr
    )
  } else {
    force(expr)
  }
}
