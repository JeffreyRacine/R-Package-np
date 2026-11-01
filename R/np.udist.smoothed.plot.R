# CDF public bridge to the unchanged CDS1 statistical construction.
.np_udist_smoothed_plot <- function(state, grid, alpha, B, bands) {
  selection <- .np_bw_outer_selection(state$bws, 'distribution')
  streams <- .np_udens_calibrated_streams(get('.Random.seed', .GlobalEnv, inherits = FALSE))
  blocklen <- if (state$plot.law == 'inid') NULL else state$plot.block.selection$length
  value <- .np_udist_smoothed_percentile(list(bws = state$bws), state$x, grid,
    streams$original, streams$outer,
    outer.fit = selection$fit,
    B = B, alpha = alpha, law = state$plot.law, blocklen = blocklen,
    max.bytes = state$max.bytes, warn = FALSE, prepared = state,
    progress = function(done, total) .np_plot_progress_notify(
      stage = 'Calibrating CDF percentile outer samples', done = done, total = total, unit = 'samples'))
  center <- original <- se <- numeric(nrow(grid))
  for (first in seq.int(1L, nrow(grid), 25L)) {
    ii <- first:min(nrow(grid), first + 24L)
    z <- .np_udist_bc_rows(state, grid[ii, , drop = FALSE])
    center[ii] <- z$center; original[ii] <- z$original.center; se[ii] <- z$se
  }
  intervals <- lapply(value$inference$bands[bands], function(z)
    array(z, c(nrow(grid), 2L, 1L), list(NULL, c('lower', 'upper'), as.character(alpha))))
  resolution <- .np_band_resolution(B, rep('cdf', nrow(grid)), alpha, rep(TRUE, nrow(grid)),
    bands, signed = TRUE, simultaneous.index = matrix(value$inference$order.index[['simultaneous']], 1L))
  resolution <- .np_band_resolution_use_indices(resolution,
    as.list(value$inference$order.index[resolution$band]))
  messages <- unique(c(selection$messages, value$warnings))
  withCallingHandlers(.np_band_resolution_warn(resolution), warning = function(w) {
    messages <<- c(messages, conditionMessage(w)); invokeRestart('muffleWarning')
  })
  if (length(messages)) warning(structure(list(message = paste(unique(messages), collapse = '; '),
    call = NULL, diagnostics = resolution), class = c('np_calibrated_percentile_warning', 'warning', 'condition')))
  value$settings$bandwidth.procedure <- selection$diagnostics
  value$settings$calibration.requested <- 'smoothed-percentile'
  value$settings$scope <- 'experimental q1 fixed Gaussian2; useful population coverage unassessed'
  list(center = center, original.center = original, se = se, covariance = NULL,
    bands = intervals, intervals = intervals, calibration = value$inference,
    resolution = resolution, settings = value$settings, pilot = value$pilot,
    rng = value$rng, calibration.streams = streams, outer.diagnostics = value$outer,
    warnings = unique(messages), budget = value$budget)
}
