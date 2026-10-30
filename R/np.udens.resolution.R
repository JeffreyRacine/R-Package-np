# Location and scale resolution are diagnostics, not arithmetic validity gates.
.np_udens_resolution <- function(B, m, alpha, bands, calibration = "absolute") {
  signed <- calibration != "absolute"
  z <- .np_band_resolution(B, rep("density", m), alpha, rep(TRUE, m), bands,
    signed = signed, simultaneous.index = if (signed) matrix(1 + floor(B*alpha/2), 1L) else NULL)
  z$calibration <- calibration
  z
}

.np_udens_resolution_scale <- function(B, m, tau) {
  z <- .np_band_resolution(B, rep("density", m), .05, tau > 0,
                           c("pointwise", "simultaneous"))
  z$band <- paste0("scale-", z$band)
  z$calibration <- "one-sided log-scale"
  z
}

# Aggregate across curve slices only after the complete public call succeeds.
# Error paths retain their original condition and do not emit partial results.
.np_udens_resolution_run <- function(code) {
  diagnostics <- list()
  percentile.messages <- character()
  value <- withVisible(withCallingHandlers(force(code),
    np_calibrated_percentile_warning = function(w) {
      percentile.messages <<- c(percentile.messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    np_band_resolution_warning = function(w) {
      diagnostics[[length(diagnostics) + 1L]] <<- w$diagnostics
      invokeRestart("muffleWarning")
    }))
  if (length(percentile.messages)) {
    rows <- if (length(diagnostics)) do.call(rbind, diagnostics) else NULL
    warning(structure(list(message = paste(unique(percentile.messages), collapse = "; "),
      call = NULL, diagnostics = rows),
      class = c("np_calibrated_percentile_warning", "warning", "condition")))
  } else .np_band_resolution_warn(diagnostics)
  if (value$visible) value$value else invisible(value$value)
}
