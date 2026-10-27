# Private IID hybrid: empirical bootstrap-t pivots, smoothed observed scale.
# Internal calibration stage; public density plotting applies the F4 scale map.
.np_udens_b3_inid <- function(state, grid, alpha = c(.05, .1), B = 9999L,
                               tile = 25L, corrected = TRUE) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (existed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind(); success <- FALSE
  on.exit({
    if (!success) {
      do.call(RNGkind, as.list(kind))
      if (existed) assign(".Random.seed", previous, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!identical(state[["metadata", exact = TRUE]][["type", exact = TRUE]], "fixed"))
    stop("hybrid density inference requires a fixed-bandwidth state", call. = FALSE)
  .np_udens_bc_grid(state, grid)
  budget <- .np_udens_b2_budget(state, nrow(grid), B, tile, TRUE)
  model <- .np_udens_b2_smoothed_covariance(state, grid, tile, corrected)
  base <- .np_udens_b2_bootstrap(state, grid, law = "inid", blocklen = 1L,
                                alpha = alpha, B = B, tile = tile,
                                corrected = corrected)
  if (!identical(model$center, base$center))
    stop("hybrid density center disagrees with the empirical bootstrap", call. = FALSE)
  ratio <- model$se/base$se
  if (any(!is.finite(ratio)) || any(ratio <= 0))
    stop("hybrid density inference requires positive finite scale ratios", call. = FALSE)
  rescale <- function(intervals) {
    lapply(intervals, function(z) {
      out <- sweep(sweep(z, 1L, base$center, "-"), 1L, ratio, "*")
      out <- sweep(out, 1L, base$center, "+")
      if (anyNA(out) || any(out[, 1L, ] > out[, 2L, ]))
        stop("invalid hybrid density interval ordering", call. = FALSE)
      out
    })
  }
  base$intervals <- rescale(base$intervals)
  base$absolute.intervals <- rescale(base$absolute.intervals)
  base$empirical.se <- base$se
  base$se <- model$se
  base$covariance <- model$covariance
  base$settings$construction <- "density-bc-v3-hybrid-inid"
  base$settings$primary <- "absolute"
  base$settings$observed.scale <- "same-bandwidth positive-KDE covariance"
  base$settings$bootstrap.scale <- "empirical complete replicate studentizer"
  base$settings$budget <- budget
  success <- TRUE
  base
}
