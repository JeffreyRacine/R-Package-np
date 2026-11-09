# Public bridge to the calibrated unshifted corrected-curve owner.
# Original fit state is immutable. Fresh outer-sample CV is explicitly authorized.
.np_udens_calibrated_streams <- function(start) {
  seed <- 104729
  for (v in as.double(start)) seed <- (65599 * seed + v %% 2147483646) %% 2147483646
  .np_bias_stream(start, {
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(as.integer(seed + 1))
    current <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
    next.state <- function() {
      value <- current
      current <<- parallel::nextRNGStream(current)
      value
    }
    list(original = list(count = next.state(), jitter = next.state()),
      outer = lapply(seq_len(39L), function(i)
        list(data = next.state(), cv = next.state(), count = next.state(), jitter = next.state())))
  })$value
}

.np_udens_calibrated_family <- function(state, grid, design) {
  if (is.null(design)) stop("smoothed-percentile requires an explicit plot family design", call. = FALSE)
  # State numeric masks are absent on the pure-continuous legacy owner.
  nm <- if (is.null(state$numeric)) rep(TRUE, ncol(state$x)) else state$numeric
  vary <- design$vary
  if (!is.numeric(vary) || !length(vary) || anyNA(vary) ||
      any(vary != as.integer(vary)) || any(vary < 1L | vary > ncol(grid)))
    stop("invalid percentile plot family", call. = FALSE)
  points <- vapply(which(nm), function(j) length(unique(grid[[j]])), integer(1L))
  # Exact original owner for full sample-range Cartesian numeric families.
  full.range <- if (is.list(design$xtrim))
    all(vapply(design$xtrim[which(nm)], function(x)
      identical(as.numeric(x), c(0, 1)), logical(1L))) else
    all(design$xtrim[which(nm)] == 0)
  if (all(which(nm) %in% vary) && full.range)
    return(list(points = points, builder = NULL))
  prototype <- grid
  indices <- lapply(seq_len(ncol(grid)), function(j) match(grid[[j]], unique(grid[[j]])))
  count <- vapply(grid, function(x) length(unique(x)), integer(1L))
  builder <- function(sy) {
    out <- prototype
    for (j in which(nm)) {
      if (j %in% vary) {
        limits <- trim.quantiles(sy$x[[j]], design$xtrim[j])
        values <- seq(limits[1L], limits[2L], length.out = count[j])
        out[[j]] <- values[indices[[j]]]
      } else out[[j]] <- rep(uocquantile(sy$x[[j]], prob = design$xq[j]), nrow(out))
    }
    # Keep requested physical profiles, expressed in the outer fit's registry.
    # Ordered outer draws may add lattice levels beyond those in the original.
    coordinates <- toMatrix(prototype)
    for (category in sy$category) {
      j <- category$index
      at <- match(coordinates[, j], category$support)
      if (anyNA(at)) stop("percentile profile absent from outer category registry", call. = FALSE)
      out[[j]] <- factor(category$levels[at], levels = category$levels,
                          ordered = category$ordered)
    }
    out
  }
  list(points = points, builder = builder)
}

.np_udens_calibrated_plot <- function(state, grid, alpha, B, bands, design) {
  if (!requireNamespace("ks", quietly = TRUE))
    stop("smoothed-percentile requires the suggested package ks", call. = FALSE)
  start <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  streams <- .np_udens_calibrated_streams(start)
  family <- .np_udens_calibrated_family(state, grid, design)
  progress <- function(done, total)
    .np_plot_progress_notify(stage = "Calibrating smoothed-percentile outer samples",
      done = done, total = total, unit = "samples")
  value <- .np_udens_calibrated_percentile(list(bws = state$plot.bws), state$x, grid,
    streams$original, streams$outer, B = B, points = family$points,
    law = state$plot.law, blocklen = state[["plot.blocklen", exact = TRUE]],
    max.bytes = state$max.bytes, alpha = alpha, warn = FALSE, progress = progress,
    family.builder = family$builder, prepared = state)
  reduced <- value$inference
  # Descriptive center/SD are returned for the standard plot payload only;
  # neither studentizes, shifts nor alters any percentile endpoint.
  center <- original <- se <- numeric(nrow(grid))
  for (first in seq.int(1L, nrow(grid), by = 25L)) {
    ii <- first:min(nrow(grid), first + 24L)
    row <- .np_udens_bc_contribution(state, grid[ii, , drop = FALSE])
    center[ii] <- row$center; original[ii] <- row$original.center; se[ii] <- row$se
  }
  intervals <- lapply(reduced$bands[bands], function(z)
    array(z, c(nrow(grid), 2L, 1L), list(NULL, c("lower", "upper"), as.character(alpha))))
  resolution <- .np_band_resolution(B, rep("density", nrow(grid)), alpha,
    rep(TRUE, nrow(grid)), bands, signed = TRUE,
    simultaneous.index = matrix(reduced$order.index[["simultaneous"]], 1L))
  resolution <- .np_band_resolution_use_indices(resolution,
    as.list(reduced$order.index[resolution$band]))
  messages <- unique(c(reduced$messages, value$pilot.warnings,
    unlist(lapply(value$outer, function(z) z$warnings), use.names = FALSE)))
  if (length(messages)) warning(structure(list(message = paste(messages, collapse = "; "),
    call = NULL, diagnostics = resolution),
    class = c("np_calibrated_percentile_warning", "warning", "condition")))
  settings <- value$settings
  settings$construction <- "density-calibrated-unshifted-smoothed-percentile-v2"
  settings$calibration.requested <- "smoothed-percentile"
  settings$calibration.method <- "OP1+SC1 with MX1/NPA1/OLR1"
  settings$alpha <- alpha; settings$budget <- value$budget
  settings$family <- design; settings$available.bands <- names(intervals)
  settings$coverage.status <- "Gaussian q1 paired evidence; multivariate/mixed capability proof and stopped bounded assessments; general coverage unqualified"
  settings$bootstrap.statistic <- "complete unshifted bias-corrected smoothed curves"
  settings$endpoint.construction <- "outer-calibrated empirical ranks with companion hulls and certified full-support completion"
  settings$block.selection <- if (state$plot.law %in% c("fixed", "geom"))
    list(blocklen = state$plot.blocklen) else NULL
  list(center = center, original.center = original, se = se,
    empirical.se = se, covariance = NULL, bands = intervals, intervals = intervals,
    calibration = reduced, resolution = resolution, settings = settings,
    pilot = value$pilot, rng = value$rng, calibration.streams = streams,
    outer.diagnostics = value$outer, warnings = messages)
}
