# Private distribution-specific inference controls; no regression pilot state.
np_udist_bias_control <- function(construction = c("joint-process", "cdf", "legacy"),
                                   calibration = c("default", "score-envelope")) {
  construction <- match.arg(construction)
  calibration <- match.arg(calibration)
  if (identical(calibration, "default"))
    return(structure(list(construction = construction, version = 1L),
                     class = "np_udist_bias_control"))
  if (!identical(construction, "joint-process"))
    stop("score-envelope requires the joint-process CDF construction", call. = FALSE)
  structure(list(construction = construction, calibration = calibration, version = 2L),
            class = "np_udist_bias_control")
}

.np_udist_plot_prepare <- function(control, bws, x, center, method, bootstrap, band, blocklen = NULL, block.automatic = FALSE, calibration = NULL) {
  requested <- calibration
  calibration <- .np_plot_calibration_method(calibration, "distribution")
  if (identical(calibration, "default")) calibration <- NULL
  if (!is.null(calibration) && (!is.character(calibration) || length(calibration) != 1L ||
      is.na(calibration) || !calibration %in% c("absolute", "separate-tails", "tail-envelope", "smoothed-percentile") ||
      !identical(center, "bias-corrected") || !identical(method, "bootstrap") ||
      !bootstrap %in% c("inid", "fixed", "geom") ||
      !band %in% c("pointwise", "simultaneous", "bonferroni", "all")))
    stop("explicit CDF calibration requires corrected joint-process bootstrap empirical bands", call. = FALSE)
  if (is.null(control)) {
    if (!identical(center, "bias-corrected")) return(NULL)
    control <- np_udist_bias_control()
  }
  construction <- if (is.list(control)) control[["construction", exact = TRUE]] else NULL
  control.calibration <- if (is.list(control)) control[["calibration", exact = TRUE]] else NULL
  if (is.null(control.calibration)) control.calibration <- "default"
  if (!is.character(control.calibration) || length(control.calibration) != 1L ||
      is.na(control.calibration) || !control.calibration %in% c("default", "score-envelope"))
    stop("bias.control must be created by np_udist_bias_control()", call. = FALSE)
  if (!inherits(control, "np_udist_bias_control") || length(construction) != 1L ||
      !is.character(construction) || is.na(construction) ||
      !construction %in% c("joint-process", "cdf", "legacy") ||
      !identical(control, np_udist_bias_control(construction, control.calibration)))
    stop("bias.control must be created by np_udist_bias_control()", call. = FALSE)
  if (!identical(center, "bias-corrected") || !method %in% c("asymptotic", "bootstrap"))
    stop("distribution bias.control requires center='bias-corrected' and asymptotic or bootstrap errors", call. = FALSE)
  if (!is.null(calibration) && !identical(construction, "joint-process"))
    stop("explicit CDF calibration requires the joint-process construction", call. = FALSE)
  if (identical(control.calibration, "score-envelope") &&
      (!identical(method, "asymptotic") || !is.null(calibration) ||
       !band %in% c("pointwise", "simultaneous", "bonferroni", "all")))
    stop("CDF score-envelope requires asymptotic empirical bands without bootstrap calibration", call. = FALSE)
  if (identical(construction, "legacy")) {
    if (!identical(method, "bootstrap"))
      stop("legacy distribution bias correction requires errors='bootstrap'", call. = FALSE)
    return(NULL)
  }
  nearest <- !identical(bws[["type", exact = TRUE]], "fixed")
  if (nearest && !identical(construction, "joint-process"))
    stop("explicit cdf construction for nearest-neighbor bandwidths remains unqualified", call. = FALSE)
  if (nearest && identical(method, "bootstrap") && !bootstrap %in% c("inid", "fixed", "geom"))
    stop("selected nearest-neighbor CDF bootstrap law remains unqualified", call. = FALSE)
  if (identical(method, "bootstrap") &&
      !bootstrap %in% if (construction == "joint-process") c("inid", "fixed", "geom") else "inid")
    stop("selected bootstrap law is not qualified for this distribution construction", call. = FALSE)
  if (!band %in% c("pmzsd", "pointwise", "simultaneous", "bonferroni", "all"))
    stop("unsupported bias-corrected distribution band", call. = FALSE)
  if (identical(calibration, "smoothed-percentile") &&
      (ncol(x) != 1L || !is.numeric(x[[1L]]) || !identical(bws$type, "fixed")))
    stop("CDF smoothed-percentile requires a fixed univariate numeric fit", call.=FALSE)
  state <- .np_udist_bc_prepare(bws, x)
  if (identical(construction, "joint-process") &&
      (!state$q %in% 0:2 || !state$s %in% 0:2))
    stop("joint CDF process currently admits at most two numeric and two ordered variables; use an explicit qualified construction", call. = FALSE)
  if (!is.null(calibration) && calibration != "absolute" &&
      !identical(state$type, "fixed") && state$q > 0L)
    stop("asymmetric numeric nearest-neighbor CDF calibration remains unqualified", call. = FALSE)
  if (identical(control.calibration, "score-envelope")) {
    if (!identical(state$type, "fixed") || !state$q %in% 1:2 || state$s != 0L)
      stop("CDF score-envelope requires fixed numeric Gaussian2 geometry with one or two variables", call. = FALSE)
    calibration <- "score-envelope"
  }
  if (!is.null(calibration)) state$plot.calibration <- calibration
  if (!identical(requested, calibration)) state$plot.calibration.requested <- requested
  state$plot.construction <- construction
  state$plot.law <- if (identical(method, "asymptotic")) "gaussian" else bootstrap
  if (state$plot.law %in% c("fixed", "geom")) {
    state$plot.block.selection <- .np_udist_bc_block_request(x, state$plot.law, blocklen)
    if (isTRUE(block.automatic) && !is.null(blocklen)) {
      state$plot.block.selection$automatic <- TRUE
      state$plot.block.selection$source <- "b.star(xdat, round=TRUE)[1,1]"
      state$plot.block.selection$owner <- "shared plot normalization"
    }
  }
  if (identical(calibration, "smoothed-percentile") &&
      state$plot.law != "inid" && state$plot.block.selection$length != 1L)
    stop("CDF smoothed-percentile admits fixed/geom only with blocklen=1", call.=FALSE)
  state
}

.np_udist_plot_errors <- function(state, grid, alpha, B, band) {
  names(grid) <- names(state$x)
  provider.band <- if (band == "pmzsd") "pointwise" else band
  infer <- if (identical(state$plot.construction, "joint-process"))
    .np_udist_joint_inference else .np_udist_bc_inference
  inference <- .np_plot_activity_run("Constructing bias-corrected distribution intervals",
    if (identical(state[["plot.calibration", exact=TRUE]], "smoothed-percentile"))
      .np_udist_smoothed_plot(state, grid, alpha, B,
        if(provider.band=="all") c("pointwise","simultaneous","bonferroni") else provider.band) else
    do.call(infer, c(list(state = state, grid = grid, law = state$plot.law,
      alpha = alpha, B = B,
      bands = if (provider.band == "all") c("pointwise", "simultaneous", "bonferroni") else provider.band),
      if (identical(state$plot.construction, "joint-process"))
        list(selection = state$plot.block.selection,
             calibration = state[["plot.calibration", exact = TRUE]]) else list())))
  if (!is.null(state[["plot.calibration.requested", exact = TRUE]]))
    inference$settings$calibration.requested <- state$plot.calibration.requested
  inference$bands <- inference$intervals
  if (band == "pmzsd") inference <- .np_plot_bias_normal(inference, alpha)
  if (band != "pmzsd" && !identical(state[["plot.calibration", exact=TRUE]], "smoothed-percentile")) .np_band_resolution_warn(inference$resolution)
  all <- lapply(inference$bands, function(z)
    cbind(inference$center - z[, 1L, 1L], z[, 2L, 1L] - inference$center))
  width <- all[[if (band == "all") "pointwise" else band]]
  list(err = cbind(width, inference$center),
       all.err = if (band == "all") all else NULL,
       center = inference$center, inference = inference, bxp = list())
}
