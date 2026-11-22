# The public boundary selects an existing provider. Keep literal density
# calibration names until the kernel-specific selector has inspected them.
.np_plot_interval_families <- c(
  plot.npdensity = "npudens", plot.bandwidth = "npudens",
  plot.npdistribution = "npudist", plot.dbandwidth = "npudist",
  plot.condensity = "npcdens", plot.conbandwidth = "npcdens",
  plot.condistribution = "npcdist", plot.condbandwidth = "npcdist",
  plot.plregression = "npplreg", plot.plbandwidth = "npplreg",
  plot.singleindex = "npindex", plot.sibandwidth = "npindex",
  plot.smoothcoefficient = "npscoef", plot.scbandwidth = "npscoef",
  plot.qregression = "npqreg", plot.lsqregression = "nplsqreg",
  plot.conmode = "npconmode", plot.npcopula = "npcopula"
)

.np_plot_family_interval_request <- function(dots, family) {
  supplied <- as.character(intersect(names(dots), c("interval", "band", "boot.control", "random.seed")))
  interval <- dots[["interval", exact = TRUE]]
  band <- dots[["band", exact = TRUE]]
  control <- dots[["boot.control", exact = TRUE]]
  if ("plot.errors.boot.calibration" %in% names(dots))
    stop("replace plot.errors.boot.calibration with the plot interval argument", call. = FALSE)
  if (inherits(control, "np_boot_control") && "calibration" %in% names(control)) {
    old <- .np_plot_boot_calibration(control)
    action <- switch(old,
      default = "omit interval to retain the default",
      legacy = "select construction='legacy' in the family bias.control",
      `tail-envelope` = "use interval='asymmetric-sd-envelope'",
      `smoothed-percentile` = "use interval='outer-calibrated-percentile'",
      sprintf("select the corresponding interval (previous calibration '%s')", old))
    stop(sprintf("rebuild boot.control without its calibration field and %s", action), call. = FALSE)
  }
  if (is.character(band) && length(band) == 1L && !is.na(band) && band == "pmzsd")
    stop("use interval='normal-sd', band='pointwise' in place of band='pmzsd'", call. = FALSE)
  if (!is.null(interval) && is.character(interval)) interval <- as.character(interval)
  dots$interval <- NULL
  dots$.np.interval.request <- list(interval = interval, supplied = supplied)
  if (is.null(interval)) return(dots)
  choices <- c("normal-sd", "percentile", "symmetric-sd", "asymmetric-sd",
               "asymmetric-sd-envelope", "outer-calibrated-percentile")
  if (family == "npudist") choices <- c(choices, "score-envelope")
  if (family == "npudens")
    choices <- c(choices, "absolute", "separate-tails", "percentile-tails")
  if (!is.character(interval) || length(interval) != 1L || is.na(interval) ||
      !interval %in% choices)
    stop(sprintf("unsupported interval construction for %s", family), call. = FALSE)
  method <- .np_plot_scalar_match(.np_plot_data_field(dots, "errors", "none"),
                                  c("none", "asymptotic", "bootstrap"), "errors")
  if (method == "none")
    stop("an explicit interval requires bootstrap or asymptotic errors", call. = FALSE)
  scope <- if (!"band" %in% supplied) "pointwise" else
    .np_plot_scalar_match(band, c("pointwise", "simultaneous", "bonferroni", "all"), "band")
  if (!"band" %in% supplied && isTRUE(dots[["bootstrap", exact = TRUE]] %in%
      c("smoothed-q6", "smoothed-beta-binomial")))
    stop("this smoothed density bootstrap requires an explicit band scope", call. = FALSE)
  if (is.null(scope)) stop("an explicit interval requires a valid band scope", call. = FALSE)
  corrected <- identical(.np_plot_scalar_match(.np_plot_data_field(dots, "center", "estimate"),
                                                c("estimate", "bias-corrected"), "center"), "bias-corrected")
  bias <- dots[["bias.control", exact = TRUE]]
  if (!is.null(bias) && !is.list(bias))
    stop("bias.control must be created by the corresponding family constructor", call. = FALSE)
  specialized <- corrected && (family %in% c("npudens", "npudist") ||
                    inherits(bias, "np_conditional_bias_control")) &&
    !identical(bias[["construction", exact = TRUE]], "legacy")
  calibration <- NULL
  if (interval == "normal-sd") {
    available <- scope == "pointwise" || (method == "asymptotic" && !specialized && scope == "bonferroni")
    if (scope == "pointwise") scope <- "pmzsd"
  } else if (interval == "percentile") {
    available <- method == "bootstrap" && !specialized
  } else if (interval == "score-envelope") {
    available <- corrected && method == "asymptotic" &&
      (is.null(bias) || identical(bias, np_udist_bias_control()) ||
         identical(bias, np_udist_bias_control(calibration = "score-envelope")))
    if (available) dots$bias.control <- np_udist_bias_control(calibration = "score-envelope")
  } else if (method == "asymptotic") {
    available <- interval == "symmetric-sd" && (specialized || scope == "simultaneous")
    if (identical(bias[["calibration", exact = TRUE]], "score-envelope")) available <- FALSE
  } else {
    available <- specialized
    calibration <- switch(interval, `asymmetric-sd-envelope` = "tail-envelope",
                          `outer-calibrated-percentile` = "smoothed-percentile", interval)
  }
  if (!available)
    stop(sprintf("interval='%s' is unavailable for this %s center, errors and band", interval, family),
         call. = FALSE)
  dots$band <- scope
  if (!is.null(calibration)) dots[["plot.errors.boot.calibration"]] <- calibration
  dots
}

# Publish components alongside provider-specific diagnostics. Historical result
# readers operate on the original objects and never enter this construction path.
# Use the emitted provider identity and the band actually returned. The literal
# selector can have two different owners across density kernels.
.np_plot_family_components <- function(meta, intervals, request, family) {
  parts <- list()
  actual <- meta$calibration
  if (is.null(actual)) actual <- meta$primary
  construction <- if (is.null(meta$construction)) "" else meta$construction
  for (key in unique(intervals$band)) {
    public <- if (key == "pmzsd") "pointwise" else key
    states <- unique(intervals$status[intervals$band == key])
    normal <- key == "pmzsd" || (identical(request$errors, "asymptotic") &&
      key %in% c("pointwise", "bonferroni") && !nzchar(construction) &&
      is.null(meta$correction))
    method <- if (identical(states, "undefined")) {
      NA_character_
    } else if (normal) {
      "normal-sd"
    } else if (identical(actual, "score-envelope")) {
      "score-envelope"
    } else if (construction %in% c("density-calibrated-unshifted-smoothed-percentile-v2",
                                    "cdf-calibrated-unshifted-smoothed-percentile-v1")) {
      "outer-calibrated-percentile"
    } else if (identical(request$errors, "asymptotic")) {
      if (identical(actual, "normal-quantile") ||
          (identical(actual, "Gaussian maximum and normal quantiles") && key != "simultaneous"))
        "normal-sd" else "symmetric-sd"
    } else if (identical(actual, "tail-envelope")) {
      "asymmetric-sd-envelope"
    } else if (isTRUE(actual %in% c("absolute", "symmetric-sd"))) {
      "symmetric-sd"
    } else if (isTRUE(actual %in% c("separate-tails", "percentile-tails", "asymmetric-sd"))) {
      if (identical(actual, "percentile-tails") ||
          (family == "npudens" && identical(actual, "separate-tails"))) actual else "asymmetric-sd"
    } else if (construction %in% c("cdf-absolute-process-v1", "cdf-bc-v1")) {
      "symmetric-sd"
    } else if (construction == "density-bc-v2-signed") {
      "asymmetric-sd"
    } else if (!nzchar(construction)) {
      "percentile"
    } else {
      # Some historical controls identify their construction only by provider.
      # Retain that identity rather than inferring a named selector.
      construction
    }
    resolution <- meta$resolution
    if (is.data.frame(resolution) && "band" %in% names(resolution))
      resolution <- resolution[resolution$band == key, , drop = FALSE]
    endpoint <- meta$endpoint_construction
    calibration <- actual
    studentization <- meta$studentization
    if (normal && !is.na(method)) {
      calibration <- "normal critical value"
      endpoint <- "center plus or minus a normal critical value times the marginal scale"
    } else if (is.null(endpoint) && !is.na(method)) {
      endpoint <- switch(method,
        percentile = "empirical bootstrap curve quantiles",
        `symmetric-sd` = "symmetric absolute-process bounds",
        `separate-tails` =, `asymmetric-sd` = "signed process limits",
        `asymmetric-sd-envelope` = "hull of signed and absolute process bounds",
        `score-envelope` = "process bounds with marginal score protection",
        NULL)
    }
    if (is.null(studentization) && identical(meta$critical_units, "CDF"))
      studentization <- "absolute CDF displacement without a studentizing divisor"
    parts[[public]] <- list(interval = method, calibration = calibration, status = states,
      endpoint_construction = endpoint, studentization = studentization, resolution = resolution,
      critical_value = if (key == "simultaneous") meta$critical_value else NULL)
  }
  meta$components <- parts
  meta
}

.np_plot_family_request_record <- function(request, selection) {
  if (identical(request$band, "pmzsd")) request$band <- "pointwise"
  request["interval"] <- list(selection$interval)
  request$supplied <- union(setdiff(request$supplied, c("band", "interval", "boot.control", "random.seed")), selection$supplied)
  request
}
