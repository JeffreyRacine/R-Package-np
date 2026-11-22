# Translate only the two regression plot entries. The existing engine owns
# law selection, kernel eligibility, draws, finite-B behavior and endpoints.
.np_reg_interval_request <- function(dots) {
  supplied <- names(dots)
  interval <- dots[["interval", exact=TRUE]]
  band <- dots[["band", exact=TRUE]]
  control <- dots[["boot.control", exact=TRUE]]
  if (inherits(control, "np_boot_control") && "calibration" %in% names(control)) {
    previous <- .np_plot_boot_calibration(control)
    if (!is.null(interval)) stop("conflicting regression constructions: use interval and rebuild boot.control without calibration",call.=FALSE)
    replacement <- switch(previous, "default"="omit interval", "legacy"="use bias.control=np_bias_control(construction='legacy')",
      "absolute"="use interval='symmetric-sd'", "separate-tails"="use interval='asymmetric-sd'",
      "tail-envelope"="use interval='asymmetric-sd-envelope'", "smoothed-percentile"="use interval='outer-calibrated-percentile'",
      "percentile-tails"="transported density tails are unavailable for regression", paste0("use interval='",previous,"'"))
    stop(paste0("rebuild boot.control without its calibration field and ",replacement),call.=FALSE)
  }
  if (is.character(band) && length(band) == 1L && !is.na(band) && band == "pmzsd")
    stop("for regression use interval='normal-sd', band='pointwise' in place of band='pmzsd'", call.=FALSE)
  request <- list(interval=interval, supplied=intersect(supplied, c("interval","band","boot.control")),
                  band=band)
  dots$interval <- NULL
  if (is.null(interval)) return(list(dots=dots, request=request))

  choices <- c("normal-sd","symmetric-sd","asymmetric-sd","asymmetric-sd-envelope",
               "basic-envelope","basic-local-envelope","percentile","outer-calibrated-percentile")
  if (!is.character(interval) || length(interval) != 1L || is.na(interval) ||
      !interval %in% choices)
    stop("interval must be a single supported regression interval construction", call.=FALSE)
  interval <- as.character(interval)
  method <- dots[["errors", exact=TRUE]]
  if (is.null(method)) method <- "none"
  method <- .np_plot_scalar_match(method, c("none","bootstrap","asymptotic"), "errors")
  if (method == "none") stop("an explicit interval requires bootstrap or asymptotic errors", call.=FALSE)
  scope <- if (is.null(band) && !"band" %in% supplied) "pointwise" else
    .np_plot_scalar_match(band,c("pointwise","bonferroni","simultaneous","all"),"band")
  if (is.null(scope)) stop("an explicit interval requires a valid band scope",call.=FALSE)
  center <- dots[["center", exact=TRUE]]
  corrected <- !is.null(center) && identical(.np_plot_scalar_match(center,
    c("estimate","bias-corrected"),"center"), "bias-corrected")
  bias <- dots[["bias.control", exact=TRUE]]
  if (!is.null(bias) && !inherits(bias, "np_bias_control"))
    stop("regression bias.control must be created by np_bias_control()", call.=FALSE)
  legacy <- identical(bias[["construction", exact=TRUE]], "legacy")
  iterated <- corrected && !legacy
  if (corrected && legacy && method=="asymptotic")
    stop("legacy corrected intervals require bootstrap errors",call.=FALSE)
  empirical <- identical(method, "bootstrap")
  calibration <- NULL
  available <- switch(interval,
    "normal-sd"=scope=="pointwise" || (!empirical && !corrected && scope=="bonferroni"),
    "percentile"=empirical && !iterated,
    "symmetric-sd"=iterated || (!empirical && !corrected && scope=="simultaneous"),
    empirical && iterated)
  if (!available)
    stop(sprintf("interval='%s' is unavailable for this regression center, errors and band='%s'",interval,scope),call.=FALSE)
  if (interval == "normal-sd" && scope == "pointwise") scope <- "pmzsd"
  if (empirical && iterated && interval != "normal-sd")
    calibration <- switch(interval, "symmetric-sd"="absolute", "asymmetric-sd"="separate-tails",
      "asymmetric-sd-envelope"="tail-envelope", "outer-calibrated-percentile"="smoothed-percentile", interval)
  dots$band <- scope
  if (!is.null(calibration)) dots[["plot.errors.boot.calibration"]] <- calibration
  list(dots=dots, request=request)
}

.np_reg_interval_components <- function(meta, recipe, keys, requested, raw.inference, intervals, gradient, grid.size) {
  components <- list()
  for (key in keys) {
    public <- if (key=="pmzsd") "pointwise" else key
    states <- unique(intervals$status[intervals$band==key])
    if (recipe$iterated && !length(raw.inference) && identical(states, "undefined")) {
      components[[public]] <- list(interval=NA_character_, calibration=NULL, status=states,
        endpoint_construction=NULL, studentization=NULL)
      next
    }
    normal <- key=="pmzsd" || (recipe$method=="asymptotic" && !recipe$iterated &&
                             key %in% c("pointwise","bonferroni"))
    outer <- identical(recipe$calibration,"smoothed-percentile")
    interval <- if (normal) "normal-sd" else if (outer) "outer-calibrated-percentile" else
      if (recipe$iterated) switch(meta$calibration,
        "separate-tails"="asymmetric-sd", "asymmetric-sd"="asymmetric-sd",
        "tail-envelope"="asymmetric-sd-envelope", "absolute"="symmetric-sd",
        "symmetric-sd"="symmetric-sd", "basic-envelope"="basic-envelope",
        "basic-local-envelope"="basic-local-envelope",
        stop("regression interval publication has an unknown calibration",call.=FALSE)) else
      if (recipe$method=="asymptotic") "symmetric-sd" else "percentile"
    if (!is.null(requested) && !identical(interval,as.character(requested)))
      stop("regression interval provider returned a different construction",call.=FALSE)
    resolution <- meta$resolution
    if (is.data.frame(resolution)) resolution <- resolution[resolution$band==key,,drop=FALSE]
    calibration <- if (normal) "normal critical value" else if (interval=="percentile")
      "empirical curve ranks" else meta$calibration
    endpoint <- if (normal) "center plus or minus a normal critical value times SD" else
      if (!is.null(meta$endpoint_construction)) meta$endpoint_construction else switch(interval,
        "percentile"="empirical bootstrap curve quantiles",
        "symmetric-sd"="absolute Gaussian or bootstrap root bounds",
        "asymmetric-sd"="signed standardized root inversion",
        "asymmetric-sd-envelope"="envelope of signed and absolute root intervals",
        "basic-envelope"="basic root envelope", "basic-local-envelope"="local-scale basic root envelope",
        "outer-calibrated-percentile"="outer-calibrated empirical residual percentile bounds")
    if (outer) {
      # The mean provider keeps nominal marginal quantiles on multirow grids.
      # Both providers collapse all scopes to the joint hull for a singleton.
      endpoint <- if (grid.size == 1L) {
        if (gradient) "single-target E6 outer-calibrated marginal and family hull" else
          "single-target E1 outer-calibrated companion interval"
      } else if (gradient) switch(key,
        pointwise="unshifted corrected curves with E6 outer-calibrated marginal ranks",
        simultaneous="E6 outer-calibrated family interval enclosing its marginal intervals",
        bonferroni="nominal Bonferroni interval enclosing the E6 outer-calibrated family interval")
      else switch(key,
        pointwise="nominal marginal quantiles of unshifted corrected curves",
        simultaneous="E1 outer-calibrated family interval enclosing nominal marginal intervals",
        bonferroni="nominal Bonferroni interval enclosing the E1 outer-calibrated family interval")
      calibration <- if (!gradient && grid.size > 1L && key == "pointwise")
        "nominal marginal ranks" else if (gradient) "E6 outer rank calibration" else
        "E1 outer rank calibration"
    }
    student <- if (normal) "normal critical value applied to the reported SD or standard error" else meta$studentization
    if (is.null(student) && !normal)
      student <- raw.inference[["covariance",exact=TRUE]][["studentization",exact=TRUE]]
    # Scalar Gaussian-maxima metadata belongs to the joint component.
    critical <- if (!normal && (key=="simultaneous" || length(keys)==1L)) meta$critical_value else NULL
    components[[public]] <- list(interval=interval,calibration=calibration,status=states,
      endpoint_construction=endpoint,studentization=student,critical_value=critical,
      resolution=resolution,calibration_status=meta$calibration_status)
    if (!is.null(meta$resolution_limited) && (key=="simultaneous" || length(keys)==1L))
      components[[public]]$resolution_limited <- meta$resolution_limited
    if (!is.null(meta$critical_value_raw) && (key=="simultaneous" || length(keys)==1L))
      components[[public]]$critical_value_raw <- meta$critical_value_raw
  }
  meta[c("calibration","endpoint_construction","studentization","critical_value","resolution",
         "calibration_status","resolution_limited","critical_value_raw")] <- NULL
  if (identical(recipe$method,"bootstrap") && identical(recipe$bootstrap,"wild"))
    meta$wild_multiplier <- recipe$wild
  meta$components <- components
  meta
}
