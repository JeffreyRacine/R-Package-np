# Translate only the two regression plot entries. The existing engine owns
# law selection, kernel eligibility, draws, finite-B behavior and endpoints.
.np_reg_interval_request <- function(dots) {
  supplied <- names(dots)
  interval <- dots[["interval", exact=TRUE]]
  band <- dots[["band", exact=TRUE]]
  control <- dots[["boot.control", exact=TRUE]]
  if (inherits(control, "np_boot_control") && "calibration" %in% names(control))
    stop("regression calibration has moved to interval; rebuild boot.control using only nonfixed, wild and blocklen, and select the construction with interval", call.=FALSE)
  if (is.character(band) && length(band) == 1L && !is.na(band) && band == "pmzsd")
    stop("for regression use interval='normal-sd', band='pointwise' in place of band='pmzsd'", call.=FALSE)
  request <- list(interval=interval, supplied=intersect(supplied, c("interval","band")),
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
  legacy <- identical(dots[["bias.control", exact=TRUE]][["construction", exact=TRUE]], "legacy")
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
