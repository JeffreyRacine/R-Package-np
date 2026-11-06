# Pre-seed refusal only. This is not a general RNG capability certificate.
# Unknown data reconstruction and other producers retain their own validation.
.np_plot_replay_training_n <- function(dots, bws = NULL) {
  # Match the engine's suppliedness and retained-value precedence without
  # evaluating calls, formulas, NA policies or historical row counts.
  if (!any(c("xdat", "ydat") %in% names(dots)) && !is.null(bws)) {
    if (!is.null(bws[["formula", exact = TRUE]])) {
      training <- bws[[".np.formula.training", exact = TRUE]]
      if (is.object(training) || !is.list(training) ||
          !"na.action" %in% names(training)) return(NULL)
      frame <- training[["frame", exact = TRUE]]
      if (!identical(class(frame), "data.frame") || length(frame) < 2L)
        return(NULL)
      tt <- attr(frame, "terms", exact = TRUE)
      if (!identical(class(tt), c("terms", "formula")) ||
          !identical(attr(tt, "response", exact = TRUE), 1L)) return(NULL)
      labels <- attr(tt, "term.labels", exact = TRUE)
      columns <- names(frame)[-1L]
      if (!identical(labels, columns)) {
        if (is.object(labels) || !is.character(labels) || anyNA(labels))
          return(NULL)
        # terms quotes simple non-syntactic symbols; model.frame does not.
        # Encode known names, without parsing or evaluating supplied terms.
        symbols <- vapply(columns, function(name)
          paste(deparse(as.name(name), backtick = TRUE), collapse = ""),
          character(1L), USE.NAMES = FALSE)
        if (!identical(labels, symbols)) return(NULL)
      }
      dots <- list(xdat = frame[-1L], ydat = frame[[1L]])
    } else {
      training <- bws[[".np.native.training", exact = TRUE]]
      if (is.object(training) || !is.list(training) ||
          !all(c("xdat", "ydat") %in% names(training))) return(NULL)
      dots <- training
    }
  }
  x <- dots[["xdat", exact = TRUE]]
  y <- dots[["ydat", exact = TRUE]]
  plain <- function(z) {
    if (is.object(z)) return((identical(class(z), "factor") ||
      identical(class(z), c("ordered", "factor"))) && !anyNA(unclass(z)))
    is.numeric(z) && is.null(dim(z)) && all(is.finite(z))
  }
  if (identical(class(x), "data.frame")) {
    if (!length(x) || !all(vapply(x, plain, logical(1L)))) return(NULL)
    n <- nrow(x)
  } else if (!is.object(x) && is.matrix(x) && is.numeric(x) && all(is.finite(x))) {
    n <- nrow(x)
  } else if (!is.object(x) && is.numeric(x) && is.null(dim(x)) && all(is.finite(x))) {
    n <- length(x)
  } else return(NULL)
  if (is.object(y) || !is.numeric(y) || !is.null(dim(y)) ||
      !all(is.finite(y)) || length(y) != n || n < 1L) return(NULL)
  n
}

.np_plot_replay_cdf_targets <- function(bws, dots) {
  if (!identical(bws[["type", exact = TRUE]], "fixed") ||
      !isTRUE(bws[["ncon", exact = TRUE]] %in% 1:2) ||
      !isTRUE(bws[["nuno", exact = TRUE]] == 0) ||
      !isTRUE(bws[["nord", exact = TRUE]] == 0)) return(NULL)
  m <- dots[["neval", exact = TRUE]]
  if (is.null(m)) m <- 50L
  if (is.object(m) || !is.numeric(m) || length(m) != 1L ||
      !is.finite(m) || m < 1 || m != floor(m) || m > .Machine$integer.max) return(NULL)
  perspective <- dots[["perspective", exact = TRUE]]
  if (is.null(perspective)) perspective <- TRUE
  if (!identical(perspective, TRUE) && !identical(perspective, FALSE)) return(NULL)
  if (bws[["ncon", exact = TRUE]] == 2L && perspective) m <- as.double(m)^2
  m
}

.np_plot_replay_uniform_state_absent <- function(kind) {
  # A header without payload establishes missing replay state, not the converse.
  # Do not inspect active bindings or try to discover arbitrary native state.
  if (!identical(kind[1L], "user-supplied") ||
      !exists(".Random.seed", .GlobalEnv, inherits = FALSE) ||
      bindingIsActive(".Random.seed", .GlobalEnv)) return(FALSE)
  seed <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  is.integer(seed) && length(seed) == 1L
}

# Explicit factor providers use a different covariance cutoff from joint CDF.
# Only exact, already-qualified controls and geometry are classified here.
.np_plot_replay_factor_preflight <- function(method, bws, dots) {
  if (!identical(dots[["plot.errors.method", exact = TRUE]], "asymptotic"))
    return(invisible(NULL))
  control <- dots[["bias.control", exact = TRUE]]
  calibration <- dots[["plot.errors.boot.calibration", exact = TRUE]]
  if (!is.null(calibration) && !identical(calibration, "default"))
    return(invisible(NULL))
  cdf <- inherits(bws, "dbandwidth") &&
    (identical(method, .np_plot_dbandwidth_engine) || identical(method, .np_plot_compat_dispatch)) &&
    identical(control, structure(list(construction = "cdf", version = 1L),
                                 class = "np_udist_bias_control"))
  density <- inherits(bws, "bandwidth") &&
    (identical(method, .np_plot_bandwidth_engine) || identical(method, .np_plot_compat_dispatch)) &&
    (identical(control, structure(list(construction = "joint-scale", calibration = "absolute", version = 3L),
                                  class = "np_udens_bias_control")) ||
     identical(control, structure(list(construction = "positive-covariance", calibration = "absolute", version = 2L),
                                  class = "np_udens_bias_control")))
  band <- dots[["plot.errors.type", exact = TRUE]]
  if ((!cdf || !isTRUE(band %in% c("simultaneous", "all"))) &&
      (!density || !isTRUE(band %in% c("pointwise", "bonferroni", "pmzsd", "simultaneous", "all"))))
    return(invisible(NULL))
  # Explicit density calibration is only valid for bootstrap inference.
  if (density && !is.null(calibration)) return(invisible(NULL))
  if (!identical(bws[["ckertype", exact = TRUE]], "gaussian") ||
      !isTRUE(bws[["ckerorder", exact = TRUE]] == 2L) ||
      !identical(bws[["ckerbound", exact = TRUE]], "none")) return(invisible(NULL))
  m <- .np_plot_replay_cdf_targets(bws, dots)
  if (is.null(m) || m <= if (density) 512L else 32L) return(invisible(NULL))
  kind <- RNGkind()
  if (.np_plot_replay_uniform_state_absent(kind))
    stop("unsafe RNG replay: this explicit factor plot requires custom uniform state absent from .Random.seed; explicitly select a qualified built-in generator before plotting", call. = FALSE)
  builtin <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister", "Knuth-TAOCP", "Knuth-TAOCP-2002", "L'Ecuyer-CMRG")
  if (kind[1L] %in% builtin && identical(kind[2L], "Box-Muller"))
    stop("unsafe RNG replay: this explicit factor plot requires Box-Muller state not retained by .Random.seed; explicitly select a qualified normal generator such as RNGkind(normal.kind = 'Inversion') before plotting", call. = FALSE)
  invisible(NULL)
}

.np_plot_rng_replay_preflight <- function(method, bws, dots) {
  if (!identical(dots[["plot.errors.center", exact = TRUE]], "bias-corrected") ||
      !isTRUE(dots[["plot.errors.method", exact = TRUE]] %in% c("asymptotic", "bootstrap")))
    return(invisible(NULL))
  .np_plot_replay_factor_preflight(method, bws, dots)
  regression <- inherits(bws, "rbandwidth") &&
    (identical(method, .np_plot_rbandwidth_engine) || identical(method, .np_plot_compat_dispatch))
  distribution <- inherits(bws, "dbandwidth") &&
    (identical(method, .np_plot_dbandwidth_engine) || identical(method, .np_plot_compat_dispatch))
  if (!regression && !distribution) return(invisible(NULL))
  calibration <- dots[["plot.errors.boot.calibration", exact = TRUE]]
  if (isTRUE(calibration %in% c("legacy", "smoothed-percentile"))) return(invisible(NULL))
  control <- dots[["bias.control", exact = TRUE]]
  if (!is.null(control) && (!is.list(control) ||
      !inherits(control, if (regression) "np_bias_control" else "np_udist_bias_control") ||
      !identical(control[["construction", exact = TRUE]],
                 if (regression) "iterated" else "joint-process"))) return(invisible(NULL))
  kind <- RNGkind()
  if (.np_plot_replay_uniform_state_absent(kind)) {
    asymptotic <- identical(dots[["plot.errors.method", exact = TRUE]], "asymptotic")
    law <- dots[["plot.errors.boot.method", exact = TRUE]]
    if (is.null(law)) law <- if (regression) "wild" else "inid"
    unsafe <- regression && (asymptotic ||
      isTRUE(law %in% c("wild", "inid", "fixed", "geom")))
    if (distribution) {
      m <- .np_plot_replay_cdf_targets(bws, dots)
      band <- dots[["plot.errors.type", exact = TRUE]]
      unsafe <- !is.null(m) && m > 25L &&
        ((asymptotic && isTRUE(band %in% c("bonferroni", "simultaneous", "all"))) ||
         (!asymptotic && identical(law, "inid")))
    }
    if (unsafe)
      stop("unsafe RNG replay: this bias-corrected plot requires custom uniform state absent from .Random.seed; explicitly select a qualified built-in generator before plotting", call. = FALSE)
  }
  if (!identical(dots[["plot.errors.method", exact = TRUE]], "asymptotic")) return(invisible(NULL))
  builtin <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister", "Knuth-TAOCP", "Knuth-TAOCP-2002", "L'Ecuyer-CMRG")
  if (!kind[1L] %in% builtin || !identical(kind[2L], "Box-Muller")) return(invisible(NULL))
  unsafe <- FALSE
  if (regression) {
    n <- .np_plot_replay_training_n(dots, bws)
    B <- dots[["gaussian.draws", exact = TRUE]]
    if (is.null(B)) B <- 9999L
    if (!is.null(n) && !is.object(B) && is.numeric(B) && length(B) == 1L &&
        is.finite(B) && B >= 1 && B == floor(B) && B <= .Machine$integer.max)
      # Same bounded-value boundary as .np_bias_gaussian_cache; no allocation.
      unsafe <- 8 * as.double(n) * B > 16 * 1024^2
  } else {
    band <- dots[["plot.errors.type", exact = TRUE]]
    m <- .np_plot_replay_cdf_targets(bws, dots)
    # Joint CDF factors use target tiles of 25; pointwise/pmzsd are deterministic.
    unsafe <- !is.null(m) && m > 25L &&
      isTRUE(band %in% c("bonferroni", "simultaneous", "all"))
  }
  if (unsafe)
    stop("unsafe RNG replay: this bias-corrected plot requires Box-Muller state not retained by .Random.seed; explicitly select a qualified normal generator such as RNGkind(normal.kind = 'Inversion') before plotting", call. = FALSE)
  invisible(NULL)
}
