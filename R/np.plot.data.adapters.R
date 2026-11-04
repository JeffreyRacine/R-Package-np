# Only plot-engine payloads enter this adapter. It is deliberately separate from
# fitted-object methods: several historical error slots have plot-only meanings.
.np_plot_data_field <- function(x, name, default = NULL) {
  value <- x[[name, exact = TRUE]]
  if (is.null(value)) default else value
}

# Publish one already-computed scalar/component payload before legacy slots can
# overwrite another component or omit undefined rows. This does no inference.
.np_plot_data_publish_target <- function(object, key, corrected, std_error,
                                         bands, boxplot) {
  object[[".plot_targets"]][[key]] <- list(corrected = corrected,
    std_error = std_error, bands = bands, boxplot = boxplot)
  object
}

.np_plot_data_publish_display <- function(value, coordinates, roles, response = NULL,
                                          overlay = FALSE, rug = FALSE, style = list(),
                                          points_args = list(), boxplot = FALSE, interval_style = list()) {
  observation <- NULL
  if (isTRUE(overlay) || isTRUE(rug)) {
    if (isTRUE(overlay) && is.data.frame(response)) {
      if (ncol(response) != 1L) stop("plot overlay requires a scalar response", call. = FALSE)
      response <- response[[1L]]
    }
    observation <- list(coordinates = coordinates, roles = roles,
      response = if (isTRUE(overlay)) as.double(response) else NULL)
  }
  attr(value, "np_plot_display") <- list(observation = observation,
    overlay = isTRUE(overlay), rug = isTRUE(rug),
    style = style[!vapply(style, is.null, logical(1L))],
    points_args = points_args, boxplot = isTRUE(boxplot), interval_style = interval_style)
  value
}

.np_plot_data_display <- function(raw, coordinates, default) {
  published <- attr(raw, "np_plot_display", exact = TRUE)
  if (is.null(published)) return(default)
  varying <- which(vapply(coordinates$grid, function(x) length(unique(x)) > 1L, logical(1L)))
  if (!length(varying)) varying <- 1L
  default$columns <- varying
  default$style <- published$style
  default$points_args <- published$points_args
  default$boxplot_enabled <- published$boxplot
  default$interval_style <- published$interval_style
  observation <- published$observation
  if (!is.null(observation)) {
    # Providers retain input order within X/Y/Z. Display names can be
    # syntactically repaired by their legacy data.frame calls; they are not IDs.
    role_keys <- function(roles) paste(roles, ave(seq_along(roles), roles, FUN = seq_along), sep = "\r")
    key <- role_keys(coordinates$roles)[varying]
    available <- role_keys(observation$roles)
    index <- match(key, available)
    if (anyNA(index)) stop("plot observations lack coordinate identity", call. = FALSE)
    geometry <- observation$coordinates[index]
    names(geometry) <- paste0("eval_", varying)
    if (published$overlay) default$overlay <- list(grid = geometry, value = observation$response)
    if (published$rug) default$rug <- geometry
  }
  default
}

.np_plot_data_contrast <- function(grid, index) {
  frames <- npCategoricalFirstDifferenceFrames(grid, index, "plot contrast metadata")
  list(coordinate = paste0("eval_", index), operation = "upper minus lower",
    lower_level = as.character(frames$lower[[index]]),
    upper_level = as.character(frames$upper[[index]]))
}

.np_plot_data_vector <- function(x, n, unavailable = FALSE) {
  if (is.null(x) || (unavailable && (length(x) == 0L || (length(x) == 1L && is.na(x)))))
    return(rep(NA_real_, n))
  if (!is.numeric(x) || length(x) != n)
    stop("plot target does not match its evaluation rows", call. = FALSE)
  as.double(x)
}

.np_plot_data_endpoints <- function(error, reference, encoding) {
  encoding <- match.arg(encoding, c("widths", "offsets", "endpoints"))
  if (!is.matrix(error) || ncol(error) != 2L || nrow(error) != length(reference))
    stop("plot interval dimensions do not match their reference", call. = FALSE)
  switch(encoding,
    widths = cbind(reference - error[, 1L], reference + error[, 2L]),
    offsets = cbind(reference + error[, 1L], reference + error[, 2L]),
    endpoints = unname(error))
}

.np_plot_data_band_matrix <- function(x, target = 1L) {
  dims <- dim(x)
  if (length(dims) == 3L) {
    if (target > dims[3L]) stop("missing plot band target", call. = FALSE)
    return(matrix(x[, , target, drop = FALSE], nrow = dims[1L], ncol = dims[2L]))
  }
  if (length(dims) != 2L) stop("invalid plot interval array", call. = FALSE)
  x
}

.np_plot_data_metadata <- function(raw, request, component = 1L, tau = 1L, quantile = FALSE) {
  inf <- .np_plot_data_field(raw, "inference", list())
  settings <- .np_plot_data_field(inf, "settings", list())
  calibration <- .np_plot_data_field(raw, "band.calibration", list())
  if (quantile && length(calibration) && is.null(calibration$method)) {
    key <- if (request$gradients) paste0("gradient", component) else "quantile"
    calibration <- .np_plot_data_field(calibration, key, list())
    calibration <- if (length(calibration) >= tau) calibration[[tau]] else list()
  } else if (isTRUE(request$gradients) && is.null(calibration$method) && length(calibration)) {
    calibration <- if (length(calibration) >= component) calibration[[component]] else list()
  }
  scalar <- function(x) if (is.atomic(x) && length(x) == 1L) x else NULL
  get <- function(name, fallback = NULL) {
    val <- scalar(.np_plot_data_field(inf, name))
    if (is.null(val)) val <- scalar(.np_plot_data_field(settings, name))
    if (is.null(val)) val <- scalar(.np_plot_data_field(calibration, name))
    if (is.null(val)) fallback else val
  }
  # Compact, value-only provenance. Never export the provider object itself.
  select_values <- function(x, keys) {
    out <- list()
    for (key in keys) {
      z <- .np_plot_data_field(x, key)
      if (is.atomic(z) && is.null(dim(z))) out[[gsub(".", "_", key, fixed = TRUE)]] <- z
    }
    out
  }
  calibrated <- .np_plot_data_field(inf, "calibration", list())
  if (!is.list(calibrated)) calibrated <- list()
  if (length(calibrated) == 1L && is.list(calibrated[[1L]]) &&
      is.null(calibrated$selection)) calibrated <- calibrated[[1L]]
  detail <- select_values(calibrated, c("attainable", "diagnostic.only", "boundary.used",
    "support.used", "outer.completions", "order.index", "endpoint.source", "calibration"))
  detail$selection <- select_values(.np_plot_data_field(calibrated, "selection", list()),
    c("q", "k", "status"))
  selection <- select_values(.np_plot_data_field(settings, "bandwidth.procedure", list()),
    c("status", "rule", "objective", "bandwidth.units", "original.rule.replayed", "descriptor.version"))
  errors <- request$errors
  actual.law <- get("law", if (identical(errors, "bootstrap")) request$bootstrap else NULL)
  list(errors = errors, center = request$center,
    projection = list(requested = scalar(raw$proper.requested),
      applied = scalar(raw$proper.applied), method = scalar(raw$proper.method),
      status = scalar(raw$proper.status)),
    construction = get("construction"), correction = get("correction"),
    calibration = get("calibration.method", get("calibration.requested", get("calibration", get("method")))),
    resampling_law = if (identical(errors, "bootstrap")) actual.law else NULL,
    pilot_model = if (identical(actual.law, "conditional-pilot")) actual.law else NULL,
    reference_law = if (identical(errors, "asymptotic")) "Gaussian" else NULL,
    scope = get("scope", get("family", "finite grid within panel and scalar target")),
    alpha = request$alpha, draws = get("draws", get("B")),
    primary = get("primary"), critical_units = get("critical.units"),
    selection = selection, calibration_details = detail,
    outer_draws = get("outer"), coverage_status = get("coverage.status"),
    endpoint_construction = get("endpoint.construction"),
    studentizer_ratio_range = select_values(inf, "studentizer.ratio.range")[["studentizer_ratio_range"]],
    qualification = get("scope"), calibration_status = get("calibration.status"),
    critical_value = get("critical", get("gaussian_critical")),
    resolution_limited = get("resolution.limited"),
    resolution = .np_plot_data_field(inf, "resolution"),
    rank = get("rank"), grid_size = get("grid.size"),
    positive_variance = get("positive.variance"),
    critical_value_raw = get("critical.raw"),
    target = get("target"), marginal_scale = get("marginal.scale"),
    dependence = get("dependence"),
    conditioning = get("conditioning"),
    studentization = get("studentization"),
    covariance_provenance = scalar(.np_plot_data_field(inf, "covariance.provenance")))
}

.np_plot_data_request <- function(options, method = NULL) {
  defaults <- if (is.function(method)) formals(method) else list()
  get <- function(public, engine, fallback) {
    value <- options[[public, exact = TRUE]]
    if (is.null(value)) value <- options[[engine, exact = TRUE]]
    if (is.null(value) && !is.null(defaults[[engine]])) {
      expr <- defaults[[engine]]
      # Defaults here are scalar literals or c() alternatives, never fit calls.
      if (is.atomic(expr)) value <- expr
      else if (is.call(expr) && identical(expr[[1L]], as.name("c")))
        value <- eval(expr, baseenv())[1L]
    }
    if (is.null(value)) fallback else value[1L]
  }
  list(errors = get("errors", "plot.errors.method", "none"),
    band = get("band", "plot.errors.type", "pmzsd"),
    alpha = get("alpha", "plot.errors.alpha", .05),
    center = get("center", "plot.errors.center", "estimate"),
    bootstrap = get("bootstrap", "plot.errors.boot.method", "inid"),
    B = get("B", "plot.errors.boot.num", 9999L),
    gaussian_draws = get("gaussian_draws", "gaussian.draws", 9999L),
    gradients = isTRUE(options$gradients),
    gradient_order = .np_plot_data_field(options, "gradient.order", 1L),
    coef = isTRUE(options$coef), coef_index = options$coef.index,
    perspective = get("perspective", "perspective", TRUE),
    supplied = .np_plot_data_supplied_names(names(options)))
}

.np_plot_data_engine_family <- function(bws) {
  classes <- c("bandwidth", "dbandwidth", "rbandwidth", "conbandwidth",
    "condbandwidth", "plbandwidth", "sibandwidth", "scbandwidth", "lsqregressionbandwidth")
  families <- c("npudens", "npudist", "npreg", "npcdens", "npcdist",
                "npplreg", "npindex", "npscoef", "nplsqreg")
  i <- match(class(bws), classes, nomatch = 0L)
  i <- i[i != 0L]
  if (!length(i)) stop("unknown plot-engine family", call. = FALSE)
  families[i[1L]]
}

.np_plot_data_grid <- function(raw, family) {
  frames <- roles <- NULL
  add <- function(x, role) {
    if (is.null(x)) return(invisible(NULL))
    if (!is.data.frame(x)) x <- as.data.frame(x)
    if (!ncol(x)) return(invisible(NULL))
    if (is.null(frames)) frames <<- x else {
      if (nrow(frames) != nrow(x)) stop("plot coordinate row counts differ", call. = FALSE)
      frames <<- cbind(frames, x)
    }
    roles <<- c(roles, rep(role, ncol(x)))
  }
  if (family == "npindex") {
    add(data.frame(index = raw$index), "index")
  } else if (family == "npplreg") {
    add(raw$evalx, "X"); add(raw$evalz, "Z")
  } else if (family == "npscoef" && !is.data.frame(raw$eval)) {
    add(raw$eval$exdat, "X"); add(raw$eval$ezdat, "Z")
  } else if (family %in% c("npcdens", "npcdist", "npqreg", "nplsqreg")) {
    add(raw$xeval, "X"); add(raw$yeval, "Y")
  } else add(raw$eval, "X")
  if (is.null(frames)) stop("plot result has no evaluation grid", call. = FALSE)
  list(grid = frames, roles = roles)
}

.np_plot_data_from_engine <- function(raw, bws, options, method = NULL) {
  request <- .np_plot_data_request(options, method)
  family <- .np_plot_data_engine_family(bws)
  if (family == "npcdist" && isTRUE(options$quantreg)) family <- "npqreg"
  if (inherits(bws, "lsqregressionbandwidth")) family <- "nplsqreg"
  panels <- list(); targets <- .np_plot_data_targets(); inference <- list(); families <- list()
  if (family == "npplreg" && request$coef) {
    n <- length(raw$coefficients); id <- "target_1"
    values <- .np_plot_data_values(raw$coefficients, id, std_error = raw$coefficient.stderr)
    intervals <- .np_plot_data_empty_intervals()
    if (request$errors != "none" && all(is.finite(raw$coefficient.stderr)))
      intervals <- .np_plot_data_intervals(list(std_error = cbind(
        raw$coefficients - raw$coefficient.stderr, raw$coefficients + raw$coefficient.stderr)),
        id, NA_real_, status = "display_only")
    terms <- names(raw$fit$xcoef)
    if (is.null(terms)) terms <- bws$xnames
    panels$panel_1 <- .np_plot_data_panel(data.frame(term = terms), "term", values,
      intervals, list(kind = "coefficients"))
    targets <- .np_plot_data_target(id, "coefficient", "Parametric coefficient")
    return(.np_plot_data(family, request, panels, targets))
  }
  for (p in seq_along(raw)) {
    o <- raw[[p]]; panel_id <- paste0("panel_", p)
    coordinates <- .np_plot_data_grid(o, family); n <- nrow(coordinates$grid)
    gradient <- request$gradients
    quantile <- family %in% c("npqreg", "nplsqreg")
    gc <- grep("^gc[0-9]+err$", names(o), value = TRUE)
    components <- if (gradient && family %in% c("npcdens", "npcdist"))
      seq_len(ncol(o$congrad)) else p
    values <- intervals <- NULL
    boxes <- list()
    for (variable.index in components) {
    selected.gc <- paste0("gc", variable.index, "err")
    payload.key <- if (gradient) paste0("gradient", variable.index) else "value"
    publication <- o[[".plot_targets", exact = TRUE]][[payload.key]]
    field <- if (quantile) "quantile" else switch(family,
      npudens = "dens", npudist = "dist", npcdens = "condens", npcdist = "condist", "mean")
    estimate <- o[[field, exact = TRUE]]
    if (gradient) {
      if (quantile) {
        d <- dim(o$quantgrad)
        estimate <- if (length(d) == 3L) matrix(o$quantgrad[, variable.index, , drop = FALSE], nrow = n) else
          matrix(o$quantgrad[, variable.index], nrow = n)
      } else if (family %in% c("npcdens", "npcdist")) {
        estimate <- o$congrad[, variable.index]
      } else estimate <- o$grad
    }
    rows <- o[[".plot_rows", exact = TRUE]]
    if (!is.null(rows)) estimate <- rows$estimate
    if (is.null(estimate)) stop("plot result is missing its selected target", call. = FALSE)
    if (is.null(dim(estimate))) estimate <- matrix(estimate, ncol = 1L)
    if (!is.matrix(estimate) || nrow(estimate) != n)
      stop("plot target row loss requires a provider-boundary repair", call. = FALSE)
    nt <- ncol(estimate)
    inf <- .np_plot_data_field(o, "inference", list())
    for (j in seq_len(nt)) {
      meta <- .np_plot_data_metadata(o, request, variable.index, j, quantile)
      id <- paste0("target_", nrow(targets) + 1L); owner <- paste0("inference_", length(inference) + 1L)
      box <- if (!is.null(publication)) publication$boxplot else .np_plot_data_field(o, "bxp")
      if (quantile && nt > 1L && is.list(box) && is.null(box$stats))
        box <- if (length(box) >= j) box[[j]] else NULL
      if (is.list(box) && !is.null(box$stats)) boxes[[id]] <- box
      quantity <- if (gradient) "derivative" else if (quantile) "quantile" else
        if (family == "npscoef" && request$coef) "coefficient" else switch(family,
          npudens = "density", npudist = "cdf", npcdens = "conditional_density",
          npcdist = "conditional_cdf", "mean")
      variable <- NA_character_
      if (gradient) {
        idx <- if (family == "npindex") j else variable.index
        if (family == "npscoef") idx <- .np_plot_data_field(options, "coef.index", 1L)
        vn <- bws$xnames
        variable <- vn[idx]
        if (length(variable) != 1L || is.na(variable))
          stop("plot derivative has no variable identity", call. = FALSE)
      }
      tau <- if (quantile) o$tau[j] else NA_real_
      label <- if (quantile) paste0(quantity, " (tau=", format(tau), ")") else
        if (gradient) paste0("Derivative: ", variable) else quantity
      multi_index <- integer(); derivative.order <- NA_integer_
      if (gradient) {
        coordinate <- which(names(coordinates$grid) == variable & coordinates$roles == "X")
        categorical <- length(coordinate) == 1L && is.factor(coordinates$grid[[coordinate]])
        if (categorical) {
          quantity <- "contrast"
          label <- paste("Contrast:", variable)
          meta$contrast <- .np_plot_data_contrast(coordinates$grid, coordinate)
        } else {
          continuous <- bws$xnames
          numeric.x <- coordinates$roles == "X" & vapply(coordinates$grid, is.numeric, logical(1L))
          if (any(numeric.x)) continuous <- names(coordinates$grid)[numeric.x]
          derivative.order <- rep_len(as.integer(request$gradient_order), length(continuous))[match(variable, continuous)]
          if (family == "npindex") derivative.order <- 1L
          multi_index <- integer(length(bws$xnames))
          multi_index[match(variable, bws$xnames)] <- derivative.order
        }
      }
      term <- if (family == "npscoef" && request$coef)
        bws$xnames[min(.np_plot_data_field(options, "coef.index", 1L), bws$xndim)] else NA_character_
      targets <- rbind(targets, .np_plot_data_target(id, quantity, label, term = term,
        variable = variable, inference_id = owner, tau = tau,
        order = derivative.order, multi_index = multi_index))
      ref.name <- if (gradient && selected.gc %in% gc) sub("err$", "bias.corrected", selected.gc) else
        if (gradient) "gradient.bias.corrected" else "bias.corrected"
      ref.raw <- o[[ref.name, exact = TRUE]]
      if (!is.null(publication)) ref.raw <- publication$corrected
      if (!is.null(rows)) ref.raw <- rows$corrected
      if (request$center != "bias-corrected") ref.raw <- NULL
      if (length(dim(ref.raw)) == 2L) ref.raw <- ref.raw[, j]
      if (is.null(ref.raw) && quantile && !is.null(o$bias) && request$center == "bias-corrected")
        ref.raw <- estimate[, j] - if (is.matrix(o$bias)) o$bias[, j] else o$bias
      corrected <- .np_plot_data_vector(ref.raw, n, unavailable = TRUE)
      se.raw <- .np_plot_data_field(inf, "se", if (!is.null(publication)) publication$std_error else
        o[[".plot_std_error", exact = TRUE]])
      if (is.null(se.raw) && gradient && family %in% c("npcdens", "npcdist") &&
          is.matrix(o$congerr)) se.raw <- o$congerr[, variable.index]
      # The incumbent smooth-coefficient provider exposes mean SEs, not beta
      # SEs. Never label a mean SE as uncertainty for a coefficient curve.
      if (family == "npscoef" && request$coef) se.raw <- NULL
      if (is.matrix(se.raw)) se.raw <- se.raw[, j]
      se <- .np_plot_data_vector(se.raw, n, unavailable = TRUE)
      meta$std_error_reference <- if (!is.null(inf$se)) request$center else "estimate"
      value <- .np_plot_data_values(estimate[, j], id, corrected, se)
      reference <- if (request$center == "bias-corrected") corrected else estimate[, j]
      bands <- list()
      if (request$errors != "none") {
        requested <- if (request$band == "all") c("pointwise", "simultaneous", "bonferroni") else request$band
        available <- inf$bands
        encoding <- "endpoints"
        if (is.null(available)) {
          available <- if (!is.null(publication)) publication$bands else
            .np_plot_data_field(o, ".plot_band_widths", o$band.errors)
          encoding <- "widths"
          if (is.null(publication) && gradient && family %in% c("npcdens", "npcdist") &&
              !is.null(available) && !all(requested %in% names(available)))
            available <- available[[variable.index]]
          if (quantile && !is.null(available)) {
            if (!is.null(publication) || !is.null(o[[".plot_band_widths", exact = TRUE]])) {
              if (!all(requested %in% names(available))) available <- available[[j]]
            } else available <- available[[if (gradient) paste0("gradient", variable.index) else "quantile"]][[j]]
          }
        }
        if (!is.null(available) && all(requested %in% names(available))) {
          bands <- lapply(available[requested], function(a)
            .np_plot_data_endpoints(.np_plot_data_band_matrix(a, j), reference, encoding))
        } else {
          if (length(requested) != 1L)
            stop("all-band plot publication is incomplete at the provider boundary", call. = FALSE)
          error.field <- if (selected.gc %in% gc && gradient) selected.gc else if (quantile) "quanterr" else
            if (gradient) "gerr" else switch(family,npudens="derr",npudist="derr",
              npcdens="conderr",npcdist="conderr","merr")
          error <- o[[error.field, exact = TRUE]]
          encoding <- if (family %in% c("npplreg", "npscoef") && !gradient) "widths" else "offsets"
          # Surfaces historically use positive widths except modern signed UD payloads.
          surface <- isTRUE(request$perspective) && ncol(coordinates$grid) == 2L &&
            all(vapply(coordinates$grid, function(z) length(unique(z)) > 1L, TRUE))
          if (surface && family == "npreg") encoding <- "widths"
          if (surface && family %in% c("npudens", "npudist") && is.null(inf$bands)) encoding <- "widths"
          if (family == "npindex" && gradient) {
            error <- cbind(o$glerr[, j], o$gherr[, j]); encoding <- "endpoints"
          } else if (!is.null(rows)) {
            error <- rows$widths; encoding <- "widths"
          } else if (is.null(error) && all(is.na(estimate[, j]))) {
            error <- matrix(NA_real_, n, 2L)
          } else error <- .np_plot_data_band_matrix(error, j)
          bands[[requested]] <- .np_plot_data_endpoints(error, reference, encoding)
        }
      }
      family_id <- paste0("family_", length(families) + 1L)
      members <- data.frame(panel_id = panel_id, row_id = seq_len(n), target_id = id)
      families[[family_id]] <- list(event = "coverage of the requested finite target rows",
        scope = meta$scope, members = members)
      interval <- .np_plot_data_intervals(bands, id, request$alpha, family_id)
      if (isTRUE(meta$calibration_details$diagnostic_only) ||
          identical(meta$calibration_details$attainable, FALSE))
        interval$status <- rep("diagnostic_only", nrow(interval))
      values <- rbind(values, value); intervals <- rbind(intervals, interval)
      inference[[owner]] <- meta
    }
    }
    display <- list(kind = if (isTRUE(request$perspective) && ncol(coordinates$grid) == 2L &&
      all(vapply(coordinates$grid, function(z) length(unique(z)) > 1L, TRUE))) "surface" else "curve",
      boxplots = boxes)
    display <- .np_plot_data_display(raw, coordinates, display)
    panels[[panel_id]] <- .np_plot_data_panel(coordinates$grid, coordinates$roles,
      values, intervals, display)
  }
  .np_plot_data(family, request, panels, targets, inference, families)
}

.np_plot_data_finish <- function(value, bws, options, method, random_seed = 42L) {
  output <- .np_plot_data_field(options, "plot.behavior", "plot")
  if (identical(output, "plot")) {
    if (identical(options$renderer, "rgl")) return(value)
    return(invisible(NULL))
  }
  if (identical(method, .np_plot_compat_dispatch)) method <- .np_plot_engine_for_bws(bws)
  result <- .np_plot_data_from_engine(value, bws, options, method)
  result$request$random_seed <- random_seed
  result
}

.np_plot_data_frame_bands <- function(frame, request) {
  if (request$errors == "none") return(list())
  requested <- if (request$band == "all")
    c("pointwise", "simultaneous", "bonferroni") else request$band
  out <- setNames(vector("list", length(requested)), requested)
  for (band in requested) {
    keys <- paste0(band, c(".lower", ".upper"))
    if (!all(keys %in% names(frame)) && length(requested) == 1L)
      keys <- c("lower", "upper")
    if (!all(keys %in% names(frame)))
      stop("plot frame is missing requested endpoints", call. = FALSE)
    out[[band]] <- as.matrix(frame[keys])
  }
  out
}

.np_plot_data_from_conmode <- function(raw, request, gradients = FALSE) {
  panels <- list(); targets <- .np_plot_data_targets(); inference <- families <- list()
  surface <- !is.null(raw$surface) && !is.null(raw$grid)
  frames <- if (surface) list(raw$surface) else raw
  for (i in seq_along(frames)) {
    d <- frames[[i]]; panel_id <- paste0("panel_", i); id <- paste0("target_", i)
    owner <- paste0("inference_", i); family_id <- paste0("family_", i)
    if (!is.data.frame(d)) stop("invalid class-probability plot frame", call. = FALSE)
    if (surface) {
      grid <- data.frame(raw$grid$values[[1L]][d$grid_index1],
                         raw$grid$values[[2L]][d$grid_index2])
      names(grid) <- raw$grid$variables
      if (ncol(raw$held)) grid <- cbind(grid, raw$held[rep(1L, nrow(d)), , drop = FALSE])
      display <- list(kind = "surface", dimensions = raw$grid$dims, coordinates = c("eval_1", "eval_2"))
    } else {
      grid <- attr(d, "plot.grid", exact = TRUE)
      if (is.null(grid)) {
        grid <- data.frame(d$x); names(grid) <- as.character(d$variable[1L])
      }
      display <- list(kind = "curve", coordinates = paste0("eval_", match(as.character(d$variable[1L]), names(grid))))
    }
    estimate <- d[[if (gradients) "effect" else "probability"]]
    corrected <- if (request$center == "bias-corrected") d$center else rep(NA_real_, nrow(d))
    se <- .np_plot_data_vector(d$stderr, nrow(d), unavailable = TRUE)
    values <- .np_plot_data_values(estimate, id, corrected, se)
    intervals <- .np_plot_data_intervals(.np_plot_data_frame_bands(d, request),
      id, request$alpha, family_id)
    families[[family_id]] <- list(event = "coverage of the requested target rows",
      scope = "finite grid within panel and response level",
      members = data.frame(panel_id = panel_id, row_id = seq_len(nrow(d)), target_id = id))
    level <- if (!is.null(d$level)) as.character(d$level[1L]) else as.character(raw$level)
    quantity <- if (gradients && is.factor(d$x)) "probability_contrast" else
      if (gradients) "probability_derivative" else "class_probability"
    targets <- rbind(targets, .np_plot_data_target(id, quantity,
      paste(quantity, level), variable = if (gradients) as.character(d$variable[1L]) else NA_character_,
      level = level, order = if (gradients && !is.factor(d$x)) 1L else NA_integer_,
      inference_id = owner))
    inference[[owner]] <- .np_plot_data_metadata(
      list(band.calibration = attr(d, "band.calibration", exact = TRUE)), request)
    if (identical(quantity, "probability_contrast"))
      inference[[owner]]$contrast <- .np_plot_data_contrast(grid,
        match(as.character(d$variable[1L]), names(grid)))
    panels[[panel_id]] <- .np_plot_data_panel(grid, rep("X", ncol(grid)), values, intervals, display)
  }
  .np_plot_data("npconmode", request, panels, targets, inference, families)
}

.np_plot_data_from_copula <- function(raw, request, density = FALSE, view = "surface") {
  frames <- if (is.data.frame(raw)) setNames(list(raw), view) else raw
  panels <- list(); targets <- .np_plot_data_targets(); inference <- families <- list()
  for (i in seq_along(frames)) {
    d <- frames[[i]]; panel_id <- paste0("panel_", i); id <- paste0("target_", i)
    owner <- paste0("inference_", i); family_id <- paste0("family_", i)
    if (!is.data.frame(d) || !all(c("u1", "u2") %in% names(d)))
      stop("invalid copula plot coordinates", call. = FALSE)
    empirical <- identical(names(frames)[i], "empirical")
    this.density <- density || identical(names(frames)[i], "density")
    grid <- d[c("u1", "u2")]
    estimate <- if (empirical) rep(NA_real_, nrow(d)) else d$copula
    corrected <- if (!empirical && request$center == "bias-corrected") d$center else rep(NA_real_, nrow(d))
    values <- .np_plot_data_values(estimate, id, corrected)
    if (empirical) values$status <- "display_only"
    bands <- if (empirical) list() else .np_plot_data_frame_bands(d, request)
    intervals <- .np_plot_data_intervals(bands, id, request$alpha, family_id)
    quantity <- if (empirical) "empirical_probability_pair" else
      if (this.density) "copula_density" else "copula_cdf"
    targets <- rbind(targets, .np_plot_data_target(id, quantity, inference_id = owner))
    families[[family_id]] <- list(event = "coverage of the requested target rows",
      scope = "finite probability grid", members = data.frame(panel_id = panel_id,
        row_id = seq_len(nrow(d)), target_id = id))
    inference[[owner]] <- .np_plot_data_metadata(
      list(band.calibration = attr(d, "band.calibration", exact = TRUE)), request)
    if (empirical) inference[[owner]]$errors <- "none"
    display <- list(kind = if (empirical || view == "sample") "scatter" else "surface",
      view = if (view == "all") if (this.density) "surface" else if (empirical) "empirical" else "contour" else view)
    panels[[panel_id]] <- .np_plot_data_panel(grid, rep("probability", 2L), values, intervals, display)
  }
  .np_plot_data("npcopula", request, panels, targets, inference, families)
}

.np_plot_data_from_iv <- function(spec, variable) {
  grid <- data.frame(spec$curve.args$x)
  names(grid) <- variable
  id <- "target_1"
  target <- .np_plot_data_target(id, if (spec$gradients) "structural_derivative" else "structural_mean",
    variable = if (spec$gradients) names(grid) else NA_character_,
    order = if (spec$gradients) 1L else NA_integer_)
  display <- list(kind = "iv", curve_args = spec$curve.args[setdiff(names(spec$curve.args), c("x", "y"))],
    frame_args = spec$frame.args, continuous = spec$continuous,
    overlay = if (spec$overlay) data.frame(x = spec$overlay.x, y = spec$overlay.y) else NULL,
    overlay_args = spec$overlay.args, rug = if (spec$rug) spec$rug.x else NULL)
  panel <- .np_plot_data_panel(grid, "Z", .np_plot_data_values(spec$curve.args$y, id),
    .np_plot_data_empty_intervals(), display)
  .np_plot_data(spec$family, list(errors = "none", gradients = spec$gradients),
    list(panel_1 = panel), target)
}
