## This function will compute the cumulative integral at each sample
## realization using the trapezoidal rule and the cumsum function as
## we need to compute this in a computationally efficient manner.

# integrate.trapezoidal <- function(x,y) {
#   n <- length(x)
#   rank.x <- rank(x)
#   order.x <- order(x)
#   y <- y[order.x]
#   x <- x[order.x]
#   int.vec <- numeric(length(x))
#   ## Use a correction term at the boundary: -cx^2/12*(f'(b)-f'(a)),
#   ## check for NaN case
#   cx  <- x[2]-x[1]
#   ca <- (y[2]-y[1])/cx
#   cb <- (y[n]-y[n-1])/cx
#   cf <- cx^2/12*(cb-ca)
#   if(!is.finite(cf)) cf <- 0
#   int.vec[1] <- 0
#   int.vec[2:n] <- cumsum((x[2:n]-x[2:n-1])*(y[2:n]+y[2:n-1])/2)
#   return(int.vec[rank.x]-cf)
# }

# Benches 1.3-1.7 times faster (mean/median) produces identical results

integrate.trapezoidal <- function(x, y) {
  n <- length(x)
  order.x <- order(x)
  x <- x[order.x]
  y <- y[order.x]
  dx <- diff(x)
  dy <- diff(y)
  cx <- dx[1]
  ca <- dy[1] / cx
  cb <- dy[n - 1] / cx
  cf <- cx^2 / 12 * (cb - ca)
  if (!is.finite(cf)) cf <- 0
  int.vec <- c(0, cumsum(dx * (y[-n] + y[-1]) / 2))
  int.vec <- int.vec - cf
  int.vec[order(order.x)]  # inverse permutation == rank(x) when no ties
}

## No Zero Denominator, used in C code for kernel estimation...

## Original, better, best

# NZD <- function(a) {
#   ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
# }
  
# NZD <- function(a) {
#   if(length(a) == 1) {
#     if(is.na(a)) return(a)
#     if(a < 0) return(min(-.Machine$double.eps, a))
#     return(max(.Machine$double.eps, a))
#   }
#   ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
# }

# Benches 5.5-7.9 times faster (mean/median) produces identical results

NZD <- function(a) {
  eps <- .Machine$double.eps
  if (length(a) == 1) {
    if (a >= 0) {
      if (a < eps) return(eps)
    } else {
      if (a > -eps) return(-eps)
    }
    return(a)
  }
  idx <- which(abs(a) < eps)
  if (length(idx) > 0) {
    small <- a[idx]
    small[small >= 0] <- eps
    small[small < 0] <- -eps
    a[idx] <- small
  }
  a
}

## New function when only positive values are guaranteed, better, best

# NZD_pos <- function(a) {
#   if(length(a) == 1) {
#     if(is.na(a)) return(a)
#     return(max(.Machine$double.eps, a))
#   }
#   pmax(.Machine$double.eps,a)
# }

# Benches 1.3-1.8 times faster (mean/median) produces identical results

NZD_pos <- function(a) {
  eps <- .Machine$double.eps
  if (length(a) == 1)
    return(if (a < eps) eps else a)
  idx <- which(a < eps)
  if (length(idx) > 0)
    a[idx] <- eps
  a
}

npRidgeSequenceAdditive <- function(n.train, cap = 1.0) {
  if (!is.numeric(n.train) || length(n.train) != 1L || is.na(n.train) ||
      !is.finite(n.train) || n.train < 1 || n.train != floor(n.train))
    stop("'n.train' must be a positive integer")
  if (!is.numeric(cap) || length(cap) != 1L || is.na(cap) ||
      !is.finite(cap) || cap <= 0)
    stop("'cap' must be a positive finite numeric scalar")

  step <- 1.0 / as.double(n.train)
  seq.out <- seq.int(from = 0.0, to = as.double(cap), by = step)
  if (utils::tail(seq.out, 1L) < cap)
    seq.out <- c(seq.out, as.double(cap))
  unique(as.double(seq.out))
}

npRidgeSequenceFromBase <- function(n.train, ridge.base = 0.0, cap = 1.0) {
  if (!is.numeric(ridge.base) || length(ridge.base) != 1L || is.na(ridge.base) ||
      !is.finite(ridge.base) || ridge.base < 0)
    stop("'ridge.base' must be a non-negative finite numeric scalar")

  base <- as.double(ridge.base)
  step <- 1.0 / as.double(n.train)
  cap <- as.double(cap)

  if (base >= cap)
    return(unique(c(base, cap)))

  seq.out <- seq.int(from = base, to = cap, by = step)
  if (utils::tail(seq.out, 1L) < cap)
    seq.out <- c(seq.out, cap)
  unique(as.double(seq.out))
}

.np_glp_degree_hard_max <- 100L
.np_glp_degree_warn_threshold <- 25L

.np_warn_high_glp_degree <- function(value, argname = "degree") {
  value <- as.integer(value)
  if (!length(value))
    return(invisible(NULL))

  threshold <- .np_glp_degree_warn_threshold
  if (any(value > threshold)) {
    warning(
      sprintf(
        "%s contains unusually large polynomial degree values above %d; numerical instability is possible",
        argname,
        threshold
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

npValidateGlpDegree <- function(regtype, degree, ncon, argname = "degree") {
  degree.max <- .np_glp_degree_hard_max

  if (!identical(regtype, "lp"))
    return(NULL)

  if (ncon == 0L) {
    if (is.null(degree))
      stop(sprintf("%s must be 0 when regtype='lp' has no continuous predictors",
                   argname),
           call. = FALSE)
    if (!length(degree))
      return(integer(0))
    if (!is.numeric(degree) ||
        length(degree) != 1L ||
        anyNA(degree) ||
        any(!is.finite(degree)) ||
        degree != 0L)
      stop(sprintf("%s must be 0 when regtype='lp' has no continuous predictors",
                   argname),
           call. = FALSE)
    return(integer(0))
  }

  if (is.null(degree)) {
    stop(sprintf("%s must be supplied explicitly when regtype='lp'", argname))
  }

  if (length(degree) != ncon)
    stop(sprintf("%s must have one entry per continuous predictor (%d expected, got %d)",
                 argname, ncon, length(degree)))

  if (!is.numeric(degree))
    stop(sprintf("%s must contain finite non-negative integers", argname))

  if (anyNA(degree) || any(!is.finite(degree)))
    stop(sprintf("%s must contain finite non-negative integers", argname))

  if (any(degree < 0))
    stop(sprintf("%s must contain finite non-negative integers", argname))

  if (any(degree != floor(degree)))
    stop(sprintf("%s must contain finite non-negative integers", argname))

  if (any(degree > degree.max))
    stop(sprintf("%s must contain finite non-negative integers in [0,%d]",
                 argname, degree.max))

  .np_warn_high_glp_degree(degree, argname = argname)

  as.integer(degree)
}

npSetupGlpDegree <- function(regtype,
                             degree,
                             ncon,
                             degree.select = c("manual", "coordinate", "exhaustive")) {
  if (!identical(regtype, "lp"))
    return(degree)

  degree.select <- match.arg(degree.select, c("manual", "coordinate", "exhaustive"))
  if (!is.null(degree) || identical(degree.select, "manual"))
    return(degree)

  if (ncon == 0L)
    return(integer(0))

  rep.int(1L, ncon)
}

npValidateGlpBernstein <- function(regtype, bernstein.basis, argname = "bernstein.basis") {
  if (!identical(regtype, "lp"))
    return(FALSE)

  if (is.null(bernstein.basis))
    bernstein.basis <- FALSE

  if (!is.logical(bernstein.basis) || length(bernstein.basis) != 1L || is.na(bernstein.basis))
    stop(sprintf("%s must be TRUE or FALSE", argname))

  isTRUE(bernstein.basis)
}

npValidateScalarLogical <- function(value, argname) {
  value <- as.logical(value)
  if (length(value) != 1L || is.na(value))
    stop(sprintf("'%s' must be TRUE or FALSE", argname))
  value
}

npValidateNomadControl <- function(value, argname = "nomad") {
  if (is.logical(value))
    return(if (npValidateScalarLogical(value, argname)) "true" else "false")

  if (!is.character(value) || length(value) != 1L || is.na(value))
    stop(sprintf("'%s' must be TRUE, FALSE, or \"auto\"", argname), call. = FALSE)

  token <- tolower(trimws(value))
  if (token %in% c("true", "false", "auto"))
    return(token)

  stop(sprintf("'%s' must be TRUE, FALSE, or \"auto\"", argname), call. = FALSE)
}

npValidateNewdataColumns <- function(newdata, required, argname = "newdata") {
  nd <- toFrame(newdata)
  required <- unique(as.character(required))
  required <- required[nzchar(required)]
  missing.names <- setdiff(required, names(nd))
  if (length(missing.names) > 0L) {
    stop(sprintf("%s must contain columns: %s",
                 argname, paste(shQuote(required), collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

npValidateNewdataFormula <- function(newdata, formula, include.response = TRUE,
                                     argname = "newdata") {
  tt <- terms(formula)
  if (!isTRUE(include.response))
    tt <- delete.response(tt)
  npValidateNewdataColumns(newdata, all.vars(tt), argname = argname)
}

.np_prepare_nomad_shortcut <- function(nomad,
                                       call_names,
                                       preset,
                                       values,
                                       where = "npregbw") {
  call_names <- if (is.null(call_names)) character(0) else call_names[nzchar(call_names)]
  nomad.mode <- if ("nomad" %in% call_names) {
    npValidateNomadControl(nomad, "nomad")
  } else {
    "false"
  }
  nomad.enabled <- nomad.mode %in% c("true", "auto")

  metadata <- list(
    enabled = isTRUE(nomad.enabled),
    where = where,
    preset = "lp_nomad",
    source = if (identical(nomad.mode, "auto")) "auto" else if (identical(nomad.mode, "true")) "explicit" else "default",
    nomad = nomad.mode,
    auto.filled = character(0),
    user.supplied = character(0),
    normalized.values = preset
  )

  if (!isTRUE(nomad.enabled))
    return(list(enabled = FALSE, values = values, metadata = metadata))

  auto.filled <- character(0)
  user.supplied <- intersect(names(preset), call_names)

  for (arg in names(preset)) {
    if (!(arg %in% call_names)) {
      values[[arg]] <- preset[[arg]]
      auto.filled <- c(auto.filled, arg)
    }
  }

  metadata$auto.filled <- auto.filled
  metadata$user.supplied <- user.supplied
  metadata$normalized.values <- values[names(preset)]

  list(enabled = TRUE, values = values, metadata = metadata)
}

.np_attach_nomad_shortcut <- function(obj, metadata) {
  if (isTRUE(metadata$enabled))
    obj$nomad.shortcut <- metadata
  obj
}

npValidateNonNegativeInteger <- function(value, argname) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0 || value != floor(value))
    stop(sprintf("'%s' must be a non-negative integer", argname))
  as.integer(value)
}

npDefaultNmulti <- function(p) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) ||
      !is.finite(p) || p < 1 || p != floor(p))
    stop("'p' must be a positive integer")
  min(2L, as.integer(p))
}

npValidatePositiveInteger <- function(value, argname) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value))
    stop(sprintf("'%s' must be a positive integer", argname))
  as.integer(value)
}

npValidateNmulti <- function(value, argname = "nmulti") {
  npValidatePositiveInteger(value, argname)
}

npValidateBwsolver <- function(value, argname = "bwsolver") {
  match.arg(value, c("powell", "mads", "mads+powell"))
}

npBwsolverUsesMads <- function(value) {
  npValidateBwsolver(value) %in% c("mads", "mads+powell")
}

npRejectUnsupportedBwsolver <- function(dots, where) {
  if (!is.list(dots) || !("bwsolver" %in% names(dots)))
    return(invisible(FALSE))
  stop(sprintf(
    "bwsolver is not supported by %s; use that family's documented optimizer controls",
    where
  ), call. = FALSE)
}

npValidatePositiveFiniteNumeric <- function(value, argname) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value <= 0)
    stop(sprintf("'%s' must be a positive finite numeric scalar", argname))
  as.double(value)
}

npScaleFactorSearchOldNames <- c(
  "cfac.init",
  "lbc.init",
  "hbc.init",
  "scale.factor.lower.bound"
)

npScaleFactorSearchNewNames <- c(
  cfac.init = "scale.factor.init",
  lbc.init = "scale.factor.init.lower",
  hbc.init = "scale.factor.init.upper",
  scale.factor.lower.bound = "scale.factor.search.lower"
)

npRejectRenamedScaleFactorSearchArgs <- function(arg.names,
                                                 where = "bandwidth search") {
  old <- intersect(npScaleFactorSearchOldNames, arg.names)
  if (!length(old))
    return(invisible(NULL))

  messages <- sprintf("'%s' has been renamed to '%s'",
                      old,
                      unname(npScaleFactorSearchNewNames[old]))
  stop(sprintf("%s: %s", where, paste(messages, collapse = "; ")),
       call. = FALSE)
}

npRejectUnsupportedLpDegreeSearchArgs <- function(arg.names,
                                                  where = "bandwidth search") {
  if (is.null(arg.names) || !length(arg.names))
    return(invisible(NULL))

  unsupported <- c(
    "nomad", "nomad.nmulti", "nomad.remin",
    "search.engine",
    "degree", "degree.select", "degree.min", "degree.max", "degree.start",
    "degree.restarts", "degree.max.cycles", "degree.verify",
    "random.seed",
    "regtype", "basis", "bernstein.basis"
  )
  bad <- intersect(unsupported, arg.names[nzchar(arg.names)])
  if (!length(bad))
    return(invisible(NULL))

  labels <- paste(sprintf("'%s'", bad), collapse = ", ")
  stop(sprintf("unused argument%s in %s: %s",
               if (length(bad) == 1L) "" else "s",
               where,
               labels),
       call. = FALSE)
}

npValidateScaleFactorLowerBound <- function(value,
                                            argname = "scale.factor.search.lower") {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0) {
    stop(sprintf("'%s' must be a nonnegative finite numeric scalar", argname),
         call. = FALSE)
  }

  as.double(value)
}

npResolveScaleFactorLowerBound <- function(value,
                                           fallback = 0.1,
                                           argname = "scale.factor.search.lower") {
  if (is.null(value))
    return(as.double(fallback))

  npValidateScaleFactorLowerBound(value, argname = argname)
}

npResolveScaleFactorSearchLower <- function(value,
                                            fallback = 0.1,
                                            argname = "scale.factor.search.lower") {
  npResolveScaleFactorLowerBound(value, fallback = fallback, argname = argname)
}

npGetScaleFactorSearchLower <- function(object,
                                        fallback = 0.1,
                                        argname = "scale.factor.search.lower") {
  new <- object[["scale.factor.search.lower"]]
  old <- object[["scale.factor.lower.bound"]]

  if (!is.null(new) && !is.null(old)) {
    new.value <- npResolveScaleFactorSearchLower(new, fallback = fallback,
                                                 argname = "scale.factor.search.lower")
    old.value <- npResolveScaleFactorSearchLower(old, fallback = fallback,
                                                 argname = "scale.factor.lower.bound")
    if (!identical(new.value, old.value)) {
      stop(sprintf(
        "%s conflicts with legacy 'scale.factor.lower.bound' (%.15g versus %.15g)",
        argname,
        new.value,
        old.value
      ), call. = FALSE)
    }
    return(new.value)
  }

  if (!is.null(new))
    return(npResolveScaleFactorSearchLower(new, fallback = fallback,
                                           argname = argname))

  npResolveScaleFactorSearchLower(old, fallback = fallback,
                                  argname = "scale.factor.lower.bound")
}

npSetScaleFactorSearchLower <- function(object, value, fallback = 0.1) {
  object[["scale.factor.search.lower"]] <-
    npResolveScaleFactorSearchLower(value, fallback = fallback)
  object[["scale.factor.lower.bound"]] <- NULL
  object
}

npEffectiveContinuousStartLower <- function(scale.factor.init.lower,
                                            scale.factor.search.lower) {
  max(as.double(scale.factor.init.lower), as.double(scale.factor.search.lower))
}

npEffectiveContinuousStartPoint <- function(scale.factor.init,
                                            scale.factor.search.lower) {
  max(as.double(scale.factor.init), as.double(scale.factor.search.lower))
}

npContinuousSearchStartControls <- function(scale.factor.init.lower,
                                            scale.factor.init.upper,
                                            scale.factor.init,
                                            scale.factor.search.lower,
                                            where = "bandwidth search") {
  lbc.raw <- npValidatePositiveFiniteNumeric(scale.factor.init.lower,
                                             "scale.factor.init.lower")
  scale.factor.init.upper <- npValidatePositiveFiniteNumeric(scale.factor.init.upper,
                                             "scale.factor.init.upper")
  cfac.raw <- npValidatePositiveFiniteNumeric(scale.factor.init,
                                             "scale.factor.init")
  scale.factor.search.lower <- npValidateScaleFactorLowerBound(
    scale.factor.search.lower,
    "scale.factor.search.lower"
  )

  lbc.effective <- npEffectiveContinuousStartLower(
    lbc.raw,
    scale.factor.search.lower
  )
  cfac.effective <- npEffectiveContinuousStartPoint(
    cfac.raw,
    scale.factor.search.lower
  )

  if (scale.factor.init.upper < lbc.effective) {
    stop(sprintf(
      "%s: 'scale.factor.init.upper' must be greater than or equal to max('scale.factor.init.lower', 'scale.factor.search.lower') (effective lower %.15g; scale.factor.init.upper %.15g)",
      where,
      lbc.effective,
      scale.factor.init.upper
    ), call. = FALSE)
  }

  list(
    scale.factor.init.lower = lbc.effective,
    scale.factor.init.upper = scale.factor.init.upper,
    scale.factor.init = cfac.effective
  )
}

npValidateLpBasis <- function(regtype, basis, argname = "basis") {
  if (!identical(regtype, "lp"))
    return("glp")

  if (is.null(basis))
    basis <- "glp"

  basis <- match.arg(as.character(basis), c("glp", "additive", "tensor"))
  basis
}

npLpBasisCode <- function(basis) {
  basis.value <- if (is.null(basis) || !length(basis)) "glp" else basis
  switch(tolower(basis.value),
         additive = 0L,
         glp = 1L,
         tensor = 2L,
         1L)
}

npValidateGlpGradientOrder <- function(regtype,
                                       gradient.order,
                                       ncon,
                                       argname = "gradient.order") {
  degree.max <- 12L

  if (!identical(regtype, "lp"))
    return(NULL)

  if (ncon == 0L)
    return(npValidateCategoricalFirstDifferenceGradientOrder(
      gradient.order = gradient.order,
      argname = argname
    ))

  if (is.null(gradient.order))
    gradient.order <- rep.int(1L, ncon)

  if (!length(gradient.order) && ncon == 0L)
    return(integer(0))

  if (length(gradient.order) == 1L && ncon > 1L)
    gradient.order <- rep.int(gradient.order, ncon)

  if (length(gradient.order) != ncon)
    stop(sprintf("%s must have one entry per continuous predictor (%d expected, got %d)",
                 argname, ncon, length(gradient.order)))

  if (!is.numeric(gradient.order))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (anyNA(gradient.order) || any(!is.finite(gradient.order)))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(gradient.order <= 0))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(gradient.order != floor(gradient.order)))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(gradient.order > degree.max))
    stop(sprintf("%s must contain integers in [1,%d]", argname, degree.max))

  as.integer(gradient.order)
}

npValidateCategoricalFirstDifferenceGradientOrder <- function(
    gradient.order,
    argname = "gradient.order",
    where = "categorical gradients/effects") {
  if (is.null(gradient.order) || !length(gradient.order))
    return(integer(0))

  if (!is.numeric(gradient.order))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (anyNA(gradient.order) || any(!is.finite(gradient.order)))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(gradient.order <= 0))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(gradient.order != floor(gradient.order)))
    stop(sprintf("%s must contain finite positive integers", argname))

  if (any(as.integer(gradient.order) != 1L)) {
    stop(sprintf("%s support only first differences for ordered/unordered predictors",
                 where),
         call. = FALSE)
  }

  integer(0)
}

.np_singleindex_reject_higher_gradient_order <- function(
    dots,
    where = "npindex") {
  order.names <- intersect(c("gradient.order", "gradient_order"), names(dots))
  for (argname in order.names) {
    value <- dots[[argname]]
    if (is.null(value))
      next
    if (!is.numeric(value) || !length(value) || anyNA(value) ||
        any(!is.finite(value)) || any(value <= 0) ||
        any(value != floor(value))) {
      stop(sprintf("%s must contain finite positive integers", argname),
           call. = FALSE)
    }
    if (any(as.integer(value) != 1L)) {
      stop(sprintf(
        "%s supports only first-order gradients; %s > 1 is not supported",
        where,
        argname
      ), call. = FALSE)
    }
  }
  invisible(NULL)
}

npValidateLcGradientOrder <- function(regtype,
                                      gradient.order,
                                      ncon,
                                      argname = "gradient.order",
                                      where = "local-constant regression") {
  if (!identical(regtype, "lc"))
    return(invisible(NULL))

  gradient.order <- npValidateGlpGradientOrder(
    regtype = "lp",
    gradient.order = gradient.order,
    ncon = ncon,
    argname = argname
  )
  if (any(gradient.order != 1L)) {
    stop(sprintf("%s supports only first derivatives for regtype='lc'; use regtype='lp' with sufficient degree for higher-order derivatives",
                 where),
         call. = FALSE)
  }

  gradient.order
}

npGlpDegree0FirstDerivativeLcOk <- function(regtype.engine,
                                            degree.engine,
                                            gradient.order,
                                            ncon) {
  identical(regtype.engine, "lp") &&
    ncon > 0L &&
    length(degree.engine) == ncon &&
    all(degree.engine == 0L) &&
    length(gradient.order) == ncon &&
    all(gradient.order == 1L)
}

npGlpGradientAvailability <- function(regtype.engine,
                                      degree.engine,
                                      gradient.order,
                                      ncon,
                                      allow.lp0.lc.first = TRUE,
                                      where = "local-polynomial estimator") {
  if (!identical(regtype.engine, "lp") || ncon == 0L)
    return(rep.int(TRUE, ncon))

  degree.engine <- as.integer(degree.engine)
  gradient.order <- as.integer(gradient.order)

  if (length(degree.engine) != ncon || length(gradient.order) != ncon) {
    stop(sprintf("%s received incoherent local-polynomial derivative metadata", where),
         call. = FALSE)
  }

  available <- gradient.order <= degree.engine
  if (isTRUE(allow.lp0.lc.first) &&
      npGlpDegree0FirstDerivativeLcOk(
        regtype.engine = regtype.engine,
        degree.engine = degree.engine,
        gradient.order = gradient.order,
        ncon = ncon
      )) {
    available[] <- TRUE
  }

  available
}

npGlpGradientUnavailableSummary <- function(degree.engine,
                                            gradient.order,
                                            available,
                                            con.names = NULL) {
  bad <- which(!available)
  if (!length(bad))
    return(character(0))

  if (is.null(con.names) || length(con.names) != length(available)) {
    con.names <- paste0("continuous[", seq_along(available), "]")
  }

  sprintf("%s (requested order %d, degree %d)",
          con.names[bad],
          as.integer(gradient.order)[bad],
          as.integer(degree.engine)[bad])
}

npStopGlpGradientNoneAvailable <- function(where,
                                           action = "provide",
                                           degree.engine,
                                           gradient.order,
                                           available,
                                           con.names = NULL) {
  if (!length(available) || any(available)) {
    stop(sprintf("%s received an invalid all-unavailable derivative condition", where),
         call. = FALSE)
  }

  unavailable <- npGlpGradientUnavailableSummary(
    degree.engine = degree.engine,
    gradient.order = gradient.order,
    available = available,
    con.names = con.names
  )
  if (!length(unavailable)) {
    stop(sprintf("%s received incoherent all-unavailable derivative metadata", where),
         call. = FALSE)
  }

  stop(sprintf(
    paste0(
      "%s cannot %s the requested derivatives because no requested component ",
      "is available (no available derivative components) at the fitted polynomial ",
      "degrees: %s. Lower gradient.order ",
      "for at least one continuous predictor or refit with sufficient polynomial degree."
    ),
    where,
    action,
    paste(unavailable, collapse = "; ")
  ), call. = FALSE)
}

npCategoricalFirstDifferenceFrames <- function(exdat, index, where) {
  exdat <- toFrame(exdat)
  index <- as.integer(index)[1L]
  if (is.na(index) || index < 1L || index > ncol(exdat)) {
    stop(sprintf("%s received an invalid categorical coordinate", where),
         call. = FALSE)
  }

  x <- exdat[[index]]
  if (!is.factor(x)) {
    stop(sprintf("%s received a non-categorical coordinate", where),
         call. = FALSE)
  }

  lev <- levels(x)
  if (!length(lev)) {
    stop(sprintf("%s received a categorical coordinate without levels", where),
         call. = FALSE)
  }

  code <- as.integer(x)
  lower.code <- rep.int(1L, length(code))
  upper.code <- code
  if (is.ordered(x)) {
    if (length(lev) < 2L) {
      lower.code[] <- 1L
      upper.code[] <- 1L
    } else {
      lower.code <- pmax.int(code - 1L, 1L)
      first <- !is.na(code) & code == 1L
      upper.code[first] <- 2L
    }
  }

  lower <- upper <- exdat
  lower[[index]] <- factor(
    lev[lower.code], levels = lev, ordered = is.ordered(x)
  )
  upper[[index]] <- factor(
    lev[upper.code], levels = lev, ordered = is.ordered(x)
  )

  list(lower = lower, upper = upper)
}

npWarnGlpGradientPartialAvailability <- function(where,
                                                 degree.engine,
                                                 gradient.order,
                                                 available,
                                                 con.names = NULL) {
  unavailable <- npGlpGradientUnavailableSummary(
    degree.engine = degree.engine,
    gradient.order = gradient.order,
    available = available,
    con.names = con.names
  )
  if (!length(unavailable))
    return(invisible(NULL))

  .np_warning(
    where,
    ": requested derivative order would exceed polynomial degree or is unavailable for ",
    paste(unavailable, collapse = ", "),
    "; returning NA for unavailable gradient component(s)",
    call. = FALSE
  )

  invisible(NULL)
}

npValidateGlpGradientDegree <- function(regtype.engine,
                                        degree.engine,
                                        gradient.order,
                                        ncon,
                                        where = "local-polynomial estimator",
                                        allow.lp0.lc.first = TRUE) {
  if (!identical(regtype.engine, "lp") || ncon == 0L)
    return(invisible(gradient.order))

  degree.engine <- as.integer(degree.engine)
  gradient.order <- as.integer(gradient.order)

  if (!length(degree.engine) || !length(gradient.order))
    return(invisible(gradient.order))

  if (length(degree.engine) != ncon || length(gradient.order) != ncon) {
    stop(sprintf("%s received incoherent local-polynomial derivative metadata", where),
         call. = FALSE)
  }

  if (isTRUE(allow.lp0.lc.first) &&
      npGlpDegree0FirstDerivativeLcOk(
        regtype.engine = regtype.engine,
        degree.engine = degree.engine,
        gradient.order = gradient.order,
        ncon = ncon
      ))
    return(invisible(gradient.order))

  bad <- which(gradient.order > degree.engine)
  if (length(bad)) {
    stop(sprintf(
      "%s supports derivative orders only up to the fitted polynomial degree; requested order %s for degree %s",
      where,
      paste(gradient.order[bad], collapse = ","),
      paste(degree.engine[bad], collapse = ",")
    ), call. = FALSE)
  }

  invisible(gradient.order)
}

npConditionalRegEngineSpec <- function(bws, where = "conditional estimator") {
  reg.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  basis.engine <- if (is.null(bws$basis.engine)) {
    if (is.null(bws$basis)) "glp" else bws$basis
  } else {
    bws$basis.engine
  }
  degree.engine <- if (is.null(bws$degree.engine)) {
    if (bws$xncon > 0L) {
      if (identical(reg.engine, "lc")) {
        rep.int(0L, bws$xncon)
      } else {
        npValidateGlpDegree(
          regtype = "lp",
          degree = bws$degree,
          ncon = bws$xncon
        )
      }
    } else {
      integer(0)
    }
  } else {
    as.integer(bws$degree.engine)
  }
  bernstein.engine <- if (is.null(bws$bernstein.basis.engine)) {
    isTRUE(bws$bernstein.basis)
  } else {
    isTRUE(bws$bernstein.basis.engine)
  }

  if (!identical(reg.engine, "lp") && bws$xncon > 0L)
    degree.engine <- rep.int(if (identical(reg.engine, "ll")) 1L else 0L, bws$xncon)

  list(
    reg.engine = reg.engine,
    basis.engine = basis.engine,
    degree.engine = degree.engine,
    bernstein.engine = bernstein.engine
  )
}

npConditionalGradientOrder <- function(bws,
                                       reg.engine,
                                       gradient.order,
                                       where = "conditional estimator") {
  if (bws$xncon == 0L)
    return(integer(0))

  gorder <- if (is.null(gradient.order)) {
    rep.int(1L, bws$xncon)
  } else {
    npValidateGlpGradientOrder(
      regtype = "lp",
      gradient.order = gradient.order,
      ncon = bws$xncon
    )
  }

  if (!identical(reg.engine, "lp") && any(gorder > 1L)) {
    stop(sprintf(
      "%s supports gradient.order > 1 only for regtype='lp' continuous explanatory predictors",
      where
    ), call. = FALSE)
  }

  gorder
}

npCheckRegressionDesignCondition <- function(reg.code,
                                             xcon,
                                             basis = "glp",
                                             degree = NULL,
                                             bernstein.basis = FALSE,
                                             where = "npregbw") {
  kappa.warn <- 1e8
  kappa.stop <- 1e12

  if (!(reg.code %in% c(REGTYPE_LL, REGTYPE_LP)))
    return(invisible(NULL))

  xcon <- as.data.frame(xcon)
  n <- nrow(xcon)
  if (is.null(n) || n <= 0L)
    return(invisible(NULL))

  B <- if (identical(reg.code, REGTYPE_LP)) {
    if (is.null(degree))
      stop(sprintf("%s: LP degree vector missing for design-conditioning check", where))
    W.lp(xdat = xcon,
         degree = degree,
         basis = basis,
         bernstein.basis = isTRUE(bernstein.basis))
  } else {
    cbind(1, as.matrix(xcon))
  }

  p <- ncol(B)
  if (is.null(p) || p <= 0L)
    return(invisible(NULL))

  sv <- suppressWarnings(tryCatch(svd(B, nu = 0L, nv = 0L)$d, error = function(e) NULL))
  if (is.null(sv) || !length(sv))
    stop(sprintf("%s: unable to compute singular values for design-conditioning check", where))
  tol.rank <- max(dim(B)) * max(sv) * .Machine$double.eps
  r <- sum(sv > tol.rank)
  if (r < p) {
    stop(sprintf("%s: regression design matrix is rank deficient (rank=%d < p=%d). Reduce polynomial degree or remove collinear continuous predictors.",
                 where, r, p))
  }

  kB <- suppressWarnings(tryCatch(kappa(B), error = function(e) Inf))
  if (!is.finite(kB))
    kB <- Inf

  if (kB > kappa.stop) {
    stop(sprintf("%s: regression design matrix is severely ill-conditioned (kappa(B)=%.3e > %.1e). Reduce polynomial degree or remove collinear continuous predictors.",
                 where, kB, kappa.stop))
  }

  if (kB > kappa.warn) {
    .np_warning(sprintf("%s: regression design matrix is ill-conditioned (kappa(B)=%.3e > %.1e). Estimation may rely heavily on ridging; consider lower degree or less collinear predictors.",
                    where, kB, kappa.warn),
            call. = FALSE, immediate. = TRUE)
  }

  invisible(NULL)
}

npRegtypeToC <- function(regtype, degree, ncon, context = "npreg") {
  # Internal regression-type codes:
  # lc -> REGTYPE_LC, ll -> REGTYPE_LL, lp -> REGTYPE_LP.
  if (identical(regtype, "lc"))
    return(list(code = REGTYPE_LC, degree = NULL))

  if (identical(regtype, "ll"))
    return(list(code = REGTYPE_LL, degree = NULL))

  degree <- npValidateGlpDegree(regtype, degree, ncon)

  list(code = REGTYPE_LP, degree = degree)
}

npCanonicalConditionalRegSpec <- function(regtype = c("lc", "ll", "lp"),
                                          basis = c("glp", "additive", "tensor"),
                                          degree = NULL,
                                          bernstein.basis = FALSE,
                                          ncon,
                                          where = "npc*") {
  regtype <- match.arg(regtype)
  ncon <- npValidateNonNegativeInteger(ncon, "ncon")
  basis <- npValidateLpBasis(regtype = "lp", basis = basis)
  bernstein.basis <- npValidateGlpBernstein(regtype = "lp",
                                            bernstein.basis = bernstein.basis)

  if (identical(regtype, "lc")) {
    degree <- rep.int(0L, ncon)
    return(list(
      regtype = "lc",
      basis = basis,
      degree = degree,
      bernstein.basis = FALSE,
      regtype.engine = "lc",
      basis.engine = basis,
      degree.engine = degree,
      bernstein.basis.engine = FALSE
    ))
  }

  if (identical(regtype, "ll")) {
    if (ncon == 0L)
      stop(sprintf("%s: regtype='ll' requires at least one continuous predictor; use regtype='lc' for categorical-only predictors",
                   where),
           call. = FALSE)

    if (ncon > 0L) {
      degree <- rep.int(1L, ncon)
      basis <- "glp"
    } else {
      degree <- integer(0)
      basis <- "glp"
    }

    return(list(
      regtype = "ll",
      basis = basis,
      degree = degree,
      bernstein.basis = FALSE,
      regtype.engine = if (ncon > 0L) "lp" else "lc",
      basis.engine = basis,
      degree.engine = degree,
      bernstein.basis.engine = FALSE
    ))
  }

  degree <- npValidateGlpDegree(regtype = "lp",
                                degree = degree,
                                ncon = ncon)
  if (ncon == 0L && length(degree) == 0L)
    bernstein.basis <- FALSE
  list(
    regtype = "lp",
    basis = basis,
    degree = degree,
    bernstein.basis = bernstein.basis,
    regtype.engine = if (ncon > 0L) "lp" else "lc",
    basis.engine = basis,
    degree.engine = degree,
    bernstein.basis.engine = bernstein.basis
  )
}

npWithLocalLinearRawBasisSearchError <- function(expr,
                                                 where,
                                                 spec,
                                                 bwmethod,
                                                 ncon) {
  expected <- switch(where,
    npcdensbw = "C_np_density_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective",
    npcdistbw = "C_np_distribution_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective",
    NULL
  )
  tryCatch(
    force(expr),
    error = function(e) {
      msg <- conditionMessage(e)
      degree <- if (is.null(spec$degree.engine)) integer(0) else as.integer(spec$degree.engine)
      ncon <- as.integer(ncon)
      targeted <- !is.null(expected) &&
        identical(msg, expected) &&
        identical(spec$regtype, "ll") &&
        identical(bwmethod, "cv.ls") &&
        identical(spec$regtype.engine, "lp") &&
        !isTRUE(spec$bernstein.basis.engine) &&
        length(degree) == ncon &&
        ncon > 0L &&
        all(degree == 1L)

      if (targeted) {
        stop(sprintf(
          "%s() local-linear cv.ls failed while using the canonical raw degree-1 basis. Try regtype = \"lp\", degree = 1, bernstein.basis = TRUE, or center/scale continuous regressors.",
          where
        ), call. = FALSE)
      }
      stop(e)
    }
  )
}

npResolveCanonicalConditionalRegSpec <- function(mc.names,
                                                 regtype = c("lc", "ll", "lp"),
                                                 basis = c("glp", "additive", "tensor"),
                                                 degree = NULL,
                                                 bernstein.basis = FALSE,
                                                 ncon,
                                                 where = "npreg") {
  mc.names <- if (is.null(mc.names)) character(0) else as.character(mc.names)
  regtype.named <- any(mc.names == "regtype")
  basis.named <- any(mc.names == "basis")
  degree.named <- any(mc.names == "degree")
  bernstein.named <- any(mc.names == "bernstein.basis")

  regtype <- if (regtype.named) match.arg(regtype) else "lc"

  if (identical(regtype, "lc") && (basis.named || degree.named || bernstein.named)) {
    stop("regtype='lc' does not accept basis/degree/bernstein.basis; use regtype='lp' for local-polynomial controls")
  }

  if (identical(regtype, "ll")) {
    if (degree.named) {
      degree <- npValidateGlpDegree(regtype = "lp", degree = degree, ncon = ncon)
      if (!identical(as.integer(degree), rep.int(1L, ncon)))
        stop("regtype='ll' uses canonical LP(degree=1, basis='glp'); remove 'degree' or use regtype='lp'")
    }
    if (basis.named) {
      basis <- match.arg(basis)
      if (!identical(basis, "glp"))
        stop("regtype='ll' uses canonical basis='glp'; use regtype='lp' for alternate LP bases")
    }
    if (bernstein.named) {
      bernstein.basis <- npValidateGlpBernstein(regtype = "lp",
                                                bernstein.basis = bernstein.basis)
      if (isTRUE(bernstein.basis))
        stop("regtype='ll' uses canonical bernstein.basis=FALSE; use regtype='lp' for Bernstein LP")
    }
  }

  npCanonicalConditionalRegSpec(
    regtype = regtype,
    basis = basis,
    degree = degree,
    bernstein.basis = bernstein.basis,
    ncon = ncon,
    where = where
  )
}

npRejectLegacyLpArgs <- function(dotnames, where = "npreg") {
  if (is.null(dotnames) || !length(dotnames))
    return(invisible(NULL))
  bad <- intersect(dotnames, c("glp.degree", "glp.bernstein", "glp.basis"))
  if (length(bad))
    stop(sprintf("%s: legacy arguments %s are no longer supported; use basis, degree and bernstein.basis",
                 where, paste(sprintf("'%s'", bad), collapse = ", ")))
  if ("remin" %in% dotnames)
    stop(sprintf("%s: argument 'remin' has been replaced by 'powell.remin' and 'nomad.remin'",
                 where))
  invisible(NULL)
}

npDimBasisCapacityError <- function() {
  stop("dim_basis: basis dimension exceeds supported capacity", call. = FALSE)
}

npCheckDimBasisNativeCapacity <- function(degree, segments, include, categories, kernel) {
  limit <- .Machine$integer.max

  if ((length(degree) && any(degree > limit)) ||
      (length(segments) && any(segments > limit)))
    npDimBasisCapacityError()

  if (length(degree)) {
    rows <- as.double(degree) + as.double(segments)
    if (any(!is.finite(rows)) || any(rows > limit))
      npDimBasisCapacityError()
  }

  if (!isTRUE(kernel) &&
      ((length(include) && any(include > limit)) ||
       (length(categories) && any(categories > limit))))
    npDimBasisCapacityError()

  if (!isTRUE(kernel) && length(include)) {
    categorical <- as.double(include) * as.double(categories)
    if (any(!is.finite(categorical)) || any(categorical > limit))
      npDimBasisCapacityError()
  }

  invisible(NULL)
}

dim_basis <- function(basis = c("glp", "additive", "tensor"),
                      kernel = TRUE,
                      degree = NULL,
                      segments = NULL,
                      include = NULL,
                      categories = NULL) {
  basis <- match.arg(basis)

  if (is.null(degree))
    degree <- integer(0)
  if (is.null(segments))
    segments <- rep.int(1L, length(degree))

  if (length(degree) != length(segments))
    stop("degree and segments must have the same length")
  if (!is.numeric(degree) ||
      anyNA(degree) ||
      any(!is.finite(degree)) ||
      any(degree < 0L) ||
      any(degree != floor(degree)))
    stop("degree must contain finite non-negative integers")
  if (!is.numeric(segments) ||
      anyNA(segments) ||
      any(!is.finite(segments)) ||
      any(segments <= 0L) ||
      any(segments != floor(segments)))
    stop("segments must contain finite positive integers")

  degree <- as.integer(degree)
  segments <- as.integer(segments)

  if (is.null(include))
    include <- integer(0)
  if (is.null(categories))
    categories <- integer(0)

  if (length(include) != length(categories))
    stop("include and categories must have the same length")
  if (!is.numeric(include) ||
      anyNA(include) ||
      any(!is.finite(include)) ||
      any(include < 0L) ||
      any(include != floor(include)))
    stop("include must contain finite non-negative integers")
  if (!is.numeric(categories) ||
      anyNA(categories) ||
      any(!is.finite(categories)) ||
      any(categories < 0L) ||
      any(categories != floor(categories)))
    stop("categories must contain finite non-negative integers")

  npCheckDimBasisNativeCapacity(degree, segments, include, categories, kernel)

  include <- as.integer(include)
  categories <- as.integer(categories)

  basis.code <- switch(basis, additive = 0L, glp = 1L, tensor = 2L)

  result <- .Call("C_np_dim_basis",
                  as.integer(basis.code),
                  as.integer(isTRUE(kernel)),
                  degree,
                  segments,
                  include,
                  categories,
                  PACKAGE = "np")

  if (!is.finite(result) || result > .Machine$integer.max)
    npDimBasisCapacityError()

  result
}

dimBS <- function(basis = "additive",
                  kernel = TRUE,
                  degree = NULL,
                  segments = NULL,
                  include = NULL,
                  categories = NULL) {
  dim_basis(basis = basis,
            kernel = kernel,
            degree = degree,
            segments = segments,
            include = include,
            categories = categories)
}

## Function to test for monotone increasing vector

is.monotone.increasing <- function(x) {
  ## Sorted and last value > first value
  !is.unsorted(x) && x[length(x)] > x[1]
}
  
## what is a badord? ... an ordered factor of numeric values to treat
## them properly one must preserve the numeric value, ie. scale not
## just their sorted order Actually, the ord/badord paradigm must go,
## in place of levels caching

matrix.sd <- function(x, na.rm=FALSE) {
  if(is.matrix(x)) apply(x, 2, sd, na.rm=na.rm)
  else if(is.vector(x)) sd(x, na.rm=na.rm)
  else if(is.data.frame(x)) sapply(x, sd, na.rm=na.rm)
  else sd(as.vector(x), na.rm=na.rm)
}

.np_validate_seed_scalar <- function(seed) {
  if (is.null(seed) || length(seed) != 1L || !is.numeric(seed) || is.na(seed) || !is.finite(seed))
    stop("'seed' must be a single finite numeric value", call. = FALSE)

  normalized <- abs(as.double(seed)[1L])

  if (normalized > .Machine$integer.max || normalized != floor(normalized))
    stop("'seed' must be representable as a non-negative integer after abs()", call. = FALSE)

  as.integer(normalized)
}

npseed <- function(seed){
  .Call("C_np_set_seed", .np_validate_seed_scalar(seed), PACKAGE = "np")
  invisible()
}

erf <- function(z) { 2 * pnorm(z*sqrt(2)) - 1 }

nptgauss <- function(b){

  rel.tol <- sqrt(.Machine$double.eps)

  b.max <- sqrt(-2*log(.Machine$double.eps))

  if((b < 0) || (b > b.max))
    stop(paste("b must be between 0 and",b.max))

  alpha <- 1.0/(pnorm(b)-pnorm(-b)-2*b*dnorm(b))

  tgauss <- function(z) {
    out <- alpha * (dnorm(z) - dnorm(b))
    out[abs(z) >= b] <- 0.0
    out
  }

  c0 <- alpha*dnorm(b)

  k <- integrate(f = function(z) { tgauss(z)^2 }, -b, b)$value
  k2 <- integrate(f = function(z) { z^2*tgauss(z) }, -b, b)$value
  k22 <- integrate(f = function(z) { (z*tgauss(z))^2 }, -b, b)$value
  km <- integrate(f = function(z) { tgauss(z-0.5)*tgauss(z+0.5) }, -b+0.5, b-0.5)$value

  a0 <- (0.5 + 2*b*c0)/integrate(f = function(z){ erf(z/2 + b)*exp(-0.25*z^2) }, -2*b, 0)$value
  a2 <- (c0 + k - a0*erf(b))/erf(b/sqrt(2))
  a1 <- -(a2*erf(b/sqrt(2)) + c0)/(2*b)

  int.kernels[CKER_TGAUSS + 1] <- k
  
  invisible(.Call("C_np_set_tgauss2",
                  as.double(c(b, alpha, c0, a0, a1, a2, k, k2, k22, km)),
                  PACKAGE = "np"))

}

numNotIn <- function(x){
  x <- unique(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x) || !(0 %in% x))
    return(0)

  n.unique <- length(x)
  if (n.unique >= .Machine$integer.max)
    stop("too many category codes to construct a deterministic padding value", call. = FALSE)

  candidate.max <- as.integer(n.unique + 1)
  candidates <- c(seq_len(candidate.max), -seq_len(candidate.max))
  # The candidate set is larger than the finite unique set, so one value is absent.
  candidates[match(FALSE, candidates %in% x)]
}

dlev <- function(x){
  if(is.ordered(x))
    x.dlev <- suppressWarnings(as.numeric(levels(x)))
  if (!is.ordered(x) || any(is.na(x.dlev)))
    x.dlev <- as.numeric(seq_len(nlevels(x)))
  x.dlev
}

isNum <- function(x){
  return(!any(is.na(suppressWarnings(as.numeric(x)))))
}

untangle <- function(frame){
  if (is.null(frame))
    return(NULL)
  
  iord <- unlist(lapply(frame, is.ordered))
  iuno <- unlist(lapply(frame, is.factor)) & !iord
  icon <- unlist(lapply(frame, is.numeric))

  if(!all(iord | iuno | icon)) 
    stop("non-allowed data type in data frame")

  inumord <-
    suppressWarnings(unlist(lapply(frame,
    function (z) {
      is.ordered(z) && is.numeric(tryCatch(as.numeric(levels(z)), warning = function (y) {
        FALSE }))
    })))

  all.lev <- lapply(frame, function(y){
    t.ret <- NULL
    if (is.factor(y))
      t.ret <- levels(y)
    t.ret
  })

  all.ulev <- lapply(frame, function (y) {
    t.ret <- NULL
    if (is.factor(y))
      t.ret <- sort(unique(y))
    t.ret
  })

  all.dlev <- lapply(frame, function (y) {
    t.ret <- NULL
    if (is.factor(y))
      t.ret <- dlev(y)
    t.ret
  })

  all.nlev <- lapply(frame, function(y) {
    t.ret <- NULL
    if (is.factor(y))
      t.ret <- nlevels(y)
    t.ret
  })

  all.min <- lapply(frame, function(y) {
    t.ret <- NULL
    if (is.numeric(y))
      t.ret <- min(y, na.rm = TRUE)
    t.ret
  })

  all.max <- lapply(frame, function(y) {
    t.ret <- NULL
    if (is.numeric(y))
      t.ret <- max(y, na.rm = TRUE)
    t.ret
  })

  list(iord = iord,
       iuno = iuno,
       icon = icon,
       inumord = inumord,
       all.lev = all.lev,
       all.ulev = all.ulev,
       all.dlev = all.dlev,
       all.nlev = all.nlev,
       all.min = all.min,
       all.max = all.max)
}

npKernelBoundsResolve <- function(dati,
                                  varnames,
                                  kerbound = c("none", "range", "fixed"),
                                  kerlb = NULL,
                                  kerub = NULL,
                                  argprefix = "cker") {
  kerbound <- match.arg(kerbound)
  icon.idx <- which(dati$icon)
  ncon <- length(icon.idx)
  if (is.null(varnames) || length(varnames) != length(dati$icon))
    varnames <- paste0("V", seq_along(dati$icon))

  full.lb <- rep(-Inf, length(dati$icon))
  full.ub <- rep(Inf, length(dati$icon))

  if (ncon == 0L) {
    return(list(bound = "none", lb = full.lb, ub = full.ub))
  }

  mins <- unlist(dati$all.min[icon.idx], use.names = FALSE)
  maxs <- unlist(dati$all.max[icon.idx], use.names = FALSE)
  cnames <- varnames[icon.idx]

  recycleBounds <- function(x, nm) {
    if (is.null(x))
      stop(sprintf("'%s' requires '%s' and '%s'.",
                   paste0(argprefix, "bound"),
                   paste0(argprefix, "lb"),
                   paste0(argprefix, "ub")))
    x <- as.numeric(x)
    if (!(length(x) %in% c(1L, ncon))) {
      stop(sprintf("length(%s) must be 1 or equal to the number of continuous variables (%d).",
                   nm, ncon))
    }
    if (length(x) == 1L)
      rep(x, ncon)
    else
      x
  }

  if (kerbound == "none") {
    lb <- rep(-Inf, ncon)
    ub <- rep(Inf, ncon)
  } else if (kerbound == "range") {
    lb <- mins
    ub <- maxs
  } else {
    lb <- recycleBounds(kerlb, paste0(argprefix, "lb"))
    ub <- recycleBounds(kerub, paste0(argprefix, "ub"))
  }

  if (any(!is.finite(lb) & !is.infinite(lb)) || any(!is.finite(ub) & !is.infinite(ub))) {
    stop(sprintf("'%s' and '%s' must be finite values or +/-Inf.",
                 paste0(argprefix, "lb"),
                 paste0(argprefix, "ub")))
  }

  if (any(lb >= ub)) {
    stop(sprintf("Invalid bounds for '%s': require lower < upper for every continuous variable.",
                 paste0(argprefix, "bound")))
  }

  bad.cover <- (lb > mins) | (ub < maxs)
  if (any(bad.cover)) {
    bad.vars <- paste(cnames[bad.cover], collapse = ", ")
    stop(sprintf("Invalid bounds for '%s': require lower <= min(training) and upper >= max(training) for each continuous variable. Violations: %s",
                 paste0(argprefix, "bound"), bad.vars))
  }

  if (kerbound == "fixed" && all(is.infinite(lb) & (lb < 0)) && all(is.infinite(ub) & (ub > 0))) {
    kerbound <- "none"
  }

  full.lb[icon.idx] <- lb
  full.ub[icon.idx] <- ub

  list(bound = kerbound, lb = full.lb, ub = full.ub)
}

npKernelBoundsCheckEval <- function(evaldat,
                                    icon,
                                    kerlb,
                                    kerub,
                                    argprefix = "cker") {
  if (is.null(evaldat) || is.null(icon) || !any(icon))
    return(invisible(TRUE))

  if (is.null(kerlb) || is.null(kerub))
    return(invisible(TRUE))

  icon.idx <- which(icon)
  ncon <- length(icon.idx)
  if (ncon == 0L)
    return(invisible(TRUE))

  lb <- as.numeric(kerlb[icon.idx])
  ub <- as.numeric(kerub[icon.idx])

  if (all(is.infinite(lb) & lb < 0) && all(is.infinite(ub) & ub > 0))
    return(invisible(TRUE))

  enames <- names(evaldat)
  if (is.null(enames) || length(enames) != length(icon))
    enames <- paste0("V", seq_along(icon))
  cnames <- enames[icon.idx]

  econ <- as.matrix(evaldat[, icon.idx, drop = FALSE])
  if (!is.double(econ))
    storage.mode(econ) <- "double"

  bad.low <- vapply(seq_len(ncon), function(i) {
    is.finite(lb[i]) && any(econ[, i] < lb[i], na.rm = TRUE)
  }, logical(1L))
  bad.high <- vapply(seq_len(ncon), function(i) {
    is.finite(ub[i]) && any(econ[, i] > ub[i], na.rm = TRUE)
  }, logical(1L))

  if (any(bad.low) || any(bad.high)) {
    low.msg <- if (any(bad.low)) {
      paste(sprintf("%s < %s", cnames[bad.low], format(lb[bad.low], digits = 12)), collapse = ", ")
    } else {
      NULL
    }
    high.msg <- if (any(bad.high)) {
      paste(sprintf("%s > %s", cnames[bad.high], format(ub[bad.high], digits = 12)), collapse = ", ")
    } else {
      NULL
    }
    msg <- paste(c(low.msg, high.msg), collapse = "; ")
    stop(sprintf("Evaluation data violate '%s' bounds: %s",
                 paste0(argprefix, "bound"), msg))
  }

  invisible(TRUE)
}

npKernelBoundsMarshal <- function(kerlb, kerub) {
  lb <- as.double(kerlb)
  ub <- as.double(kerub)
  big <- .Machine$double.xmax
  lb[is.infinite(lb) & lb < 0] <- -big
  ub[is.infinite(ub) & ub > 0] <- big
  list(lb = lb, ub = ub)
}

validateBandwidth <- function(bws){
  vari <- names(bws$bandwidth)
  bchecker <- function(j){
    v <- vari[j]
    dati <- bws$dati[[v]]
    bwv <- bws$bandwidth[[j]]
    stopifnot(length(bwv) == length(dati$iord))

    cd <- function(a,b){
      (a-b)/(a+b+.Machine$double.eps) > 5.0*.Machine$double.eps
    }
    
    nn.lower.bound <- if (!identical(v, "x")) 1L else npRegressionNnLowerBound(bws)

    vb <- sapply(seq_along(bwv), function(i){
      falg <- (bwv[i] < 0)

      if (dati$icon[i] && (falg || (!is.finite(bwv[i])))){
        stop(paste("Invalid bandwidth supplied for continuous",
                   "variable", bws$varnames[[v]][i], ":",bwv[i]))
      }

      if (dati$icon[i] &&
          !identical(bws$type, "fixed") &&
          (bwv[i] > 0) &&
          (bwv[i] < nn.lower.bound)) {
        stop(sprintf(
          "%s: nearest-neighbor bandwidth must be at least %d for continuous nonfixed regression",
          bws$varnames[[v]][i],
          nn.lower.bound
        ))
      }
      
      if (dati$iord[i] &&
          (falg || cd(bwv[i],oMaxL(dati$all.nlev[[i]],
                         kertype = bws$klist[[v]]$okertype)))){
        stop(paste("Invalid bandwidth supplied for ordered",
                   "variable", bws$varnames[[v]][i], ":",bwv[i]))
      }
      
      if (dati$iuno[i] &&
          (falg || cd(bwv[i],uMaxL(dati$all.nlev[[i]],
                         kertype = bws$klist[[v]]$ukertype)))){
        stop(paste("Invalid bandwidth supplied for unordered",
                   "variable", bws$varnames[[v]][i], ":",bwv[i]))
      }
      return(TRUE)
    })
    
    return(vb)
  }
  vbl <- lapply(seq_along(vari), bchecker)
  invisible(vbl)
}

validateBandwidthTF <- function(bws){
  vari <- names(bws$bandwidth)
  bchecker <- function(j){
    v <- vari[j]
    dati <- bws$dati[[v]]
    bwv <- bws$bandwidth[[j]]

    if(length(bwv) != length(dati$iord))
      return(FALSE)

    cd <- function(a,b){
      (a-b)/(a+b+.Machine$double.eps) > 5.0*.Machine$double.eps
    }
    
    nn.lower.bound <- if (!identical(v, "x")) 1L else npRegressionNnLowerBound(bws)

    vb <- sapply(seq_along(bwv), function(i){
      falg <- (bwv[i] < 0)

      if (dati$icon[i]) {
        if(bws$type == "fixed") {
          if(falg || (!is.finite(bwv[i]))){
            return(FALSE)
          }
        } else if((((bwv[i] > 0) && (bwv[i] < nn.lower.bound))) || (!is.finite(bwv[i]))) {
          return(FALSE)
        }
      }
      
      if (dati$iord[i] &&
          (falg || cd(bwv[i],oMaxL(dati$all.nlev[[i]],
                         kertype = bws$klist[[v]]$okertype)))){
        return(FALSE)
      }
      
      if (dati$iuno[i] &&
          (falg || cd(bwv[i],uMaxL(dati$all.nlev[[i]],
                         kertype = bws$klist[[v]]$ukertype)))){
        return(FALSE)
      }
      return(TRUE)
    })
    
    return(all(vb))
  }
  vbl <- all(unlist(lapply(seq_along(vari), bchecker)))
  return(vbl)
}

.np_refresh_xy_bandwidth_metadata <- function(bws) {
  if (is.null(bws$xbw) || is.null(bws$ybw) ||
      is.null(bws$xdati) || is.null(bws$ydati) ||
      is.null(bws$nconfac) || is.null(bws$ncatfac) ||
      is.null(bws$sdev) || is.null(bws$dati) || is.null(bws$klist)) {
    return(bws)
  }

  bandwidth <- list(x = bws$xbw, y = bws$ybw)
  sfactor <- bandwidth

  if ((bws$xnuno + bws$ynuno) > 0L) {
    if (isTRUE(bws$scaling)) {
      if (bws$xnuno > 0L)
        bandwidth$x[bws$xdati$iuno] <- bandwidth$x[bws$xdati$iuno] * bws$ncatfac
      if (bws$ynuno > 0L)
        bandwidth$y[bws$ydati$iuno] <- bandwidth$y[bws$ydati$iuno] * bws$ncatfac
    } else {
      if (bws$xnuno > 0L)
        sfactor$x[bws$xdati$iuno] <- sfactor$x[bws$xdati$iuno] / bws$ncatfac
      if (bws$ynuno > 0L)
        sfactor$y[bws$ydati$iuno] <- sfactor$y[bws$ydati$iuno] / bws$ncatfac
    }
  }

  if ((bws$xnord + bws$ynord) > 0L) {
    if (isTRUE(bws$scaling)) {
      if (bws$xnord > 0L)
        bandwidth$x[bws$xdati$iord] <- bandwidth$x[bws$xdati$iord] * bws$ncatfac
      if (bws$ynord > 0L)
        bandwidth$y[bws$ydati$iord] <- bandwidth$y[bws$ydati$iord] * bws$ncatfac
    } else {
      if (bws$xnord > 0L)
        sfactor$x[bws$xdati$iord] <- sfactor$x[bws$xdati$iord] / bws$ncatfac
      if (bws$ynord > 0L)
        sfactor$y[bws$ydati$iord] <- sfactor$y[bws$ydati$iord] / bws$ncatfac
    }
  }

  if ((bws$xncon + bws$yncon) > 0L) {
    sdev <- as.numeric(bws$sdev)
    sx <- if (bws$xncon > 0L) sdev[seq_len(bws$xncon)] else numeric(0L)
    sy <- if (bws$yncon > 0L) sdev[bws$xncon + seq_len(bws$yncon)] else numeric(0L)

    if (isTRUE(bws$scaling)) {
      if (bws$xncon > 0L)
        bandwidth$x[bws$xdati$icon] <- bandwidth$x[bws$xdati$icon] * (sx * bws$nconfac)
      if (bws$yncon > 0L)
        bandwidth$y[bws$ydati$icon] <- bandwidth$y[bws$ydati$icon] * (sy * bws$nconfac)
    } else {
      if (bws$xncon > 0L)
        sfactor$x[bws$xdati$icon] <- sfactor$x[bws$xdati$icon] / (sx * bws$nconfac)
      if (bws$yncon > 0L)
        sfactor$y[bws$ydati$icon] <- sfactor$y[bws$ydati$icon] / (sy * bws$nconfac)
    }
  }

  sumNum <- lapply(c("x", "y"), function(v) {
    dati <- bws$dati[[v]]
    sfv <- sfactor[[v]]
    sapply(seq_along(sfv), function(i) {
      if (dati$icon[i])
        return(sfv[i])
      if (dati$iord[i])
        return(oMaxL(dati$all.nlev[[i]], kertype = bws$klist[[v]]$okertype))
      if (dati$iuno[i])
        return(uMaxL(dati$all.nlev[[i]], kertype = bws$klist[[v]]$ukertype))
      NA_real_
    })
  })
  names(sumNum) <- c("x", "y")

  bws$bandwidth <- bandwidth
  bws$sfactor <- sfactor
  bws$sumNum <- sumNum
  bws
}

npRegressionNnLowerBound <- function(bws) {
  if (!inherits(bws, "rbandwidth"))
    return(1L)

  if (is.null(bws$type) || identical(bws$type, "fixed"))
    return(1L)

  if (is.null(bws$icon) || !any(bws$icon))
    return(1L)

  2L
}

npExtendedNnEnabled <- function() {
  npLogicalOption("np.extendednn", FALSE)
}

npContinuousExtendedNnNomadUpper <- function(traindat,
                                          bwtype,
                                          ckertype,
                                          cont.idx,
                                          evaldat = NULL,
                                          safety.margin = 1.5,
                                          hard.upper = .Machine$integer.max / 4) {
  ncon <- length(cont.idx)
  nobs <- NROW(traindat)
  base.k <- as.integer(nobs) - 1L
  fallback <- rep.int(as.double(base.k), ncon)
  bwtype <- as.character(bwtype)[1L]

  if (!npExtendedNnEnabled() ||
      ncon < 1L ||
      !(bwtype %in% c("generalized_nn", "adaptive_nn")) ||
      base.k < 1L) {
    return(fallback)
  }

  rel.tol <- npLargehRelTol()

  kern <- as.character(ckertype)[1L]
  utol <- switch(kern,
    gaussian = sqrt(-2.0 * log(1.0 - rel.tol)),
    epanechnikov = sqrt(rel.tol),
    uniform = 1.0 - 32.0 * .Machine$double.eps,
    "truncated gaussian" = sqrt(-2.0 * log(1.0 - rel.tol)),
    0.0
  )

  if (!is.finite(utol) || utol <= 0)
    return(fallback)

  traindat <- toFrame(traindat)
  evaldat <- if (is.null(evaldat)) traindat else toFrame(evaldat)
  upper <- fallback

  for (j in seq_along(cont.idx)) {
    train.vals <- as.double(traindat[[cont.idx[j]]])
    train.vals <- train.vals[is.finite(train.vals)]
    eval.vals <- as.double(evaldat[[cont.idx[j]]])
    eval.vals <- eval.vals[is.finite(eval.vals)]
    if (length(train.vals) < 2L || length(eval.vals) < 1L)
      next

    xmin <- min(train.vals)
    xmax <- max(train.vals)
    train.range <- xmax - xmin
    if (!is.finite(train.range) || train.range <= 0)
      next

    eval.saturated <- pmax(abs(eval.vals - xmin), abs(xmax - eval.vals))
    saturated <- if (identical(bwtype, "adaptive_nn")) {
      pmax(abs(train.vals - xmin), abs(xmax - train.vals))
    } else {
      eval.saturated
    }
    base.h.min <- min(saturated[is.finite(saturated) & saturated > 0])
    if (!is.finite(base.h.min) || base.h.min <= 0)
      next

    h.large <- max(eval.saturated[is.finite(eval.saturated)], train.range) / utol
    if (!is.finite(h.large) || h.large <= base.h.min)
      next

    k.upper <- ceiling(as.double(base.k) * h.large / base.h.min * safety.margin)
    if (is.finite(k.upper))
      upper[j] <- min(hard.upper, max(as.double(base.k), k.upper))
  }

  upper
}

npRegressionExtendedNnNomadUpper <- function(xdat,
                                          template,
                                          cont.idx,
                                          safety.margin = 1.5,
                                          hard.upper = .Machine$integer.max / 4) {
  if (!inherits(template, "rbandwidth"))
    return(rep.int(as.double(as.integer(NROW(xdat)) - 1L), length(cont.idx)))

  npContinuousExtendedNnNomadUpper(
    traindat = xdat,
    bwtype = template$type,
    ckertype = template$ckertype,
    cont.idx = cont.idx,
    safety.margin = safety.margin,
    hard.upper = hard.upper
  )
}

npRegressionHasExtendedNn <- function(bws) {
  if (is.null(bws))
    return(FALSE)

  if (!is.null(bws$reg.bws) && npRegressionHasExtendedNn(bws$reg.bws))
    return(TRUE)

  if (!is.null(bws$tau.bws) && is.list(bws$tau.bws)) {
    if (any(vapply(bws$tau.bws, npRegressionHasExtendedNn, logical(1L))))
      return(TRUE)
  }

  if (!is.null(bws$bw) && is.list(bws$bw) && !is.data.frame(bws$bw)) {
    child.is.bw <- vapply(
      bws$bw,
      function(x) any(c("rbandwidth", "plbandwidth", "scbandwidth", "sibandwidth") %in% class(x)),
      logical(1L)
    )
    if (any(child.is.bw) &&
        any(vapply(bws$bw[child.is.bw], npRegressionHasExtendedNn, logical(1L))))
      return(TRUE)
    if (any(child.is.bw))
      return(FALSE)
  }

  if (is.null(bws$type) ||
      !(as.character(bws$type)[1L] %in% c("generalized_nn", "adaptive_nn")) ||
      is.null(bws$nobs)) {
    return(FALSE)
  }

  upper <- as.double(as.integer(bws$nobs) - 1L)
  if (!is.finite(upper) || upper < 1)
    return(FALSE)

  if (!is.null(bws$bw) && !is.null(bws$icon) && any(bws$icon)) {
    bw <- suppressWarnings(as.double(unlist(bws$bw, use.names = FALSE)))
    icon <- as.logical(bws$icon)
    if (length(bw) >= length(icon) && any(is.finite(bw[icon]) & bw[icon] > upper))
      return(TRUE)
  }

  if (!is.null(bws$bw) && is.null(bws$icon) &&
      !is.null(bws$ncon) && isTRUE(as.integer(bws$ncon)[1L] > 0L)) {
    bw <- suppressWarnings(as.double(unlist(bws$bw, use.names = FALSE)))
    ncon <- as.integer(bws$ncon)[1L]
    if (length(bw) >= ncon && any(is.finite(bw[seq_len(ncon)]) & bw[seq_len(ncon)] > upper))
      return(TRUE)
  }

  if (!is.null(bws$xbw) && !is.null(bws$ixcon) && any(bws$ixcon)) {
    xbw <- suppressWarnings(as.double(bws$xbw))
    ixcon <- as.logical(bws$ixcon)
    if (length(xbw) >= length(ixcon) && any(is.finite(xbw[ixcon]) & xbw[ixcon] > upper))
      return(TRUE)
  }

  if (!is.null(bws$ybw) && !is.null(bws$iycon) && any(bws$iycon)) {
    ybw <- suppressWarnings(as.double(bws$ybw))
    iycon <- as.logical(bws$iycon)
    if (length(ybw) >= length(iycon) && any(is.finite(ybw[iycon]) & ybw[iycon] > upper))
      return(TRUE)
  }

  FALSE
}

npValidateExtendedNnContinuousBandwidth <- function(bws,
                                                 where,
                                                 nobs = NULL) {
  if (is.null(bws$type) ||
      identical(as.character(bws$type)[1L], "fixed") ||
      is.null(bws$icon) ||
      !any(bws$icon)) {
    return(invisible(bws))
  }

  nobs <- if (is.null(nobs)) bws$nobs else nobs
  if (is.null(nobs))
    return(invisible(bws))

  upper <- as.integer(nobs) - 1L
  if (!is.finite(upper) || upper < 1L)
    return(invisible(bws))

  bw <- as.double(bws$bw)
  icon <- which(as.logical(bws$icon))
  offenders <- is.finite(bw[icon]) & bw[icon] > upper
  if (!any(offenders))
    return(invisible(bws))

  bwtype <- as.character(bws$type)[1L]
  if (!(bwtype %in% c("generalized_nn", "adaptive_nn"))) {
    stop(
      sprintf(
        "%s: extended nearest-neighbor bandwidths above n-1 are not enabled for bwtype='%s'",
        where,
        bwtype
      ),
      call. = FALSE
    )
  }

  if (!npExtendedNnEnabled()) {
    stop(
      sprintf(
        "%s: nearest-neighbor bandwidth exceeds n-1; set options(np.extendednn = TRUE) to allow extended generalized_nn/adaptive_nn bandwidths",
        where
      ),
      call. = FALSE
    )
  }

  invisible(bws)
}

npValidateConditionalExtendedNn <- function(bws,
                                         where,
                                         nobs = NULL) {
  if (is.null(bws$type) ||
      identical(as.character(bws$type)[1L], "fixed") ||
      is.null(bws$ixcon) ||
      is.null(bws$iycon)) {
    return(invisible(bws))
  }

  nobs <- if (is.null(nobs)) bws$nobs else nobs
  if (is.null(nobs))
    return(invisible(bws))

  upper <- as.integer(nobs) - 1L
  if (!is.finite(upper) || upper < 1L)
    return(invisible(bws))

  xbw <- as.double(bws$xbw)
  ybw <- as.double(bws$ybw)
  xicon <- which(as.logical(bws$ixcon))
  yicon <- which(as.logical(bws$iycon))
  offenders <- FALSE
  if (length(xicon))
    offenders <- offenders || any(is.finite(xbw[xicon]) & xbw[xicon] > upper)
  if (length(yicon))
    offenders <- offenders || any(is.finite(ybw[yicon]) & ybw[yicon] > upper)

  if (!isTRUE(offenders))
    return(invisible(bws))

  bwtype <- as.character(bws$type)[1L]
  if (!(bwtype %in% c("generalized_nn", "adaptive_nn"))) {
    stop(
      sprintf(
        "%s: extended nearest-neighbor bandwidths above n-1 are not enabled for bwtype='%s'",
        where,
        bwtype
      ),
      call. = FALSE
    )
  }

  if (!npExtendedNnEnabled()) {
    stop(
      sprintf(
        "%s: nearest-neighbor bandwidth exceeds n-1; set options(np.extendednn = TRUE) to allow extended generalized_nn/adaptive_nn bandwidths",
        where
      ),
      call. = FALSE
    )
  }

  invisible(bws)
}

npValidateRegressionExtendedNn <- function(bws,
                                        where,
                                        bandwidth.compute = FALSE) {
  if (!inherits(bws, "rbandwidth") ||
      is.null(bws$type) ||
      identical(as.character(bws$type)[1L], "fixed") ||
      is.null(bws$icon) ||
      !any(bws$icon) ||
      is.null(bws$nobs)) {
    return(invisible(bws))
  }

  upper <- as.integer(bws$nobs) - 1L
  if (!is.finite(upper) || upper < 1L)
    return(invisible(bws))

  bw <- as.double(bws$bw)
  icon <- which(as.logical(bws$icon))
  offenders <- is.finite(bw[icon]) & bw[icon] > upper
  if (!any(offenders))
    return(invisible(bws))

  bwtype <- as.character(bws$type)[1L]
  if (!(bwtype %in% c("generalized_nn", "adaptive_nn"))) {
    stop(
      sprintf(
        "%s: extended nearest-neighbor bandwidths above n-1 are not enabled for bwtype='%s'",
        where,
        bwtype
      ),
      call. = FALSE
    )
  }

  if (!npExtendedNnEnabled()) {
    stop(
      sprintf(
        "%s: nearest-neighbor bandwidth exceeds n-1; set options(np.extendednn = TRUE) to allow extended generalized_nn/adaptive_nn bandwidths",
        where
      ),
      call. = FALSE
    )
  }

  invisible(bws)
}

npValidateRegressionNnLowerBound <- function(bws,
                                             where,
                                             allow.zero.placeholder = FALSE) {
  lower <- npRegressionNnLowerBound(bws)
  if (lower <= 1L)
    return(invisible(bws))

  bw <- as.double(bws$bw)
  icon <- which(as.logical(bws$icon))
  if (!length(icon))
    return(invisible(bws))

  offenders <- bw[icon] < lower
  if (allow.zero.placeholder)
    offenders <- offenders & (bw[icon] > 0)
  if (!any(offenders))
    return(invisible(bws))

  upper <- max(lower, as.integer(bws$nobs) - 1L)
  stop(
    sprintf(
      "%s: nearest-neighbor bandwidth must be in [%d, %d] for continuous nonfixed regression",
      where,
      lower,
      upper
    ),
    call. = FALSE
  )
}


.np_formula_quote_name_if_needed <- function(x) {
  reserved <- c("if", "else", "repeat", "while", "function", "for", "in",
                "next", "break", "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA",
                "NA_integer_", "NA_real_", "NA_complex_", "NA_character_")
  vapply(x, function(term) {
    if (make.names(term) == term && !(term %in% reserved))
      return(term)

    paste0("`", gsub("`", "\\\\`", term), "`")
  }, character(1), USE.NAMES = FALSE)
}

explodeFormula <- function(formula, data=NULL){
  formula.terms <- if (is.null(data)) terms(formula) else terms(formula, data = data)
  response <- if (length(formula) == 3L) all.vars(formula[[2L]]) else character(0)
  term.labels.raw <- attr(formula.terms, "term.labels")
  term.labels <- ifelse(grepl("^`[^`]+`$", term.labels.raw),
                        substring(term.labels.raw, 2L, nchar(term.labels.raw) - 1L),
                        term.labels.raw)
  res <- list(response, term.labels)
  stopifnot(all(sapply(res,length) > 0))
  names(res) <- c("response","terms")
  attr(res, "formula.labels") <- list(
    response = if (length(formula) == 3L) .np_formula_quote_name_if_needed(response) else character(0),
    terms = term.labels.raw
  )
  res
}

.np_formula_quote_if_needed <- function(x) {
  vapply(x, function(term) {
    parsed <- try(parse(text = term), silent = TRUE)
    if (!inherits(parsed, "try-error"))
      return(term)

    paste0("`", gsub("`", "\\\\`", term), "`")
  }, character(1), USE.NAMES = FALSE)
}


explodePipe <- function(formula, env = parent.frame()){
  if (!inherits(formula, "formula")) {
    if (is.symbol(formula) && is.environment(env)) {
      not_found <- .np_missing_binding_sentinel
      sym_val <- get0(as.character(formula), envir = env, inherits = TRUE, ifnotfound = not_found)
      if (!identical(sym_val, not_found)) {
        formula <- sym_val
      } else {
        out <- .np_try_eval_in_frames(formula, eval_env = env, search_frames = FALSE)
        if (isTRUE(out$ok))
          formula <- out$value
        else
          formula <- out$error
      }
    } else {
      out <- .np_try_eval_in_frames(formula, eval_env = env, search_frames = FALSE)
      if (isTRUE(out$ok))
        formula <- out$value
      else
        formula <- out$error
    }
    if (inherits(formula, "error"))
      stop(conditionMessage(formula), call. = FALSE)
  }
  tf <- as.character(formula)  
  tf <- tf[length(tf)]
  lhs <- if (length(as.character(formula)) == 3) {
    strsplit(as.character(formula)[2], " *[+] *")
  } else {
    list()
  }
  rhs <- strsplit(strsplit(tf, " *[|] *")[[1]], " *[+] *")
  c(lhs, rhs)
}

"%~%" <- function(a,b) {
  identical(class(a), class(b)) && (length(a) == length(b)) &&
  all(vapply(a, coarseclass, character(1)) == vapply(b, coarseclass, character(1)))
}

coarseclass <- function(a) {
  if (inherits(a, "integer")) return("numeric")
  return(class(a)[1])
}

toFrame <- function(frame) {
  if(!is.data.frame(frame)){
    t.names <- NULL

    if(!(is.vector(frame) || is.factor(frame) || is.matrix(frame)))
      stop(deparse(substitute(frame))," must be a data frame, matrix, vector, or factor")

    if(!is.matrix(frame))
      t.names <- paste(deparse(substitute(frame)),
                       collapse = "")
    
    frame <- data.frame(frame, check.names=FALSE)
    
    if(!is.null(t.names))
      names(frame) <- t.names
  }
  return(frame)
}


cast <- function(a, b, same.levels = TRUE){
  if(is.ordered(b)){
    if(same.levels)
      ordered(a, levels = levels(b))
    else
      ordered(a)
  }   
  else if(is.factor(b)){
    if(same.levels)
      factor(a, levels = levels(b))
    else
      factor(a)
  }
  else if (coarseclass(b) == "numeric")
    as.double(a)
  else if (is.data.frame(b)) {
    if (dim(a)[2] == dim(b)[2]){
      r = data.frame(a)
      for (i in seq_along(b))
        r[,i] = cast(a[,i],b[,i], same.levels = same.levels)
      r
    } else { stop("a could not be cast as b") }
  }
}

subcol <- function(x, v, i){
  x[,i] = cast(v,x[,i])
  x
}

mcvConstruct <- function(dati){
  nuno <- sum(dati$iuno)
  nord <- sum(dati$iord)

  num.row <- max(sapply(dati$all.lev,length))

  pad.num <- numNotIn(unlist(dati$all.dlev))

  mcv <- matrix(data = pad.num, nrow = num.row, ncol = (nuno+nord))
  attr(mcv, "num.row") <- num.row
  attr(mcv, "pad.num") <- pad.num

  cnt <- 0
  if (nuno > 0)
    for (i in which(dati$iuno)) 
      mcv[seq_along(dati$all.lev[[i]]), (cnt <- cnt+1)] <- dati$all.dlev[[i]]

  cnt <- 0
  if (nord > 0)
    for (i in which(dati$iord))
      mcv[seq_along(dati$all.lev[[i]]), (cnt <- cnt+1)+nuno] <- dati$all.dlev[[i]]

  mcv
}

## when admitting new categories, adjustLevels will attempt to catch possible mistakes:
## if an unordered variable contains more than one new category, warn
## if an ordered, but scaleless variable contains a new category, error
## if an ordered, scale-possessing variable contains a new category lacking scale, error

adjustLevels <- function(data, dati, allowNewCells = FALSE){
  for (i in which(dati$iord | dati$iuno)){
    if (allowNewCells){
      newCats <- setdiff(levels(data[,i]),dati$all.lev[[i]])
      if (length(newCats) >= 1){
        if (dati$iuno[i]){
          if (length(newCats) > 1)
            .np_warning(paste("more than one 'new' category is redundant when estimating on unordered data.\n",
                          "training data categories: ", paste(dati$all.lev[[i]], collapse=" "),"\n",
                          "redundant estimation data categories: ", paste(newCats, collapse=" "), "\n", sep=""))
          data[,i] <- factor(data[,i], levels = c(dati$all.lev[[i]], newCats))
        } else {
          if (dati$inumord[i]){
            if (!isNum(newCats))
              stop(paste("estimation data contains a new qualitative category, but training data is\n",
                         "ordered, and numeric.\n",
                         "training data categories: ", paste(dati$all.lev[[i]], collapse=" "),"\n",
                         "conflicting estimation data categories: ", paste(newCats, collapse=" "), "\n", sep=""))
          } else {
            stop(paste("estimation beyond the support of training data of an ordered,\n",
                       "categorical, qualitative variable is not supported.\n"))
          }

          data[,i] <- ordered(data[,i], levels = sort(as.numeric(c(dati$all.lev[[i]], newCats))))
        }
      } else {
        data[,i] <- factor(data[,i], levels = dati$all.lev[[i]])
      }
    } else {
      if (!all(is.element(levels(data[,i]), dati$all.lev[[i]])))
        stop("data contains unknown factors (wrong dataset provided?)")
      data[,i] <- factor(data[,i], levels = dati$all.lev[[i]])
    }
  }

  data
}

toMatrix <- function(data) {
  tq <- sapply(data, function(y) {
    if(is.factor(y))
      y <- dlev(y)[as.integer(y)]
    y})
  dim(tq) <- dim(data) ## cover the corner case of single element d.f.
  tq
}

## this doesn't just strictly check for the response, but does tell you
## that evaluating with response fails ... in principle the evaluation
## could fail without the response too, but then the calling routine is about
## to die a noisy death anyhow ...
succeedWithResponse <- function(tt, frame){
  vars <- attr(tt, "variables")
  if (is.null(vars))
    return(FALSE)
  out <- .np_try_eval_in_frames(vars, eval_env = frame, enclos = NULL, search_frames = FALSE)
  isTRUE(out$ok)
}

## determine whether a bandwidth
## matches a data set
bwMatch <- function(data, dati){
  if (length(dati$icon) != ncol(data))
    stop("bandwidth vector is of improper length")

  test.dati <- untangle(data)

  if (any(xor(dati$icon,test.dati$icon)) ||
      any(xor(dati$iord,test.dati$iord)) ||
      any(xor(dati$iuno,test.dati$iuno)))
    stop(paste("supplied bandwidths do not match","data", "in type"))
}

uMaxL <- function(c, kertype = c("aitchisonaitken","liracine")){
  switch(kertype,
         aitchisonaitken = (c-1.0)/c,
         liracine = 1.0)
}

oMaxL <- function(c, kertype = c("wangvanryzin", "liracine", "nliracine", "racineliyan")){
  switch(kertype,
         wangvanryzin = 1.0,
         liracine = 1.0,
         nliracine = 1.0,
         "racineliyan" = 1.0)
}

## tested with: rbandwidth
## right now all bandwidth objects have some crusty
## vestiges of their evolution, ie. non-list metadata
## such as xnames or ynames. The new metadata system is
## for the most part list based and facilitates generic
## operations

updateBwNameMetadata <- function(nameList, bws){
  ## names of 'old' metadata in bw object
  onames <- names(nameList)
  if (length(onames) == 0L)
    return(bws)

  for (nm in onames) {
    bws[[nm]] <- nameList[[nm]]
    bws$varnames[[substr(nm, 1L, 1L)]] <- nameList[[nm]]
  }

  bws
}

## some string utility functions

pad <- function(s){
  idx <- nchar(s) > 0
  out <- rep.int(" ", length(s))
  out[idx] <- paste("", s[idx], "")
  names(out) <- names(s)
  out
}

rpad <- function(s){
  idx <- nchar(s) > 0
  out <- rep.int("", length(s))
  out[idx] <- paste(s[idx], "")
  names(out) <- names(s)
  out
}

blank <- function(len){
  sapply(len, function(nb){
    paste(rep(' ', times = nb), collapse='')
  })
}

formatv <- function(v){
  sapply(seq_along(v), function(j){ format(v[j]) })
}

## strings used in various report generating functions

genOmitStr <- function(x){
  t.str <- ''
  if(!is.null(x$rows.omit) && !identical(x$rows.omit, NA))
    t.str <- paste("\nNo. Complete Observations: ", x$nobs,
                   "\nNo. Incomplete (NA) Observations: ", x$nobs.omit,
                   "\nObservations omitted or excluded: ", paste(x$rows.omit, collapse=" "),
                   "\n")
  return(t.str)
}

## Estimation-related rgf's
genGofStr <- function(x){
  mse.str <- ""
  if (!is.na(x$MSE))
    mse.str <- paste("\nResidual standard error:", format(sqrt(x$MSE)))

  r2.str <- ""
  if (!is.na(x$R2))
    r2.str <- paste("\nR-squared:", format(x$R2))

  paste(mse.str, r2.str, sep = "")
}

.np_timing_optim_label <- function(x) {
  degree.search <- if (is.list(x$degree.search)) {
    x$degree.search
  } else if (is.list(x$bws) && is.list(x$bws$degree.search)) {
    x$bws$degree.search
  } else {
    NULL
  }

  if (is.list(degree.search)) {
    labels <- character(0)
    for (field in c("mode", "method", "engine")) {
      value <- degree.search[[field]]
      if (!is.null(value) && length(value))
        labels <- c(labels, as.character(value[1L]))
    }
    labels <- labels[!is.na(labels)]
    if ("exhaustive" %in% labels)
      return("Exhaustive Powell")
  }

  "optim"
}

genTimingStr <- function(x){
  if (is.null(x$total.time) || is.na(x$total.time))
    return("")

  nomad.time <- if (!is.null(x$nomad.time) && !is.na(x$nomad.time))
    as.double(x$nomad.time) else NA_real_
  powell.time <- if (!is.null(x$powell.time) && !is.na(x$powell.time))
    as.double(x$powell.time) else NA_real_
  fit.time <- if (!is.null(x$fit.time) && !is.na(x$fit.time))
    as.double(x$fit.time) else NA_real_

  if (is.finite(nomad.time) || is.finite(powell.time)) {
    detail <- character(0)
    if (is.finite(nomad.time))
      detail <- c(detail, paste("NOMAD ", format(nomad.time), "s", sep = ""))
    if (is.finite(powell.time)) {
      optim.method <- .np_summary_r_optim_method(x)
      refinement.label <- if (is.na(optim.method)) {
        "Powell"
      } else {
        paste0("R optim: ", optim.method)
      }
      detail <- c(detail, paste0(refinement.label, " ",
                                 format(powell.time), "s"))
    }
    if (is.finite(fit.time))
      detail <- c(detail, paste("fit ", format(fit.time), "s", sep = ""))

    return(paste("\nEstimation Time: ", format(x$total.time), " seconds (",
                 paste(detail, collapse = ", "), ")", sep = ""))
  }

  optim.label <- .np_timing_optim_label(x)
  optim.time <- if (!is.null(x$optim.time) && !is.na(x$optim.time)) {
    as.double(x$optim.time)
  } else if (!identical(optim.label, "optim")) {
    as.double(x$total.time)
  } else {
    NA_real_
  }

  if (is.finite(optim.time) && !is.null(x$fit.time) && !is.na(x$fit.time))
    return(paste("\nEstimation Time: ", format(x$total.time), " seconds (",
                 optim.label, " ", format(optim.time), "s, fit ",
                 format(x$fit.time), "s)", sep = ""))

  if (is.finite(optim.time) && !identical(optim.label, "optim"))
    return(paste("\nEstimation Time: ", format(x$total.time), " seconds (",
                 optim.label, " ", format(optim.time), "s)", sep = ""))

  paste("\nEstimation Time: ",format(x$total.time)," seconds",sep = "")
}

.np_attach_nomad_restart_summary <- function(bws, search_result, tol = 1e-10) {
  restart.results <- search_result$restart.results
  direction <- if (!is.null(search_result$direction) && length(search_result$direction)) {
    as.character(search_result$direction[1L])
  } else {
    "min"
  }
  restart.fval <- if (!is.null(restart.results) && length(restart.results)) {
    vapply(
      restart.results,
      function(x) {
        if (is.null(x$objective) || !length(x$objective)) NA_real_ else as.numeric(x$objective[1L])
      },
      numeric(1L)
    )
  } else {
    numeric(0L)
  }

  best.objective <- if (!is.null(search_result$best$objective) && length(search_result$best$objective)) {
    as.numeric(search_result$best$objective[1L])
  } else {
    NA_real_
  }

  best.restart <- NA_integer_
  if (!is.null(search_result$best.restart) &&
      length(search_result$best.restart) &&
      is.finite(search_result$best.restart[1L])) {
    best.restart <- as.integer(search_result$best.restart[1L])
  } else if (length(restart.fval) && is.finite(best.objective)) {
    diffs <- abs(restart.fval - best.objective)
    diffs[!is.finite(diffs)] <- Inf
    if (any(is.finite(diffs))) {
      idx <- which.min(diffs)
      scale <- max(1, abs(best.objective))
      if (is.finite(diffs[idx]) && diffs[idx] <= tol * scale)
        best.restart <- as.integer(idx)
    }
  } else if (length(restart.fval) && any(is.finite(restart.fval))) {
    best.restart <- as.integer(if (identical(direction, "max")) which.max(restart.fval) else which.min(restart.fval))
  }

  bws$nomad.restart.results <- restart.results
  bws$nomad.restart.fval <- restart.fval
  bws$nomad.best.restart <- best.restart
  bws$nomad.restart.starts <- search_result$restart.starts
  bws$nomad.restart.degree.starts <- search_result$restart.degree.starts
  bws$nomad.restart.bandwidth.starts <- search_result$restart.bandwidth.starts
  bws$nomad.restart.start.info <- search_result$restart.start.info
  bws$nomad.remin <- isTRUE(search_result$nomad.remin)
  bws$nomad.remin.index <- search_result$nomad.remin.index
  bws$nomad.remin.roundtrip <- search_result$nomad.remin.roundtrip
  bws
}

.np_multistart_label <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (!is.null(x$nomad.best.restart) &&
      length(x$nomad.best.restart) == 1L &&
      is.finite(x$nomad.best.restart)) {
    idx <- as.numeric(x$nomad.best.restart[1L])
    if (abs(idx - round(idx)) <= tol && idx >= 1)
      return(as.integer(round(idx)))
  }

  if (!is.null(x$ifval) &&
      length(x$ifval) == 1L &&
      is.finite(x$ifval)) {
    idx <- as.numeric(x$ifval[1L])
    if (abs(idx - round(idx)) <= tol) {
      idx <- as.integer(round(idx))
      return(if (idx <= 0L) 1L else idx)
    }
  }

  NA_integer_
}
  
pCatGofStr <- function(x){
  if(!identical(x$confusion.matrix, NA)){
    cat("\nConfusion Matrix\n")
    print(x$confusion.matrix)
  }

  if (!identical(x$CCR.overall,NA))
    cat("\nOverall Correct Classification Ratio: ", format(x$CCR.overall))

  if (!identical(x$CCR.byoutcome,NA)){
    cat("\nCorrect Classification Ratio By Outcome:\n")
    print(x$CCR.byoutcome)
  }

  if (!identical(x$fit.mcfadden,NA))
    cat("\nMcFadden-Puig-Kerschner performance measure: ", format(x$fit.mcfadden))

}

genDenEstStr <- function(x){
  loglik.str <- ""
  if (!(is.null(x$log_likelihood) || identical(x$log_likelihood, NA))) {
    loglik.str <- paste("\nLog Likelihood:", format(x$log_likelihood))
  }
  paste("\nBandwidth Type: ", x$ptype, loglik.str, sep = "")
}

genRegEstStr <- function(x){
  x_bws <- if (is.list(x)) x[["bws", exact = TRUE]] else NULL
  regtype <- if (!is.null(x$regtype)) x$regtype else if (!is.null(x_bws)) x_bws$regtype else NULL
  basis <- if (!is.null(x$basis)) x$basis else if (!is.null(x_bws)) x_bws$basis else NULL
  bern <- if (!is.null(x$bernstein.basis)) x$bernstein.basis else if (!is.null(x_bws)) x_bws$bernstein.basis else NULL
  est.label <- if (identical(regtype, "lp")) npFormatRegressionType(x) else x$pregtype
  basis.family <- if (identical(regtype, "lp")) npLpBasisFamilyLabel(basis) else NULL
  basis.rep <- if (identical(regtype, "lp")) npLpBasisRepresentationLabel(bern) else NULL
  est.prefix <- if (!is.null(x_bws) && npUsesPolynomialSummaryLabel(x_bws))
    "Polynomial Type"
  else
    "Kernel Regression Estimator"
  est.label.str <- if (is.null(est.label)) "" else paste("\n", est.prefix, ": ", est.label, sep = "")
  basis.family.str <- if (is.null(basis.family)) "" else paste("\nLP Basis Family:", basis.family)
  basis.rep.str <- if (is.null(basis.rep)) "" else paste("\nLP Basis Representation:", basis.rep)
  ptype.str <- if (is.null(x$ptype)) "" else paste("\nBandwidth Type:", x$ptype)
  tau.str <- if (is.null(x$tau)) "" else paste("\nTau:", paste(x$tau, collapse = ", "))
  paste(est.label.str, basis.family.str, basis.rep.str, ptype.str, tau.str, sep = "")
}

npLpBasisFamilyLabel <- function(basis){
  basis.value <- if (is.null(basis) || !length(basis)) "glp" else basis
  b <- tolower(basis.value)
  switch(b,
         glp = "Generalized",
         additive = "Additive",
         tensor = "Tensor",
         tools::toTitleCase(b))
}

npLpBasisRepresentationLabel <- function(bernstein){
  if (isTRUE(bernstein)) "Bernstein" else "Raw"
}

npLpBasisNcol <- function(basis = "glp", degree){
  if (is.null(degree) || !length(degree))
    return(NA_real_)
  d <- dim_basis(basis = basis,
                 kernel = TRUE,
                 degree = as.integer(degree),
                 segments = rep.int(1L, length(degree)))
  if (identical(tolower(basis), "tensor")) d else d + 1.0
}

npFormatRegressionType <- function(x){
  xbws <- if (is.list(x)) x[["bws", exact = TRUE]] else NULL
  regtype <- if (!is.null(x$regtype)) {
    x$regtype
  } else if (!is.null(xbws) && !is.null(xbws$regtype)) {
    xbws$regtype
  } else {
    NULL
  }

  pregtype <- if (!is.null(x$pregtype)) {
    x$pregtype
  } else if (!is.null(xbws) && !is.null(xbws$pregtype)) {
    xbws$pregtype
  } else {
    NULL
  }

  if (!identical(regtype, "lp"))
    return(pregtype)

  degree <- if (!is.null(x$degree)) {
    x$degree
  } else if (!is.null(xbws) && !is.null(xbws$degree)) {
    xbws$degree
  } else {
    NULL
  }

  if (is.null(degree) || length(degree) == 0)
    return("Local-Polynomial")

  basis <- if (!is.null(x$basis)) {
    x$basis
  } else if (!is.null(xbws) && !is.null(xbws$basis)) {
    xbws$basis
  } else {
    "glp"
  }
  basis.family <- npLpBasisFamilyLabel(basis)

  if (!is.null(x$child.degree.common) && !isTRUE(x$child.degree.common))
    return(sprintf("Local-Polynomial (%s basis; child-specific degree)",
                   basis.family))

  sprintf("Local-Polynomial (%s basis; degree = %s)",
          basis.family, paste(degree, collapse = ","))
}

npBandwidthSummaryLabel <- function(bwtype, bwscaling = FALSE){
  if (isTRUE(bwscaling))
    return("Scale Factor(s)")

  if (identical(bwtype, "fixed"))
    return("Bandwidth(s)")

  "NN Index(s)"
}

.np_search_param_type <- function(bwtype, is.continuous) {
  if (isTRUE(is.continuous)) {
    if (identical(as.character(bwtype)[1L], "fixed")) "Bandwidth" else "NN Index"
  } else {
    "Lambda"
  }
}

.np_nn_sample_max <- function(bws) {
  if (is.null(bws) || is.null(bws$nobs))
    return(NA_real_)
  out <- suppressWarnings(as.double(as.integer(bws$nobs)[1L] - 1L))
  if (is.finite(out) && out >= 1) out else NA_real_
}

.np_search_param_format <- function(x) {
  trimws(format(sapply(x, format), trim = TRUE), which = "both")
}

.np_text_lpad <- function(x, width) {
  x <- as.character(x)
  paste0(blank(pmax(0L, width - nchar(x))), x)
}

.np_text_rpad <- function(x, width) {
  x <- as.character(x)
  paste0(x, blank(pmax(0L, width - nchar(x))))
}

.np_search_param_table_lines <- function(mat) {
  mat <- as.matrix(mat)
  col.widths <- vapply(seq_len(ncol(mat)), function(j) {
    max(nchar(c(colnames(mat)[j], mat[, j])), na.rm = TRUE)
  }, numeric(1L))
  row.labels <- rownames(mat)
  row.width <- max(nchar(row.labels), na.rm = TRUE)

  out <- paste0(blank(row.width), " ",
                paste(.np_text_lpad(colnames(mat), col.widths),
                      collapse = "  "))
  rows <- vapply(seq_len(nrow(mat)), function(i) {
    paste0(.np_text_rpad(row.labels[i], row.width), " ",
           paste(.np_text_lpad(mat[i, ], col.widths), collapse = "  "))
  }, character(1L))
  c(out, rows)
}

.np_search_param_table_fits <- function(out, width = getOption("width", 80L)) {
  if (!length(out))
    return(TRUE)
  max(nchar(out), na.rm = TRUE) <= width
}

printSearchParameterSummary <- function(values,
                                        varnames,
                                        bws,
                                        vari = "x",
                                        role = NULL,
                                        fallback.label = NULL,
                                        digits = NULL) {
  if (is.null(bws) || is.null(values) || is.null(varnames)) {
    print(matrix(values, ncol = length(values),
                 dimnames = list(fallback.label, varnames)))
    return(invisible(NULL))
  }

  dat <- if (!is.null(bws$dati) && !is.null(bws$dati[[vari]])) {
    bws$dati[[vari]]
  } else if (identical(vari, "x") && !is.null(bws$xdati)) {
    bws$xdati
  } else if (identical(vari, "y") && !is.null(bws$ydati)) {
    bws$ydati
  } else {
    NULL
  }

  if (is.null(dat) || is.null(dat$icon) || length(dat$icon) != length(values)) {
    print(matrix(values, ncol = length(values),
                 dimnames = list(fallback.label, varnames)))
    return(invisible(NULL))
  }

  values <- as.double(values)
  varnames <- as.character(varnames)
  is.cont <- as.logical(dat$icon)
  types <- vapply(is.cont, function(z) .np_search_param_type(bws$type, z), character(1L))

  # Preserve the compact historical display when all entries are plain
  # continuous fixed bandwidths.
  if (all(types == "Bandwidth")) {
    print(matrix(values, ncol = length(values),
                 dimnames = list(fallback.label, varnames)))
    return(invisible(NULL))
  }

  max.values <- rep(NA_real_, length(values))
  sample.max <- .np_nn_sample_max(bws)
  nn.idx <- which(types == "NN Index")
  if (length(nn.idx) && is.finite(sample.max))
    max.values[nn.idx] <- sample.max

  sum.num <- if (!is.null(bws$sumNum) && !is.null(bws$sumNum[[vari]])) {
    suppressWarnings(as.double(bws$sumNum[[vari]]))
  } else {
    rep(NA_real_, length(values))
  }
  if (length(sum.num) == length(values)) {
    lambda.idx <- which(types == "Lambda")
    max.values[lambda.idx] <- sum.num[lambda.idx]
  }

  value.chr <- .np_search_param_format(values)
  max.chr <- ifelse(is.finite(max.values), .np_search_param_format(max.values), "--")
  extended <- rep(FALSE, length(values))
  if (length(nn.idx) && is.finite(sample.max))
    extended[nn.idx] <- is.finite(values[nn.idx]) & values[nn.idx] > sample.max

  mat <- rbind(Type = types, Value = value.chr)
  if (any(is.finite(max.values)))
    mat <- rbind(mat, Max = max.chr)
  if (any(extended))
    mat <- rbind(mat, Extended = ifelse(extended, "yes", "--"))
  colnames(mat) <- varnames

  title <- if (is.null(role) || !nzchar(role)) {
    "Search Parameter(s):"
  } else {
    paste0(role, " Search Parameter(s):")
  }
  out <- .np_search_param_table_lines(mat)
  if (.np_search_param_table_fits(out)) {
    cat(title, "\n", paste(out, collapse = "\n"), "\n", sep = "")
    return(invisible(NULL))
  }

  cat(title)
  for (i in seq_along(values)) {
    pieces <- c(
      paste0("Type: ", types[i]),
      paste0("Value: ", value.chr[i])
    )
    if (is.finite(max.values[i]))
      pieces <- c(pieces, paste0("Max: ", max.chr[i]))
    if (isTRUE(extended[i]))
      pieces <- c(pieces, "Extended: yes")
    cat("\n", varnames[i], " ", paste(pieces, collapse = " "), sep = "")
  }
  cat("\n")
  invisible(NULL)
}

npPolynomialSummaryLabel <- function(x){
  if (npUsesPolynomialSummaryLabel(x))
    "Polynomial Type"
  else
    "Regression Type"
}

npUsesPolynomialSummaryLabel <- function(x){
  density.classes <- c("bandwidth", "dbandwidth", "conbandwidth", "condbandwidth")
  any(class(x) %in% density.classes)
}

.np_summary_missing_number <- function(x) {
  length(x) != 1L || is.na(x) || !is.finite(x)
}

.np_summary_format_count <- function(x) {
  format(round(as.numeric(x)[1L]), big.mark = ",", trim = TRUE,
         scientific = FALSE)
}

.np_summary_format_percent <- function(hits, requests) {
  format(round(100 * as.numeric(hits)[1L] / as.numeric(requests)[1L], 1L),
         nsmall = 1L, trim = TRUE)
}

.np_summary_evaluation_cache_line <- function(stage, hits, lookups) {
  if (.np_summary_missing_number(hits) ||
      .np_summary_missing_number(lookups) ||
      hits <= 0 || lookups <= 0 || lookups < hits)
    return(NULL)

  label <- if (is.null(stage) || !length(stage) || is.na(stage[1L]) ||
               !nzchar(as.character(stage[1L]))) {
    "Evaluation cache"
  } else {
    paste0("Evaluation cache (", as.character(stage[1L]), ")")
  }

  paste0(label, ": ", .np_summary_format_count(hits), " hits / ",
         .np_summary_format_count(lookups), " lookups (",
         .np_summary_format_percent(hits, lookups), "%)")
}

.np_summary_fast_route_line <- function(fast, computed) {
  if (.np_summary_missing_number(fast) ||
      .np_summary_missing_number(computed) ||
      fast <= 0 || computed <= 0 || fast > computed)
    return(NULL)

  paste0("Fast CV route: ", .np_summary_format_count(fast),
         " of ", .np_summary_format_count(computed),
         " function evaluations (",
         .np_summary_format_percent(fast, computed), "%)")
}

.np_summary_named_number <- function(x, name) {
  if (is.null(x) || is.null(names(x)) || !(name %in% names(x)))
    return(NA_real_)
  out <- suppressWarnings(as.numeric(x[[name]][1L]))
  if (length(out) != 1L) NA_real_ else out
}

.np_summary_nomad_cache_counts <- function(x) {
  ds <- x$degree.search
  if (!is.list(ds) || is.null(ds$restart.results) ||
      !is.list(ds$restart.results))
    return(list(hits = NA_real_, lookups = NA_real_))

  hits <- 0
  lookups <- 0
  found <- FALSE

  for (res in ds$restart.results) {
    native <- if (is.list(res)) res$native else NULL
    if (!is.list(native))
      next
    h <- .np_summary_named_number(native, "cache_hits")
    l <- .np_summary_named_number(native, "total_evaluations")
    if (.np_summary_missing_number(h) || .np_summary_missing_number(l))
      next
    hits <- hits + h
    lookups <- lookups + l
    found <- TRUE
  }

  if (!isTRUE(found))
    return(list(hits = NA_real_, lookups = NA_real_))
  list(hits = hits, lookups = lookups)
}

.np_summary_evaluation_cache_counts <- function(x) {
  ds <- x$degree.search
  nn.cache <- if (is.list(ds) && !is.null(ds$nn.cache)) ds$nn.cache else x$nn.cache
  objective.hits <- .np_summary_named_number(nn.cache, "objective.hits")
  objective.lookups <- .np_summary_named_number(nn.cache, "objective.visits")
  objective.raw <- .np_summary_named_number(nn.cache, "objective.raw.evals")
  if (.np_summary_missing_number(objective.hits))
    objective.hits <- .np_summary_named_number(nn.cache, "hits")
  if (.np_summary_missing_number(objective.lookups))
    objective.lookups <- .np_summary_named_number(nn.cache, "visits")
  if (.np_summary_missing_number(objective.lookups)) {
    if (.np_summary_missing_number(objective.raw))
      objective.raw <- .np_summary_named_number(nn.cache, "raw.evals")
    if (!.np_summary_missing_number(objective.raw) &&
        !.np_summary_missing_number(objective.hits))
      objective.lookups <- objective.raw + objective.hits
  }
  list(hits = objective.hits, lookups = objective.lookups)
}

.np_summary_r_optim_method <- function(x) {
  method <- x$pomethod
  if (is.null(method) || !length(method))
    return(NA_character_)
  method <- as.character(method[1L])
  if (is.na(method) || !(method %in% c("Nelder-Mead", "BFGS", "CG")))
    return(NA_character_)
  method
}

.np_summary_evaluation_cache_stage <- function(x) {
  optim.method <- .np_summary_r_optim_method(x)
  if (!is.na(optim.method))
    return(paste0("R optim: ", optim.method))

  route <- character()
  for (field in c("bwsolver", "search.engine")) {
    value <- x[[field]]
    if (!is.null(value) && length(value) && !is.na(value[1L]))
      route <- c(route, as.character(value[1L]))
  }
  degree.search <- x$degree.search
  if (is.list(degree.search)) {
    for (field in c("mode", "method", "engine")) {
      value <- degree.search[[field]]
      if (!is.null(value) && length(value) && !is.na(value[1L]))
        route <- c(route, as.character(value[1L]))
    }
  }

  if (any(route %in% c("powell", "mads+powell")))
    return("Powell")

  semiparametric <- inherits(x, "scbandwidth") || inherits(x, "sibandwidth")
  if (any(route %in% c("cell", "nomad+powell")))
    return(if (semiparametric) NA_character_ else "Powell")
  if (any(route %in% c("mads", "nomad")))
    return("NOMAD")

  core.classes <- c("bandwidth", "dbandwidth", "rbandwidth",
                    "conbandwidth", "condbandwidth", "plbandwidth")
  if (any(vapply(core.classes, function(class.name) inherits(x, class.name),
                 logical(1L))))
    return("Powell")

  NA_character_
}

.np_bandwidth_eval_accounting_lines <- function(x) {
  lines <- character()
  cache.hits.displayed <- 0

  nomad.cache <- .np_summary_nomad_cache_counts(x)
  line <- .np_summary_evaluation_cache_line("NOMAD",
                                             nomad.cache$hits,
                                             nomad.cache$lookups)
  has.nomad.cache <- !is.null(line)
  if (has.nomad.cache)
    lines <- c(lines, line)

  ## NOMAD cache hits are whole-parameter point lookups. The package-side
  ## evaluation-cache and fast-route counts below belong to different layers,
  ## so they are not subtracted from the NOMAD counts.
  evaluation.cache <- .np_summary_evaluation_cache_counts(x)
  line <- .np_summary_evaluation_cache_line(
    .np_summary_evaluation_cache_stage(x),
    evaluation.cache$hits,
    evaluation.cache$lookups
  )
  if (!is.null(line)) {
    lines <- c(lines, line)
    cache.hits.displayed <- cache.hits.displayed + evaluation.cache$hits
  }

  if (!(is.null(x$num.feval.fast) ||
        (length(x$num.feval.fast) == 1L && is.na(x$num.feval.fast)))) {
    fast <- suppressWarnings(as.numeric(x$num.feval.fast[1L]))
    computed <- suppressWarnings(as.numeric(x$num.feval[1L]))
    if (is.finite(fast) && fast > 0) {
      fast.extra <- if (cache.hits.displayed > 0) {
        fast - cache.hits.displayed
      } else {
        fast
      }
      if (is.finite(fast.extra) && fast.extra > 0) {
        line <- .np_summary_fast_route_line(fast.extra, computed)
        if (!is.null(line))
          lines <- c(lines, line)
      }
    }
  }

  if (!(is.null(x$num.feval.guarded) ||
        (length(x$num.feval.guarded) == 1L && is.na(x$num.feval.guarded)))) {
    guarded <- suppressWarnings(as.numeric(x$num.feval.guarded[1L]))
    if (is.finite(guarded) && guarded > 0)
      lines <- c(lines, paste0("Guarded evaluations: ",
                               .np_summary_format_count(guarded)))
  }

  lines
}


## bandwidth-related report generating functions
genBwSelStr <- function(x){
  fval.str <- ""
  if (!identical(x$fval, NA)) {
    ms.label <- .np_multistart_label(x)
    fval.str <- if (is.na(ms.label)) {
      paste("\nObjective Function Value: ", format(x$fval), sep = "")
    } else {
      paste("\nObjective Function Value: ", format(x$fval),
            " (achieved on multistart ", ms.label, ")", sep = "")
    }
  }

  nfe.str <- ""
  if(!(is.null(x$num.feval) || (length(x$num.feval) == 1L && is.na(x$num.feval)))){
    nfe.str <- paste("\nNumber of Function Evaluations: ", format(x$num.feval), sep="")
    nfe.lines <- .np_bandwidth_eval_accounting_lines(x)
    if(length(nfe.lines) > 0L)
      nfe.str <- paste(c(nfe.str, nfe.lines), collapse = "\n")
  }

  pregtype <- npFormatRegressionType(x)

  pregtype.str <- if (is.null(pregtype)) "" else paste("\n", npPolynomialSummaryLabel(x), ": ", pregtype, sep = "")
  pmethod.str <- if (is.null(x$pmethod)) "" else paste("\nBandwidth Selection Method:", x$pmethod)
  degree.search.str <- .np_degree_search_summary_str(x)
  formula.str <- if (!identical(x$formula, NULL)) paste("\nFormula:", paste(deparse(x$formula), collapse = "\n")) else ""
  ptype.str <- if (is.null(x$ptype)) "" else paste("\nBandwidth Type: ", x$ptype, sep = "")
  extendednn.str <- if (npRegressionHasExtendedNn(x)) {
    "\nExtended NN: K above n-1 scales the saturated nearest-neighbor bandwidth"
  } else {
    ""
  }

  paste(pregtype.str,
        pmethod.str,
        degree.search.str,
        formula.str,
        ptype.str,
        extendednn.str,
        fval.str,
        nfe.str,
        sep = "")
}

genBwScaleStrs <- function(x){
  ## approach is to take metadata and flatten it so it then can be
  ## processed into a single string

  vari <- names(x$sumNum)

  flat_varnames <- unlist(x$varnames[vari], use.names = FALSE)
  flat_vartitles <- unlist(
    lapply(vari, function(v) rep.int(x$vartitleabb[[v]], length(x$varnames[[v]]))),
    use.names = FALSE
  )
  flat_bandwidth <- unlist(x$bandwidth[vari], use.names = FALSE)
  flat_sum <- unlist(x$sumNum[vari], use.names = FALSE)
  flat_icon <- unlist(lapply(vari, function(v) x$dati[[v]]$icon), use.names = FALSE)

  sumText <- lapply(seq_along(flat_varnames), function(i) {
    if (isTRUE(flat_icon[[i]])) {
      if (x$type == "fixed") "Scale Factor:" else "NN Max:"
    } else {
      "Lambda Max:"
    }
  })

  valueText <- lapply(seq_along(flat_varnames), function(i) {
    if (isTRUE(flat_icon[[i]])) {
      if (x$type == "fixed") "Bandwidth:" else "NN Index:"
    } else {
      "Lambda:"
    }
  })

  nn.sample.max <- .np_nn_sample_max(x)
  flat_sum_display <- flat_sum
  if (!identical(x$type, "fixed") && is.finite(nn.sample.max)) {
    for (i in seq_along(flat_sum_display)) {
      if (isTRUE(flat_icon[[i]]))
        flat_sum_display[[i]] <- nn.sample.max
    }
  }

  print.sumText <- lapply(seq_along(sumText), function(i) {
    nzchar(sumText[[i]]) &&
      !(isTRUE(flat_icon[[i]]) && !identical(x$type, "fixed") && !is.finite(nn.sample.max))
  })

  bandwidth.display <- trimws(npFormat(flat_bandwidth), which = "both")
  sum.display <- trimws(npFormat(flat_sum_display), which = "both")
  value.label.width <- max(nchar(unlist(valueText)), na.rm = TRUE)
  sum.label.width <- max(nchar(unlist(sumText)), na.rm = TRUE)
  value.width <- max(nchar(bandwidth.display), na.rm = TRUE)
  sum.width <- max(nchar(sum.display), na.rm = TRUE)

  t.nchar <- lapply(flat_varnames, nchar)
                    
  maxNameLen <- max(unlist(t.nchar))

  vatText <- lapply(seq_along(t.nchar), function(j){
    paste("\n", rpad(flat_vartitles[[j]]), "Var. Name: ",
          flat_varnames[[j]],
          sapply(t.nchar[[j]], function(nc){
            paste(rep(' ', maxNameLen - nc), collapse='')
          }), sep='')
  })

  return(sapply(seq_along(t.nchar), function(j){
    sum.str <- ""
    if (isTRUE(print.sumText[[j]]))
      sum.str <- paste(.np_text_rpad(sumText[[j]], sum.label.width),
                       .np_text_lpad(sum.display[[j]], sum.width),
                       sep = " ")
    if (isTRUE(flat_icon[[j]]) &&
        !identical(x$type, "fixed") &&
        is.finite(nn.sample.max) &&
        is.finite(flat_bandwidth[[j]]) &&
        flat_bandwidth[[j]] > nn.sample.max)
      sum.str <- paste(sum.str, " Extended: yes", sep = "")
    value.str <- paste(.np_text_rpad(valueText[[j]], value.label.width),
                       .np_text_lpad(bandwidth.display[[j]], value.width),
                       sep = " ")
    paste(vatText[[j]], " ", value.str,
          if (nzchar(sum.str)) paste0("  ", sum.str) else "",
          sep = "", collapse = "")
  }))
}

npFormat <- function(x){
  format(sapply(x,format))
}

genBwKerStrs <- function(x){
  vari <- names(x$klist)

  ncon <- sapply(vari, function(v){
    sum(x$dati[[v]]$icon)
  })

  nuno <- sapply(vari, function(v){
    sum(x$dati[[v]]$iuno)
  })

  nord <- sapply(vari, function(v){
    sum(x$dati[[v]]$iord)
  })

  cktype <- sapply(vari, function(v){
    x$klist[[v]]$ckertype
  })

  uktype <- sapply(vari, function(v){
    x$klist[[v]]$ukertype
  })

  oktype <- sapply(vari, function(v){
    x$klist[[v]]$okertype
  })

  tt <- ''

  if(any(ncon > 0)){
    ctype.str <- ""
    if (length(unique(cktype)) == 1) {
      ctype.str <- paste("\nContinuous Kernel Type:",
                         x$klist[[vari[1]]]$pckertype)
    } else {
      ctype.str <- paste(sapply(seq_along(vari), function(v){
        if (ncon[v] > 0)
          paste("\nContinuous Kernel Type (",
                x$vartitleabb[[vari[v]]],
                " Var.): ", x$klist[[vari[v]]]$pckertype, sep = "")
        else ""
      }), collapse = "")
    }
    tt <- paste("\n", ctype.str, sep = "")
    cont.str <- paste(sapply(seq_along(vari), function(i){
      if (ncon[i] > 0)
        paste("\nNo. Continuous", pad(x$vartitle[[vari[i]]]), "Vars.: ",
              ncon[i], sep = "")
      else ""
    }), collapse = "")
    tt <- paste(tt, cont.str, sep = "")
  }
                
    
  if(any(nuno > 0)) {
    utype.str <- ""
    if (length(unique(uktype)) == 1) {
      utype.str <- paste("\nUnordered Categorical Kernel Type:",
                         x$klist[[vari[1]]]$pukertype)
    } else {
      utype.str <- paste(sapply(seq_along(vari), function(i){
        if (nuno[i] > 0)
          paste("\nUnordered Categorical Kernel Type (",
                x$vartitleabb[[vari[i]]],
                " Var.): ", x$klist[[vari[i]]]$pukertype, sep = "")
        else ""
      }), collapse = "")
    }
    tt <- paste(tt, "\n", utype.str, sep = "")
    uno.str <- paste(sapply(seq_along(vari), function(i){
      if (nuno[i] > 0)
        paste("\nNo. Unordered Categorical", pad(x$vartitle[[vari[i]]]), "Vars.: ",
              nuno[i], sep = "")
      else ""
    }), collapse = "")
    tt <- paste(tt, uno.str, sep = "")

  }

  if(any(nord > 0)) {
    otype.str <- ""
    if (length(unique(oktype)) == 1) {
      otype.str <- paste("\nOrdered Categorical Kernel Type:",
                         x$klist[[vari[1]]]$pokertype)
    } else {
      otype.str <- paste(sapply(seq_along(vari), function(i){
        if (nord[i] > 0)
          paste("\nOrdered Categorical Kernel Type (",
                x$vartitleabb[[vari[i]]],
                " Var.): ", x$klist[[vari[i]]]$pokertype, sep = "")
        else ""
      }), collapse = "")
    }
    tt <- paste(tt, "\n", otype.str, sep = "")
    ord.str <- paste(sapply(seq_along(vari), function(i){
      if (nord[i] > 0)
        paste("\nNo. Ordered Categorical", pad(x$vartitle[[vari[i]]]), "Vars.: ",
              nord[i], sep = "")
      else ""
    }), collapse = "")
    tt <- paste(tt, ord.str, sep = "")

  }

  return(tt)
}

genBwKerStrsXY <- function(x){
  t.str <- ''
  cnt <- 0
  
  if (x$xncon + x$yncon > 0){
    cker.str <- ""
    if (x$pcxkertype == x$pcykertype) {
      cker.str <- paste("\n\nContinuous Kernel Type:", x$pcxkertype)
    } else {
      exp.str <- if (x$xncon > 0)
        paste("\nContinuous Kernel Type (Exp. Var.):", x$pcxkertype)
      else ""
      dep.str <- if (x$yncon > 0)
        paste("\nContinuous Kernel Type (Dep. Var.):", x$pcykertype)
      else ""
      cker.str <- paste("\n", exp.str, dep.str)
    }
    ycon.str <- if (x$yncon > 0) paste("\nNo. Continuous Dependent Vars.:", x$yncon) else ""
    xcon.str <- if (x$xncon > 0) paste("\nNo. Continuous Explanatory Vars.:", x$xncon) else ""
    t.str[cnt <- cnt + 1] <-
      paste(cker.str, ycon.str, xcon.str)
  }

  if (x$xnuno + x$ynuno > 0){
    uker.str <- ""
    if (x$puxkertype == x$puykertype) {
      uker.str <- paste("\n\nUnordered Categorical Kernel Type:", x$puxkertype)
    } else {
      exp.str <- if (x$xnuno > 0)
        paste("\nUnordered Categorical Kernel Type (Exp. Var.):", x$puxkertype)
      else ""
      dep.str <- if (x$ynuno > 0)
        paste("\nUnordered Categorical Kernel Type (Dep. Var.):", x$puykertype)
      else ""
      uker.str <- paste("\n", exp.str, dep.str)
    }
    yuno.str <- if (x$ynuno > 0) paste("\nNo. Unordered Categorical Dependent Vars.:", x$ynuno) else ""
    xuno.str <- if (x$xnuno > 0) paste("\nNo. Unordered Categorical Explanatory Vars.:", x$xnuno) else ""
    t.str[cnt <- cnt + 1] <-
      paste(uker.str, yuno.str, xuno.str)
  }

  if (x$xnord + x$ynord > 0){
    oker.str <- ""
    if (x$poxkertype == x$poykertype) {
      oker.str <- paste("\n\nOrdered Categorical Kernel Type:", x$poxkertype)
    } else {
      exp.str <- if (x$xnord > 0)
        paste("\nOrdered Categorical Kernel Type (Exp. Var.):", x$poxkertype)
      else ""
      dep.str <- if (x$ynord > 0)
        paste("\nOrdered Categorical Kernel Type (Dep. Var.):", x$poykertype)
      else ""
      oker.str <- paste("\n", exp.str, dep.str)
    }
    yord.str <- if (x$ynord > 0) paste("\nNo. Ordered Categorical Dependent Vars.:", x$ynord) else ""
    xord.str <- if (x$xnord > 0) paste("\nNo. Ordered Categorical Explanatory Vars.:", x$xnord) else ""
    t.str[cnt <- cnt + 1] <-
      paste(oker.str, yord.str, xord.str)
  }
  return(t.str)
}

genBwGOFStrs <- function(x) {
  ###paste("Residual standard error:", sqrt
}
  
## statistical functions
RSQfunc <- function(y,y.pred) {
  y.mean <- mean(y)
  (sum((y-y.mean)*(y.pred-y.mean))^2)/(sum((y-y.mean)^2)*sum((y.pred-y.mean)^2))
}

MSEfunc <- function(y,y.fit) {
  mean((y-y.fit)^2)
}

MAEfunc <- function(y,y.fit) {
  mean(abs(y-y.fit))
}

MAPEfunc <- function(y,y.fit) {
  jj = which(y != 0)
  
  mean(c(abs((y[jj]-y.fit[jj])/y[jj]), as.numeric(replicate(length(y)-length(jj),2))))
}

CORRfunc <- function(y,y.fit) {
  abs(corr(cbind(y,y.fit)))
}

SIGNfunc <- function(y,y.fit) {
  sum(sign(y) == sign(y.fit))/length(y)
}

EssDee <- function(y){
  if(any(dim(as.matrix(y)) == 0))
      return(0)
  sd.vec <- apply(as.matrix(y),2,sd)
  IQR.vec <- apply(as.matrix(y),2,IQR)/QFAC
  mad.vec <- apply(as.matrix(y),2,mad)
  a <- apply(cbind(sd.vec,IQR.vec,mad.vec),1, function(x) max(x))
  if(any(a<=0)) .np_warning(paste("variable ",which(a<=0)," appears to be constant",sep=""))
  a <- apply(cbind(sd.vec,IQR.vec,mad.vec),1, function(x) min(x[x>0]))  
  return(a)
}

.npConditionalNomadContScale <- function(ycon, xcon, iycon, ixcon, nconfac, where) {
  y_scale <- if (length(which(iycon)) > 0L) EssDee(ycon) else numeric(0L)
  x_scale <- if (length(which(ixcon)) > 0L) EssDee(xcon) else numeric(0L)
  cont_scale <- c(y_scale, x_scale) * nconfac
  expected <- length(which(iycon)) + length(which(ixcon))

  if (length(cont_scale) != expected) {
    stop(sprintf(
      "%s: internal NOMAD bandwidth scale mismatch (%d scale values for %d continuous bandwidth slots)",
      where, length(cont_scale), expected
    ))
  }

  cont_scale
}

.npAssertConditionalNomadSetup <- function(setup, where) {
  if (length(setup$cont_scale) != length(setup$cont_flat)) {
    stop(sprintf(
      "%s: internal NOMAD bandwidth transform mismatch (%d continuous scales for %d continuous slots)",
      where, length(setup$cont_scale), length(setup$cont_flat)
    ))
  }
  if (length(setup$cat_upper) != length(setup$cat_flat)) {
    stop(sprintf(
      "%s: internal NOMAD bandwidth transform mismatch (%d categorical bounds for %d categorical slots)",
      where, length(setup$cat_upper), length(setup$cat_flat)
    ))
  }
  invisible(setup)
}


##EssDee <- function(y){
##
##  sd.vec <- apply(as.matrix(y),2,sd)
##  IQR.vec <- apply(as.matrix(y),2,IQR)/QFAC
##  return(ifelse(sd.vec<IQR.vec|IQR.vec==0,sd.vec,IQR.vec))
##  
##}

## consolidating various bits of code related to converting internal settings
## to printable strings

bwmToPrint <- function(s){
  switch(s,
         manual = "Manual",
         cv.aic = "Expected Kullback-Leibler Cross-Validation",
         cv.ml = "Maximum Likelihood Cross-Validation",
         cv.cdf = "Least Squares Cross-Validation",
         cv.ls = "Least Squares Cross-Validation",
         "normal-reference" = "Normal Reference")
}

bwtToPrint <- function(s){
  switch(s,
         fixed = "Fixed",
         generalized_nn = "Generalized Nearest Neighbour",
         adaptive_nn = "Adaptive Nearest Neighbour" )
}

cktToPrint <- function(s, order = "", kerbound = "none"){
  pck <- switch(s,
                gaussian = paste(order,"Gaussian"),
                epanechnikov =  paste(order,"Epanechnikov"),
                uniform = "Uniform",
                "truncated gaussian" = "Truncated Gaussian")
  if (!is.null(kerbound) && !identical(kerbound, "none"))
    pck <- paste0(pck, " (bounded/", kerbound, ")")
  pck
}

uktToPrint <- function(s){
  switch(s,
         aitchisonaitken = "Aitchison and Aitken",
         liracine = "Li and Racine (normalized)")
}

oktToPrint <- function(s, normalized = FALSE) {
  if(normalized){
    pok <- 
      switch(s,
             wangvanryzin = "Wang and Van Ryzin", 
             liracine = "Li and Racine (normalized)",
             nliracine = "Li and Racine (normalized)",
             "racineliyan" = "Racine, Li, and Yan")
  } else {
    pok <- 
      switch(s,
             wangvanryzin = "Wang and Van Ryzin", 
             liracine = "Li and Racine",
             nliracine = "Li and Racine (normalized)",
             "racineliyan" = "Racine, Li, and Yan")
  }
  return(pok)
}

### holding place for some generic methods

se <- function(x){
  UseMethod("se",x)
}

gradients <- function(x, ...){
  UseMethod("gradients",x)
}

.np_reject_unused_dots <- function(dots, where) {
  if (length(dots) == 0L)
    return(invisible(TRUE))

  dot.names <- names(dots)
  if (is.null(dot.names))
    dot.names <- rep("", length(dots))

  labels <- ifelse(nzchar(dot.names), dot.names, "<unnamed>")
  labels <- paste(sprintf("'%s'", labels), collapse = ", ")
  stop(sprintf("unused argument%s in %s: %s",
               if (length(dots) == 1L) "" else "s",
               where,
               labels),
       call. = FALSE)
}

## From crs to avoid crs:::W.glp

mypoly <- function(x,
                   ex=NULL,
                   degree,
                   gradient.compute = FALSE,
                   r=0,
                   Bernstein = TRUE) {

  if(missing(x)) stop(" Error: x required")
  if(missing(degree)) stop(" Error: degree required")
  if(degree < 1) stop(" Error: degree must be a positive integer")
  if(!is.logical(Bernstein)) stop(" Error: Bernstein must be logical")

  if(!Bernstein) {

    ## Raw polynomials and their derivatives

    if(!is.null(ex)) x <- ex

    if(gradient.compute) {
      Z <- NULL
      for(i in seq_len(degree)) {
        if((i-r) >= 0) {
          tmp <- (factorial(i)/factorial(i-r))*x^max(0,i-r)
        } else {
          tmp <- rep(0,length(x))
        }
        Z <- cbind(Z,tmp)
      }
    } else {
      Z <- outer(x,1L:degree,"^")
    }

  } else {
    ## Bernstein polynomial basis (no interior knots) aligned with C-side LP code.
    x.train <- as.double(x)
    x.use <- if (is.null(ex)) x.train else as.double(ex)
    r <- if (gradient.compute) as.integer(r) else 0L
    if (is.na(r) || r < 0L) stop(" Error: derivative order must be a non-negative integer")

    xmin <- min(x.train)
    xmax <- max(x.train)
    xrange <- xmax - xmin
    if (!is.finite(xrange) || xrange <= 0) {
      Z <- matrix(0.0, nrow = length(x.use), ncol = degree)
    } else if (r > degree) {
      Z <- matrix(0.0, nrow = length(x.use), ncol = degree)
    } else {
      u <- (x.use - xmin) / xrange
      m <- degree - r
      coef.deriv <- if (r == 0L) {
        1.0
      } else {
        exp(lfactorial(degree) - lfactorial(m)) / (xrange^r)
      }

      z.list <- lapply(seq_len(degree), function(idx) {
        out <- rep(0.0, length(u))
        for (k in 0:r) {
          j <- idx - k
          if (j >= 0L && j <= m) {
            out <- out + ((-1.0)^(r - k)) * choose(r, k) * choose(m, j) * (u^j) * ((1.0 - u)^(m - j))
          }
        }
        coef.deriv * out
      })
      Z <- do.call(cbind, z.list)
    }
  }

  return(as.matrix(Z))

}

## W.lp accepts a vector of degrees and provides local-polynomial bases
## with selectable term structure (glp/additive/tensor).

npBuildLpTerms <- function(degree, basis = c("glp", "additive", "tensor")) {
  basis <- match.arg(basis)
  k <- length(degree)
  if (k == 0L)
    return(matrix(integer(0), nrow = 1L, ncol = 0L))

  degree <- as.integer(degree)
  degree.list <- lapply(degree, function(d) 0:d)
  z <- as.matrix(do.call(base::expand.grid, degree.list))
  s <- rowSums(z)

  if (identical(basis, "glp")) {
    ind <- (s > 0) & (s <= max(degree))
    z <- z[ind, , drop = FALSE]
    if (!all(degree == max(degree))) {
      for (j in seq_along(degree)) {
        d <- degree[j]
        if ((d < max(degree)) && (d > 0)) {
          s <- rowSums(z)
          dropj <- (s > d) & (z[, j, drop = FALSE] == matrix(d, nrow(z), 1, byrow = TRUE))
          z <- z[!dropj, , drop = FALSE]
        }
      }
    }
  } else if (identical(basis, "additive")) {
    ind <- (s > 0) & (rowSums(z > 0) == 1L)
    z <- z[ind, , drop = FALSE]
  } else if (identical(basis, "tensor")) {
    z <- z[s > 0, , drop = FALSE]
  }

  rbind(matrix(0L, nrow = 1L, ncol = k), z)
}

W.lp <- function(xdat = NULL,
                 exdat = NULL,
                 degree = NULL,
                 gradient.vec = NULL,
                 basis = c("glp", "additive", "tensor"),
                 bernstein.basis = TRUE,
                 Bernstein = bernstein.basis) {

  if(is.null(xdat)) stop(" Error: You must provide data")
  if(is.null(degree) || any(degree < 0)) stop(paste(" Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))
  basis <- match.arg(basis)

  xdat <- as.data.frame(xdat)

  xdat.col.numeric <- vapply(seq_len(ncol(xdat)), function(i) is.numeric(xdat[,i]), logical(1))
  k <- ncol(as.data.frame(xdat[,xdat.col.numeric]))

  xdat.numeric <- NULL
  if(k > 0) {
    xdat.numeric <- as.data.frame(xdat[,xdat.col.numeric])
    if(!is.null(exdat)) {
      exdat.numeric <- as.data.frame(exdat[,xdat.col.numeric])
    } else {
      exdat.numeric <- NULL
    }
  }

  if(!is.null(gradient.vec) && (length(gradient.vec) != k)) stop(paste(" Error: gradient vector and number of numeric predictors must be conformable\n",sep=""))
  if(!is.null(gradient.vec) && any(gradient.vec < 0)) stop(paste(" Error: gradient vector must contain non-negative integers\n",sep=""))

  if(!is.null(gradient.vec)) {
    gradient.compute <- TRUE
  } else {
    gradient.compute <- FALSE
    gradient.vec <- rep(NA,k)
  }

  if(length(degree) != k) stop(" Error: degree vector and number of numeric predictors incompatible")

  if(all(degree == 0) || (k == 0)) {

    ## Local constant OR no continuous variables

    if(is.null(exdat)) {
      return(matrix(1,nrow=nrow(as.data.frame(xdat)),ncol=1))
    } else {
      return(matrix(1,nrow=nrow(as.data.frame(exdat)),ncol=1))
    }

  } else {

    z <- npBuildLpTerms(degree = degree, basis = basis)
    z.noi <- z[-1L, , drop = FALSE]
    if(is.null(exdat)) {
      res <- rep.int(1,nrow(xdat.numeric))
    } else {
      res <- rep.int(1,nrow(exdat.numeric))
    }
    res.deriv <- 1
    if(degree[1] > 0) {
      res <- cbind(1, mypoly(x=xdat.numeric[,1],
                             ex=exdat.numeric[,1],
                             degree=degree[1],
                             gradient.compute=gradient.compute,
                             r=gradient.vec[1],
                             Bernstein=Bernstein))[, 1 + z.noi[, 1]]

      if(gradient.compute && gradient.vec[1] != 0) res.deriv <- cbind(1,matrix(NA,1,degree[1]))[, 1 + z.noi[, 1],drop=FALSE]
      if(gradient.compute && gradient.vec[1] == 0) res.deriv <- cbind(1,matrix(0,1,degree[1]))[, 1 + z.noi[, 1],drop=FALSE]
    }
    if(k > 1) for (i in 2:k) if(degree[i] > 0) {
      res <- res * cbind(1, mypoly(x=xdat.numeric[,i],
                                   ex=exdat.numeric[,i],
                                   degree=degree[i],
                                   gradient.compute=gradient.compute,
                                   r=gradient.vec[i],
                                   Bernstein=Bernstein))[, 1 + z.noi[, i]]
      if(gradient.compute && gradient.vec[i] != 0) res.deriv <- res.deriv * cbind(1,matrix(NA,1,degree[i]))[, 1 + z.noi[, i],drop=FALSE]
      if(gradient.compute && gradient.vec[i] == 0) res.deriv <- res.deriv *cbind(1,matrix(0,1,degree[i]))[, 1 + z.noi[, i],drop=FALSE]
    }

    if(is.null(exdat)) {
      res <- matrix(res,nrow=NROW(xdat))
    } else {
      res <- matrix(res,nrow=NROW(exdat))
    }
    if(gradient.compute) res.deriv <- matrix(res.deriv,nrow=1)
    colnames(res) <- apply(z.noi, 1L, function(x) paste(x, collapse = "."))
    if(gradient.compute) colnames(res.deriv) <- apply(z.noi, 1L, function(x) paste(x, collapse = "."))

    if(gradient.compute) {
      res[,!is.na(as.numeric(res.deriv))] <- 0
      return(cbind(0,res))
    } else {
      return(cbind(1,res))
    }

  }

}

### internal constants used in the c backend

SF_NORMAL = 0
SF_ARB = 1

BW_FIXED = 0
BW_GEN_NN = 1
BW_ADAP_NN = 2

IMULTI_TRUE = 1
IMULTI_FALSE = 0

RE_MIN_TRUE = 0
RE_MIN_FALSE = 1

IO_MIN_TRUE = 1
IO_MIN_FALSE = 0

USE_START_NO = 0
USE_START_YES = 1

NP_DO_DENS = 1
NP_DO_DIST = 0

## initially making an np-wide option via the 'options' mechanism
DO_TREE_NO = 0
DO_TREE_YES = 1

##kernel defs
CKER_GAUSS = 0
CKER_EPAN  = 4
CKER_UNI   = 8
CKER_TGAUSS = 9

UKER_AIT = 0
UKER_LR = 1

OKER_WANG = 0
OKER_LR = 1
OKER_NLR = 2
OKER_RLY = 3

##density 
BWM_CVML = 0
BWM_CVLS = 1
BWM_CVML_NP= 2

##distribution
DBWM_CVLS = 0

##regression
BWM_CVAIC = 0
RBWM_CVKS = 3L

REGTYPE_LC = 0
REGTYPE_LL = 1
REGTYPE_LP = 2
# legacy alias retained for internal/backward compatibility
REGTYPE_GLP = REGTYPE_LP

##conditional density/distribution
CBWM_CVML = 0
CBWM_CVLS = 1
CBWM_NPLS = 2
CBWM_CCDF = 3 # Added 7/2/2010 jracine

##conditional distribution
CDBWM_CVLS = 0

##integral operators on kernels
OP_NOOP        = -1
OP_NORMAL      = 0
OP_CONVOLUTION = 1
OP_DERIVATIVE  = 2
OP_INTEGRAL    = 3

ALL_OPERATORS = c(OP_NORMAL, OP_CONVOLUTION, OP_DERIVATIVE, OP_INTEGRAL)
names(ALL_OPERATORS) <- c("normal","convolution", "derivative", "integral")

PERMUTATION_OPERATORS <- c(OP_NOOP, OP_NORMAL, OP_DERIVATIVE, OP_INTEGRAL)
names(PERMUTATION_OPERATORS) <- c("none", "normal", "derivative", "integral")

## useful numerical constants of kernel integrals
int.kernels <- c(0.28209479177387814348, 0.47603496111841936711, 0.62396943688265038571, 0.74785078617543927990,
                 0.26832815729997476357, 0.55901699437494742410, 0.84658823667359826246, 1.1329342579014329689,
                 0.5, 2.90113075268188e-01)

QFAC <- qnorm(.25,lower.tail=FALSE)*2

.np_eval_bws_call_arg <- function(bws, arg) {
  if (is.null(bws$call))
    stop("bandwidth object does not contain a call component")

  expr <- bws$call[[arg]]
  if (is.null(expr))
    stop(sprintf("bandwidth call does not contain '%s'", arg))

  eval.env <- environment(bws$call)
  if (is.null(eval.env))
    eval.env <- parent.frame()

  if (!is.language(expr))
    return(expr)

  val <- .np_try_eval_in_frames(expr, eval_env = eval.env, search_frames = FALSE)
  if (isTRUE(val$ok))
    return(val$value)

  fallback <- bws[[arg]]
  if (!is.null(fallback))
    return(fallback)

  if (inherits(val$error, "error"))
    stop(conditionMessage(val$error), call. = FALSE)
  stop(sprintf("unable to evaluate '%s' in bandwidth call", arg), call. = FALSE)
}

.np_eval_call_arg <- function(call_obj, arg, caller_env = parent.frame()) {
  if (is.null(call_obj))
    stop("object does not contain a call component")

  expr <- call_obj[[arg]]
  if (is.null(expr))
    stop(sprintf("call does not contain '%s'", arg))

  eval.env <- environment(call_obj)
  if (is.null(eval.env))
    eval.env <- caller_env

  if (!is.language(expr))
    return(expr)

  if (is.symbol(expr)) {
    val <- .np_try_eval_in_frames(expr, eval_env = eval.env, search_frames = FALSE)
    if (isTRUE(val$ok))
      return(val$value)

    if (is.environment(caller_env)) {
      not_found <- .np_missing_binding_sentinel
      sym_val <- get0(as.character(expr), envir = caller_env, inherits = TRUE, ifnotfound = not_found)
      if (!identical(sym_val, not_found))
        return(sym_val)
    }
  }

  val <- .np_try_eval_in_frames(expr, eval_env = eval.env)
  if (isTRUE(val$ok))
    return(val$value)

  if (inherits(val$error, "error"))
    stop(conditionMessage(val$error), call. = FALSE)
  stop(sprintf("unable to evaluate call argument '%s'", arg), call. = FALSE)
}

.np_nn_cache_stats <- function(x) {
  if (is.null(x))
    return(NULL)
  x <- as.numeric(x)
  nms <- c("enabled", "key.length", "visits", "unique", "repeats",
           "raw.evals", "hits", "allocation.failed",
           "objective.enabled", "objective.key.length", "objective.visits",
           "objective.unique", "objective.repeats", "objective.raw.evals",
           "objective.hits", "objective.allocation.failed")
  names(x) <- nms[seq_along(x)]
  x
}

.np_r_nn_cache_new <- function(enabled, key.length = 0L) {
  cache <- new.env(parent = emptyenv(), hash = FALSE)
  cache$enabled <- isTRUE(enabled)
  cache$key.length <- if (isTRUE(cache$enabled)) as.integer(key.length) else 0L
  cache$store <- new.env(parent = emptyenv(), hash = TRUE)
  cache$visits <- 0
  cache$unique <- 0
  cache$repeats <- 0
  cache$raw.evals <- 0
  cache$hits <- 0
  cache$allocation.failed <- 0
  cache
}

.np_r_nn_cache_key <- function(key) {
  paste(as.integer(key), collapse = "\r")
}

.np_r_nn_cache_param_key <- function(doubles = numeric(0), integers = integer(0)) {
  paste(c(sprintf("%.17g", as.double(doubles)), as.integer(integers)), collapse = "\r")
}

.np_r_nn_cache_get_token <- function(cache, token) {
  if (!is.environment(cache) || !isTRUE(cache$enabled))
    return(list(hit = FALSE, token = NULL, value = NULL))
  cache$visits <- cache$visits + 1
  if (exists(token, envir = cache$store, inherits = FALSE)) {
    cache$repeats <- cache$repeats + 1
    cache$hits <- cache$hits + 1
    return(list(
      hit = TRUE,
      token = token,
      value = get(token, envir = cache$store, inherits = FALSE)
    ))
  }
  cache$unique <- cache$unique + 1
  list(hit = FALSE, token = token, value = NULL)
}

.np_r_nn_cache_get <- function(cache, key) {
  .np_r_nn_cache_get_token(cache, .np_r_nn_cache_key(key))
}

.np_r_nn_cache_put <- function(cache, token, value) {
  if (!is.environment(cache) || !isTRUE(cache$enabled) || is.null(token))
    return(invisible(FALSE))
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value))
    return(invisible(FALSE))
  assign(token, as.numeric(value), envir = cache$store)
  cache$raw.evals <- cache$raw.evals + 1
  invisible(TRUE)
}

.np_r_nn_cache_stats <- function(cache) {
  if (!is.environment(cache))
    return(NULL)
  c(
    enabled = if (isTRUE(cache$enabled)) 1 else 0,
    key.length = as.numeric(cache$key.length),
    visits = cache$visits,
    unique = cache$unique,
    repeats = cache$repeats,
    raw.evals = cache$raw.evals,
    hits = cache$hits,
    allocation.failed = cache$allocation.failed
  )
}

.np_r_nn_cache_combine_stats <- function(stats) {
  stats <- Filter(Negate(is.null), stats)
  if (!length(stats))
    return(NULL)
  nms <- c("enabled", "key.length", "visits", "unique", "repeats",
           "raw.evals", "hits", "allocation.failed",
           "objective.enabled", "objective.key.length", "objective.visits",
           "objective.unique", "objective.repeats", "objective.raw.evals",
           "objective.hits", "objective.allocation.failed")
  mat <- do.call(rbind, lapply(stats, function(x) {
    x <- as.numeric(x)
    names(x) <- nms[seq_along(x)]
    x[nms]
  }))
  out <- colSums(mat, na.rm = TRUE)
  out["enabled"] <- as.numeric(any(mat[, "enabled"] > 0))
  out["key.length"] <- if (all(is.na(mat[, "key.length"]))) {
    0
  } else {
    max(mat[, "key.length"], na.rm = TRUE)
  }
  if ("objective.enabled" %in% names(out))
    out["objective.enabled"] <- as.numeric(any(mat[, "objective.enabled"] > 0, na.rm = TRUE))
  if ("objective.key.length" %in% names(out))
    out["objective.key.length"] <- if (all(is.na(mat[, "objective.key.length"]))) {
      0
    } else {
      max(mat[, "objective.key.length"], na.rm = TRUE)
    }
  out
}

.np_objective_exact_cache_key <- function(x) {
  paste(sprintf("%a", as.double(x)), collapse = "\r")
}

.np_objective_exact_cache_new <- function(enabled) {
  cache <- new.env(parent = emptyenv(), hash = FALSE)
  cache$enabled <- isTRUE(enabled)
  cache$store <- new.env(parent = emptyenv(), hash = TRUE)
  cache$visits <- 0
  cache$unique <- 0
  cache$repeats <- 0
  cache$raw.evals <- 0
  cache$hits <- 0
  cache
}

.np_objective_exact_cache_get <- function(cache, x) {
  if (!is.environment(cache) || !isTRUE(cache$enabled))
    return(list(hit = FALSE, token = NULL, value = NULL))
  token <- .np_objective_exact_cache_key(x)
  cache$visits <- cache$visits + 1
  if (exists(token, envir = cache$store, inherits = FALSE)) {
    cache$repeats <- cache$repeats + 1
    cache$hits <- cache$hits + 1
    return(list(
      hit = TRUE,
      token = token,
      value = get(token, envir = cache$store, inherits = FALSE)
    ))
  }
  cache$unique <- cache$unique + 1
  list(hit = FALSE, token = token, value = NULL)
}

.np_objective_exact_cache_put <- function(cache, token, value) {
  if (!is.environment(cache) || !isTRUE(cache$enabled) || is.null(token))
    return(invisible(FALSE))
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value))
    return(invisible(FALSE))
  assign(token, as.numeric(value), envir = cache$store)
  cache$raw.evals <- cache$raw.evals + 1
  invisible(TRUE)
}

.np_objective_exact_cache_stats <- function(cache) {
  if (!is.environment(cache))
    return(NULL)
  c(
    enabled = if (isTRUE(cache$enabled)) 1 else 0,
    visits = cache$visits,
    unique = cache$unique,
    repeats = cache$repeats,
    raw.evals = cache$raw.evals,
    hits = cache$hits
  )
}
