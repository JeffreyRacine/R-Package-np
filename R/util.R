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

  if (is.null(degree)) {
    if (ncon == 0L)
      return(integer(0))
    stop(sprintf("%s must be supplied explicitly when regtype='lp'", argname))
  }

  if (!length(degree) && ncon == 0L)
    return(integer(0))

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

.np_prepare_nomad_shortcut <- function(nomad,
                                       call_names,
                                       preset,
                                       values,
                                       where = "npregbw") {
  call_names <- if (is.null(call_names)) character(0) else call_names[nzchar(call_names)]
  nomad.enabled <- if ("nomad" %in% call_names) {
    npValidateScalarLogical(nomad, "nomad")
  } else {
    FALSE
  }

  metadata <- list(
    enabled = isTRUE(nomad.enabled),
    where = where,
    preset = "lp_nomad",
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

.np_nomad_sum_time <- function(nomad.time = NULL,
                               powell.time = NULL,
                               fallback = NULL) {
  pieces <- c(
    if (!is.null(nomad.time) && length(nomad.time)) as.double(nomad.time[1L]) else NA_real_,
    if (!is.null(powell.time) && length(powell.time)) as.double(powell.time[1L]) else NA_real_
  )

  if (any(is.finite(pieces)))
    return(sum(pieces, na.rm = TRUE))

  if (!is.null(fallback) && length(fallback) && is.finite(as.double(fallback[1L])))
    return(as.double(fallback[1L]))

  NA_real_
}

.npRmpi_reconcile_nomad_search_timing <- function(search.result) {
  if (!is.list(search.result))
    return(search.result)

  optim.time <- .np_nomad_sum_time(
    nomad.time = search.result$nomad.time,
    powell.time = search.result$powell.time,
    fallback = search.result$optim.time
  )

  if (is.finite(optim.time))
    search.result$optim.time <- optim.time

  search.result
}

.npRmpi_reconcile_nomad_bws_timing <- function(bws) {
  if (!is.list(bws))
    return(bws)

  optim.time <- .np_nomad_sum_time(
    nomad.time = bws$nomad.time,
    powell.time = bws$powell.time,
    fallback = bws$total.time
  )

  if (is.finite(optim.time))
    bws$total.time <- optim.time

  if (is.list(bws$degree.search)) {
    degree.optim <- .np_nomad_sum_time(
      nomad.time = bws$degree.search$nomad.time,
      powell.time = bws$degree.search$powell.time,
      fallback = bws$degree.search$optim.time
    )
    if (is.finite(degree.optim))
      bws$degree.search$optim.time <- degree.optim
  }

  bws
}

.npRmpi_restore_nomad_fit_bws_metadata <- function(result, bws) {
  if (!is.list(result) || !is.list(result$bws) || !is.list(bws))
    return(result)

  bws <- .npRmpi_reconcile_nomad_bws_timing(bws)
  result$bws <- .npRmpi_reconcile_nomad_bws_timing(result$bws)

  # Fit objects should preserve the selected bandwidth object's telemetry and labels.
  bandwidth.metadata.fields <- c(
    "method",
    "pmethod",
    "fval",
    "ifval",
    "num.feval",
    "num.feval.fast",
    "fval.history",
    "eval.history",
    "invalid.history",
    "timing",
    "timing.profile",
    "degree.search",
    "nomad.shortcut"
  )
  for (field in bandwidth.metadata.fields) {
    if (!is.null(bws[[field]]))
      result$bws[[field]] <- bws[[field]]
  }

  if (!is.null(bws$nomad.time) && is.finite(bws$nomad.time))
    result$nomad.time <- as.double(bws$nomad.time)

  if (!is.null(bws$powell.time) && is.finite(bws$powell.time))
    result$powell.time <- as.double(bws$powell.time)

  optim.time <- .np_nomad_sum_time(
    nomad.time = result$nomad.time,
    powell.time = result$powell.time,
    fallback = bws$total.time
  )
  if (is.finite(optim.time))
    result$optim.time <- optim.time

  if (is.finite(result$optim.time) &&
      !is.null(result$fit.time) &&
      is.finite(result$fit.time))
    result$total.time <- as.double(result$optim.time) + as.double(result$fit.time)

  if ((is.null(result$bws$nomad.time) || !is.finite(result$bws$nomad.time)) &&
      !is.null(bws$nomad.time) && is.finite(bws$nomad.time))
    result$bws$nomad.time <- as.double(bws$nomad.time)

  if ((is.null(result$bws$powell.time) || !is.finite(result$bws$powell.time)) &&
      !is.null(bws$powell.time) && is.finite(bws$powell.time))
    result$bws$powell.time <- as.double(bws$powell.time)

  if ((is.null(result$bws$total.time) || !is.finite(result$bws$total.time)) &&
      !is.null(bws$total.time) && is.finite(bws$total.time))
    result$bws$total.time <- as.double(bws$total.time)

  result$bws <- .npRmpi_reconcile_nomad_bws_timing(result$bws)

  result
}

.np_nomad_validate_inner_multistart <- function(call_names = character(),
                                                dot.args = list(),
                                                nomad.nmulti = 0L,
                                                regtype,
                                                automatic.degree.search,
                                                search.engine) {
  if (is.null(call_names))
    call_names <- character()
  dot.names <- names(dot.args)
  if (is.null(dot.names))
    dot.names <- character()

  inner.named <- ("nomad.nmulti" %in% call_names) || ("nomad.nmulti" %in% dot.names)
  inner.raw <- if ("nomad.nmulti" %in% dot.names) dot.args[["nomad.nmulti"]] else nomad.nmulti
  inner.nmulti <- if (inner.named) {
    npValidateNonNegativeInteger(inner.raw, "nomad.nmulti")
  } else {
    0L
  }

  regtype.value <- if (length(regtype)) as.character(regtype)[1L] else ""
  search.engine.value <- if (length(search.engine)) as.character(search.engine)[1L] else ""

  if (inner.named &&
      (!identical(regtype.value, "lp") ||
       !isTRUE(automatic.degree.search) ||
       !(search.engine.value %in% c("nomad", "nomad+powell")))) {
    stop("nomad.nmulti is only supported when regtype='lp', automatic degree search is active, and search.engine is 'nomad' or 'nomad+powell'")
  }

  list(named = inner.named, nmulti = inner.nmulti)
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
      regtype.engine = "lp",
      basis.engine = basis,
      degree.engine = degree,
      bernstein.basis.engine = FALSE
    ))
  }

  degree <- npValidateGlpDegree(regtype = "lp",
                                degree = degree,
                                ncon = ncon)
  list(
    regtype = "lp",
    basis = basis,
    degree = degree,
    bernstein.basis = bernstein.basis,
    regtype.engine = "lp",
    basis.engine = basis,
    degree.engine = degree,
    bernstein.basis.engine = bernstein.basis
  )
}

npIsRawDegreeOneConditionalSpec <- function(spec, ncon) {
  degree <- if (is.null(spec$degree.engine)) integer(0) else as.integer(spec$degree.engine)
  ncon <- as.integer(ncon)
  identical(spec$regtype.engine, "lp") &&
    !isTRUE(spec$bernstein.basis.engine) &&
    length(degree) == ncon &&
    ncon > 0L &&
    all(degree == 1L)
}

npIsRawDegreeOneConditionalRequest <- function(regtype,
                                               degree = NULL,
                                               bernstein.basis = FALSE) {
  regtype <- as.character(regtype)[1L]
  if (identical(regtype, "ll"))
    return(TRUE)
  if (!identical(regtype, "lp"))
    return(FALSE)
  if (isTRUE(bernstein.basis) || is.null(degree))
    return(FALSE)
  degree <- suppressWarnings(as.numeric(degree))
  length(degree) > 0L &&
    all(is.finite(degree)) &&
    all(degree == 1)
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
      targeted <- !is.null(expected) &&
        identical(msg, expected) &&
        identical(bwmethod, "cv.ls") &&
        npIsRawDegreeOneConditionalSpec(spec, ncon)

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

  include <- as.integer(include)
  categories <- as.integer(categories)

  basis.code <- switch(basis, additive = 0L, glp = 1L, tensor = 2L)

  .Call("C_np_dim_basis",
        as.integer(basis.code),
        as.integer(isTRUE(kernel)),
        degree,
        segments,
        include,
        categories,
        PACKAGE = "npRmpi")
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
  .Call("C_np_set_seed", .np_validate_seed_scalar(seed), PACKAGE = "npRmpi")
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
                  PACKAGE = "npRmpi"))

}

numNotIn <- function(x){
  while(is.element(num <- rnorm(1),x)){}
  num
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


explodeFormula <- function(formula, data=NULL){
  if(any(grepl("\\.",deparse(formula)))) {
      if(is.null(data)) stop("'.' in formula and no 'data' argument")
      formula <- terms(formula, data=data)
  }
  res <- strsplit(strsplit(paste(deparse(formula), collapse=""),
                           " *[~] *")[[1]], " *[+] *")
  stopifnot(all(sapply(res,length) > 0))
  names(res) <- c("response","terms")
  res
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

genTimingStr <- function(x){
  if (is.null(x$total.time) || is.na(x$total.time))
    return("")

  nomad.time <- if (!is.null(x$nomad.time) && !is.na(x$nomad.time))
    as.double(x$nomad.time) else NA_real_
  powell.time <- if (!is.null(x$powell.time) && !is.na(x$powell.time))
    as.double(x$powell.time) else NA_real_
  fit.time <- if (!is.null(x$fit.time) && !is.na(x$fit.time))
    as.double(x$fit.time) else NA_real_
  if ((!is.finite(nomad.time) || !is.finite(powell.time)) &&
      is.list(x$bws)) {
    if (!is.finite(nomad.time) &&
        !is.null(x$bws$nomad.time) && !is.na(x$bws$nomad.time))
      nomad.time <- as.double(x$bws$nomad.time)
    if (!is.finite(powell.time) &&
        !is.null(x$bws$powell.time) && !is.na(x$bws$powell.time))
      powell.time <- as.double(x$bws$powell.time)
  }

  .npRmpiTimingProfileRecord <- function() {
    rec <- NULL

    if (is.list(x) && is.list(x$timing.profile) && !is.null(x$timing.profile$where)) {
      rec <- x$timing.profile
    } else if (is.list(x) && is.list(x$bws) && is.list(x$bws$timing.profile) &&
               !is.null(x$bws$timing.profile$where)) {
      rec <- x$bws$timing.profile
    }

    if (!is.list(rec) || is.null(rec$where))
      return(NULL)

    rec
  }

  .npRmpiTimingSessionStr <- function() {
    if (isFALSE(getOption("npRmpi.profile.summary", TRUE)))
      return("")

    if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
      return("")

    comm <- 1L
    size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) NA_integer_)
    rank <- tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) NA_integer_)
    if (is.na(size)) {
      comm <- 0L
      size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) NA_integer_)
      rank <- tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) NA_integer_)
    }
    if (is.na(size))
      return("")

    if (is.na(rank))
      rank <- NA_integer_

    nslaves <- max(as.integer(size) - 1L, 0L)
    autodispatch <- isTRUE(getOption("npRmpi.autodispatch", FALSE))

    rec <- .npRmpiTimingProfileRecord()
    ratio <- suppressWarnings(as.double(rec$comm_ratio)[1L])
    ratio.str <- if (is.finite(ratio)) {
      paste0(", overhead_ratio=", format(round(100 * ratio, 2), nsmall = 2), "%")
    } else {
      ""
    }

    paste0(
      "\nMPI Session: comm=", as.integer(comm),
      ", rank=", ifelse(is.na(rank), "NA", as.integer(rank)),
      ", size=", as.integer(size),
      ", nslaves=", as.integer(nslaves),
      ", autodispatch=", ifelse(autodispatch, "on", "off"),
      ratio.str
    )
  }

  .npRmpiTimingProfileStr <- function() {
    if (isFALSE(getOption("npRmpi.profile.summary", TRUE)))
      return("")

    rec <- .npRmpiTimingProfileRecord()

    if (!is.list(rec) || is.null(rec$where))
      return("")

    wall <- suppressWarnings(as.double(rec$wall_elapsed_sec)[1L])
    comm <- suppressWarnings(as.double(rec$comm_elapsed_sec)[1L])
    comp <- suppressWarnings(as.double(rec$compute_elapsed_sec)[1L])
    ratio <- suppressWarnings(as.double(rec$comm_ratio)[1L])
    calls <- suppressWarnings(as.integer(rec$comm_calls)[1L])
    method <- if (!is.null(rec$method)) as.character(rec$method)[1L] else NA_character_
    B <- suppressWarnings(as.integer(rec$B)[1L])
    where <- as.character(rec$where)[1L]
    kind <- if (!is.null(rec$profile_kind)) as.character(rec$profile_kind)[1L] else "bootstrap"

    fmt <- function(v, d = 4L) ifelse(is.finite(v), format(round(v, d), nsmall = d), "NA")
    ratio.pct <- if (is.finite(ratio)) paste0(format(round(100 * ratio, 2), nsmall = 2), "%") else "NA"

    if (identical(kind, "call")) {
      return(paste0(
        "\nMPI Call Profile: ", ifelse(is.na(where) || !nzchar(where), "call", where),
        "\n  wall=", fmt(wall),
        "s, comm=", fmt(comm),
        "s, compute=", fmt(comp),
        "s, comm_ratio=", ratio.pct,
        ", comm_calls=", ifelse(is.finite(calls), calls, "NA")
      ))
    }

    paste0(
      "\nMPI Bootstrap Profile: ", where,
      " [method=", ifelse(is.na(method), "NA", method),
      ", B=", ifelse(is.finite(B), B, "NA"), "]",
      "\n  wall=", fmt(wall),
      "s, comm=", fmt(comm),
      "s, compute=", fmt(comp),
      "s, comm_ratio=", ratio.pct,
      ", comm_calls=", ifelse(is.finite(calls), calls, "NA")
    )
  }

  if (is.finite(nomad.time) || is.finite(powell.time)) {
    detail <- character(0)
    if (is.finite(nomad.time))
      detail <- c(detail, paste("NOMAD ", format(nomad.time), "s", sep = ""))
    if (is.finite(powell.time))
      detail <- c(detail, paste("Powell ", format(powell.time), "s", sep = ""))
    if (is.finite(fit.time))
      detail <- c(detail, paste("fit ", format(fit.time), "s", sep = ""))

    return(paste("\nEstimation Time: ", format(x$total.time), " seconds (",
                 paste(detail, collapse = ", "), ")",
                 .npRmpiTimingSessionStr(), .npRmpiTimingProfileStr(), sep = ""))
  }

  if (!is.null(x$optim.time) && !is.na(x$optim.time) &&
      !is.null(x$fit.time) && !is.na(x$fit.time))
    return(paste("\nEstimation Time: ", format(x$total.time), " seconds (optim ",
                 format(x$optim.time), "s, fit ", format(x$fit.time), "s)",
                 .npRmpiTimingSessionStr(), .npRmpiTimingProfileStr(), sep = ""))

  paste("\nEstimation Time: ",format(x$total.time)," seconds",
        .npRmpiTimingSessionStr(), .npRmpiTimingProfileStr(), sep = "")
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
  regtype <- if (!is.null(x$regtype)) x$regtype else if (!is.null(x$bws)) x$bws$regtype else NULL
  basis <- if (!is.null(x$basis)) x$basis else if (!is.null(x$bws)) x$bws$basis else NULL
  bern <- if (!is.null(x$bernstein.basis)) x$bernstein.basis else if (!is.null(x$bws)) x$bws$bernstein.basis else NULL
  est.label <- if (identical(regtype, "lp")) npFormatRegressionType(x) else x$pregtype
  basis.family <- if (identical(regtype, "lp")) npLpBasisFamilyLabel(basis) else NULL
  basis.rep <- if (identical(regtype, "lp")) npLpBasisRepresentationLabel(bern) else NULL
  est.label.str <- if (is.null(est.label)) "" else paste("\nKernel Regression Estimator:", est.label)
  basis.family.str <- if (is.null(basis.family)) "" else paste("\nLP Basis Family:", basis.family)
  basis.rep.str <- if (is.null(basis.rep)) "" else paste("\nLP Basis Representation:", basis.rep)
  ptype.str <- if (is.null(x$ptype)) "" else paste("\nBandwidth Type:", x$ptype)
  tau.str <- if (is.null(x$tau)) "" else paste("\nTau:", x$tau)
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
  regtype <- if (!is.null(x$regtype)) {
    x$regtype
  } else if (!is.null(x$bws) && !is.null(x$bws$regtype)) {
    x$bws$regtype
  } else {
    NULL
  }

  pregtype <- if (!is.null(x$pregtype)) {
    x$pregtype
  } else if (!is.null(x$bws) && !is.null(x$bws$pregtype)) {
    x$bws$pregtype
  } else {
    NULL
  }

  if (!identical(regtype, "lp"))
    return(pregtype)

  degree <- if (!is.null(x$degree)) {
    x$degree
  } else if (!is.null(x$bws) && !is.null(x$bws$degree)) {
    x$bws$degree
  } else {
    NULL
  }

  if (is.null(degree) || length(degree) == 0)
    return("Local-Polynomial")

  basis <- if (!is.null(x$basis)) {
    x$basis
  } else if (!is.null(x$bws) && !is.null(x$bws$basis)) {
    x$bws$basis
  } else {
    "glp"
  }
  basis.family <- npLpBasisFamilyLabel(basis)

  sprintf("Local-Polynomial (%s basis; degree = %s)",
          basis.family, paste(degree, collapse = ","))
}

npBandwidthSummaryLabel <- function(bwtype, bwscaling = FALSE){
  if (isTRUE(bwscaling))
    return("Scale Factor(s)")

  if (identical(bwtype, "fixed"))
    return("Bandwidth(s)")

  "Bandwidth Nearest Neighbor(s)"
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
    if(!(is.null(x$num.feval.fast) || (length(x$num.feval.fast) == 1L && is.na(x$num.feval.fast)))){
      nfe.str <- paste(nfe.str, " (fast = ", format(x$num.feval.fast), ")", sep="")
    }
  }

  pregtype <- npFormatRegressionType(x)

  pregtype.str <- if (is.null(pregtype)) "" else paste("\nRegression Type:", pregtype)
  pmethod.str <- if (is.null(x$pmethod)) "" else paste("\nBandwidth Selection Method:", x$pmethod)
  formula.str <- if (!identical(x$formula, NULL)) paste("\nFormula:", paste(deparse(x$formula), collapse = "\n")) else ""
  ptype.str <- if (is.null(x$ptype)) "" else paste("\nBandwidth Type: ", x$ptype, sep = "")

  paste(pregtype.str,
        pmethod.str,
        formula.str,
        ptype.str,
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
      if (x$type == "fixed") "Scale Factor:" else ""
    } else {
      "Lambda Max:"
    }
  })

  maxNameLen <- max(nchar(unlist(sumText)))
  print.sumText <- lapply(sumText, '!=', "")

  sumText <- lapply(seq_along(sumText), function(i){
    paste(blank(maxNameLen - nchar(sumText[[i]])), sumText[[i]], sep="")
  })

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
      sum.str <- paste(sumText[[j]], " ", npFormat(flat_sum[[j]]), sep = "")
    paste(vatText[[j]], " Bandwidth: ", npFormat(flat_bandwidth[[j]]), " ",
          sum.str, sep = "", collapse = "")
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
            out <- out + ((-1.0)^k) * choose(r, k) * choose(m, j) * (u^j) * ((1.0 - u)^(m - j))
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
      return(matrix(1,nrow=nrow(xdat.numeric),ncol=1))
    } else {
      return(matrix(1,nrow=nrow(exdat.numeric),ncol=1))
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
