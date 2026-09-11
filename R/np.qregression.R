npqreg <-
  function(bws, ...){
    mc <- match.call(expand.dots = FALSE)
    .np_validate_public_dots(mc[["..."]], "npqreg")
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npqreg",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npqreg",bws$call)
        else if (!is.call(bws))
          UseMethod("npqreg",bws)
        else
          UseMethod("npqreg",NULL)
      } else {
        UseMethod("npqreg", NULL)
      }
    } else {
      UseMethod("npqreg", NULL)
    }
  }

.npqreg.fit.control.names <- c("data", "newdata", "exdat", "tau", "gradients", "se", "tol", "small", "itmax",
                             ".np_conditional_cat_se_demand")
.npqreg.removed.solver.controls <- c("ftol",
                                     "lbc.dir", "dfc.dir", "cfac.dir", "initc.dir",
                                     "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir")

.npqreg_validate_tau <- function(tau) {
  if (!is.numeric(tau) || !length(tau) || anyNA(tau) ||
      any(!is.finite(tau)) || any(tau <= 0) || any(tau >= 1))
    stop("'tau' must contain numeric values in (0,1)")
  as.double(tau)
}

.npqreg_tau_labels <- function(tau) {
  paste0("tau=", format(tau, trim = TRUE, scientific = FALSE))
}

.npqreg_napredict_eval <- function(omit, x) {
  if (is.null(x) || !length(omit))
    return(x)
  omit <- as.integer(omit)
  keep <- seq_len(NROW(x) + length(omit))[-omit]
  if (is.null(dim(x))) {
    out <- x[rep(NA_integer_, length(x) + length(omit))]
    out[keep] <- x
    attr(out, "na.action") <- NULL
    return(out)
  }
  if (is.data.frame(x)) {
    out <- x[rep(NA_integer_, nrow(x) + length(omit)), , drop = FALSE]
    out[keep, ] <- x
    attr(out, "na.action") <- NULL
    return(out)
  }
  if (length(dim(x)) <= 2L) {
    out <- x[rep(NA_integer_, nrow(x) + length(omit)), , drop = FALSE]
    out[keep, ] <- x
    attr(out, "na.action") <- NULL
    return(out)
  }
  d <- dim(x)
  dn <- dimnames(x)
  new.dim <- c(d[1L] + length(omit), d[-1L])
  if (!is.null(dn)) {
    dn[[1L]] <- NULL
    if (length(dn) != length(new.dim) ||
        any(vapply(seq_along(dn), function(i) {
          !is.null(dn[[i]]) && length(dn[[i]]) != new.dim[[i]]
        }, logical(1L))))
      dn <- NULL
  }
  out <- array(NA_real_,
               dim = new.dim,
               dimnames = dn)
  out[keep, , ] <- x
  out
}

.npqreg_validate_newdata_terms <- function(newdata, xnames) {
  nd <- toFrame(newdata)
  missing.names <- setdiff(xnames, names(nd))
  if (length(missing.names))
    stop(sprintf(
      "newdata must contain columns: %s",
      paste(shQuote(xnames), collapse = ", ")
    ), call. = FALSE)
  invisible(TRUE)
}

.npqreg_fit_dots <- function(dots, allow.bandwidth.controls = FALSE) {
  dot.names <- names(dots)
  if (is.null(dot.names))
    dot.names <- rep("", length(dots))

  stale <- intersect(dot.names[nzchar(dot.names)], .npqreg.removed.solver.controls)
  if (length(stale) && !allow.bandwidth.controls) {
    stop(sprintf(
      "'%s' %s no longer accepted by npqreg; the canonical one-dimensional quantile extractor is controlled by 'tol', 'small', and 'itmax'",
      paste(stale, collapse = "', '"),
      if (length(stale) == 1L) "is" else "are"
    ))
  }

  if (!allow.bandwidth.controls) {
    bad <- dot.names == "" | !(dot.names %in% .npqreg.fit.control.names)
    if (any(bad))
      .np_reject_unused_dots(dots[bad], "npqreg")
  }

  keep <- (!nzchar(dot.names)) | (dot.names %in% .npqreg.fit.control.names)
  dots[keep]
}

.npqreg_reject_gradient_order_dots <- function(dots) {
  dot.names <- names(dots)
  if (is.null(dot.names))
    dot.names <- rep("", length(dots))

  bad <- dot.names %in% c("gradient.order", "gradient_order")
  if (any(bad))
    .np_reject_unused_dots(dots[which(bad)[1L]], "npqreg")

  invisible(TRUE)
}

.npqreg_validate_itmax <- function(itmax) {
  if (!is.numeric(itmax) || length(itmax) != 1L || is.na(itmax) ||
      !is.finite(itmax) || itmax < 1 || itmax != floor(itmax) ||
      itmax > .Machine$integer.max)
    stop("'itmax' must be a positive integer")
  as.integer(itmax)
}

.npqreg_quantile_clamp_attr <- "npqreg.clamp"

.npqreg_quantile_clamp <- function(quantile) {
  clamp <- attr(quantile, .npqreg_quantile_clamp_attr, exact = TRUE)
  if (is.null(clamp) || length(clamp) != length(quantile))
    return(rep.int("none", length(quantile)))
  as.character(clamp)
}

.npqreg_mark_clamped_delta <- function(delta, clamp) {
  clamp <- as.character(clamp)
  bad <- which(!is.na(clamp) & clamp != "none")
  if (!length(bad))
    return(delta)

  if (!is.null(delta$quanterr) && length(delta$quanterr) >= max(bad))
    delta$quanterr[bad] <- NA_real_
  if (is.matrix(delta$quantgrad) && nrow(delta$quantgrad) >= max(bad))
    delta$quantgrad[bad, ] <- NA_real_
  if (is.matrix(delta$quantgerr) && nrow(delta$quantgerr) >= max(bad))
    delta$quantgerr[bad, ] <- NA_real_
  delta
}

.npqreg_assert_monotone_cdf_kernel <- function(bws) {
  order <- bws$cykerorder
  if (is.null(order))
    return(invisible(TRUE))
  order <- suppressWarnings(as.integer(order[1L]))
  if (!is.finite(order) || order != 2L)
    stop("npqreg requires cykerorder = 2 for conditional-quantile inversion; higher-order dependent-variable CDF kernels can be nonmonotone",
         call. = FALSE)
  invisible(TRUE)
}

.npqreg_strip_fit_controls_from_bw_call <- function(call) {
  for (nm in c("tau", "gradients", "se", "tol", "small", "itmax", "newdata", "exdat",
               ".np_conditional_cat_se_demand")) {
    if (nm %in% names(call))
      call[[nm]] <- NULL
  }
  call
}

.npqreg_empty_rows <- function(value, n, base = FALSE) {
  key <- if (base) ".np.empty.base.rows" else ".np.empty.rows"
  flags <- attr(value, key, exact = TRUE)
  if (is.null(flags)) return(NULL)
  if (!is.integer(flags) || length(flags) != n ||
      anyNA(flags) || any(flags != 0L & flags != 1L))
    stop("npqreg received malformed empty-row status", call. = FALSE)
  if (any(flags == 1L)) flags else NULL
}

.npqreg_copy_empty_rows <- function(value, source, n) {
  for (base in c(FALSE, TRUE)) {
    key <- if (base) ".np.empty.base.rows" else ".np.empty.rows"
    attr(value, key) <- .npqreg_empty_rows(source, n, base)
  }
  value
}

.npqreg_quantile_delta_from_conditional <- function(bws,
                                                    xdat,
                                                    ydat,
                                                    exdat,
                                                    quantile,
                                                    gradients = FALSE,
                                                    tau = NULL,
                                                    tol = 1.490116e-04,
                                                    small = 1.490116e-05,
                                                    itmax = 10000L,
                                                    cdf.cache = NULL,
                                                    lp.first.se.demand = NULL,
                                                    cat.se.demand = NULL,
                                                    se = TRUE,
                                                    allow.external = FALSE) {
  se <- npValidateScalarLogical(se, "se")
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  if (length(quantile) != nrow(exdat))
    stop("quantile delta helper requires one quantile per evaluation row")
  if (ncol(ydat) != 1L)
    stop("quantile delta helper requires a single response")
  empty <- .npqreg_empty_rows(quantile, nrow(exdat), base = TRUE)
  if (!is.null(empty)) {
    if (!isTRUE(allow.external))
      stop("npqreg received empty rows during a required delta evaluation", call. = FALSE)
    omit <- which(empty == 1L)
    keep <- which(empty == 0L)
    if (any(!is.na(quantile[omit])))
      stop("npqreg received finite quantiles for empty rows", call. = FALSE)
    if (length(keep)) {
      q.keep <- as.double(quantile[keep])
      attr(q.keep, .npqreg_quantile_clamp_attr) <-
        .npqreg_quantile_clamp(quantile)[keep]
      out <- .npqreg_quantile_delta_from_conditional(
        bws, xdat, ydat, exdat[keep, , drop = FALSE], q.keep,
        gradients = gradients, tau = tau, tol = tol, small = small,
        itmax = itmax, cdf.cache = cdf.cache,
        lp.first.se.demand = lp.first.se.demand, cat.se.demand = cat.se.demand,
        se = se, allow.external = allow.external)
      flags <- .npqreg_empty_rows(out, length(keep))
      if (!is.null(flags)) {
        full.flags <- integer(nrow(exdat))
        full.flags[keep] <- flags
        attr(out, ".np.empty.rows") <- full.flags
      }
      if (se) out$quanterr <- .npqreg_napredict_eval(omit, out$quanterr)
      if (gradients) {
        out$quantgrad <- .npqreg_napredict_eval(omit, out$quantgrad)
        if (se) out$quantgerr <- .npqreg_napredict_eval(omit, out$quantgerr)
      }
      # Private diagnostics remain genuine subset fits; the public delta
      # arrays above retain the complete query layout.
      out$evaluated.rows <- keep
    } else {
      grad <- if (gradients) matrix(NA_real_, nrow(exdat), ncol(exdat),
                                   dimnames = list(NULL, names(exdat))) else NULL
      out <- list(quanterr = if (se) rep.int(NA_real_, nrow(exdat)) else NULL,
                  quantgrad = if (gradients) grad else NA,
                  quantgerr = if (se) { if (gradients) grad else NA } else NULL,
                  cdf = NULL, dens = NULL, evaluated.rows = integer(0L))
    }
    attr(out, ".np.empty.rows") <- .npreg_merge_empty_rows(
      .npqreg_empty_rows(out, nrow(exdat)),
      .npreg_merge_empty_rows(.npqreg_empty_rows(quantile, nrow(exdat)), empty))
    attr(out, ".np.empty.base.rows") <- empty
    return(out)
  }
  cat.se.demand <- .np_conditional_cat_se_demand(
    cat.se.demand, bws$xnuno + bws$xnord)
  qclamp <- .npqreg_quantile_clamp(quantile)
  if (all(!is.na(qclamp) & qclamp != "none"))
    cat.se.demand[] <- FALSE
  quantile <- as.double(quantile)

  reg.spec <- npConditionalRegEngineSpec(
    bws,
    where = "npqreg categorical effects"
  )
  glp.categorical.effects <- npGlpCategoricalEffectsRequired(
    regtype.engine = reg.spec$reg.engine,
    degree.engine = reg.spec$degree.engine,
    ncat = bws$xnuno + bws$xnord,
    gradients = gradients
  )
  if (glp.categorical.effects &&
      (is.null(tau) || length(tau) != 1L || !is.finite(tau))) {
    stop("quantile delta helper requires a finite scalar tau for categorical effects")
  }

  eydat <- stats::setNames(data.frame(quantile), names(ydat)[1L])
  cdf.obj <- .np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = TRUE,
    gradients = gradients,
    se = se,
    categorical.effects = !glp.categorical.effects,
    lp.first.se.demand = lp.first.se.demand,
    cat.se.demand = if (glp.categorical.effects) FALSE else cat.se.demand,
    allow.external = allow.external,
    .np.defer.empty.rows = TRUE
  )
  dens.obj <- .np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = FALSE,
    gradients = FALSE,
    se = FALSE,
    allow.external = allow.external,
    .np.defer.empty.rows = TRUE
  )

  flags <- .npreg_merge_empty_rows(
    .npqreg_empty_rows(cdf.obj, nrow(exdat)),
    .npqreg_empty_rows(dens.obj, nrow(exdat)))

  dens <- as.double(dens.obj$condens)
  quanterr <- NULL
  if (se) {
    quanterr <- as.double(cdf.obj$conderr) / NZD(dens)
    quanterr[!is.finite(quanterr) | quanterr < 0.0] <- NA_real_
  }

  if (!gradients) {
    out <- list(
      quanterr = quanterr,
      quantgrad = NA,
      quantgerr = if (se) NA else NULL,
      cdf = cdf.obj,
      dens = dens.obj
    )
    if (!is.null(flags)) attr(out, ".np.empty.rows") <- flags
    return(out)
  }

  dens.mat <- matrix(NZD(dens),
                     nrow = nrow(cdf.obj$congrad),
                     ncol = ncol(cdf.obj$congrad))
  grad <- -cdf.obj$congrad / dens.mat
  grad[!is.finite(grad)] <- NA_real_

  gerr <- NULL
  if (se) {
    gerr <- cdf.obj$congerr / dens.mat
    gerr[!is.finite(gerr) | gerr < 0.0] <- NA_real_
  }

  if (glp.categorical.effects) {
    cat.grad <- .npqreg_categorical_first_differences(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      tau = tau,
      tol = tol,
      small = small,
      itmax = itmax,
      cdf.cache = cdf.cache,
      allow.external = allow.external
    )
    flags <- .npreg_merge_empty_rows(flags, .npqreg_empty_rows(cat.grad, nrow(exdat)))
    cat.idx <- which(bws$ixuno | bws$ixord)
    grad[, cat.idx] <- cat.grad[, cat.idx, drop = FALSE]
    if (se)
      gerr[, cat.idx] <- NA_real_
  }

  out <- list(
    quanterr = quanterr,
    quantgrad = grad,
    quantgerr = gerr,
    cdf = cdf.obj,
    dens = dens.obj
  )
  if (!is.null(flags)) attr(out, ".np.empty.rows") <- flags
  out
}

.npqreg_selected_cdf_values <- function(bws,
                                        xdat,
                                        ydat,
                                        exdat,
                                        ycand,
                                        allow.external = FALSE) {
  ydat <- toFrame(ydat)
  yname <- names(ydat)[1L]
  eydat <- stats::setNames(data.frame(as.double(ycand)), yname)
  fit <- .np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = TRUE,
    gradients = FALSE,
    se = FALSE,
    allow.external = allow.external,
    .np.defer.empty.rows = TRUE
  )
  .npqreg_copy_empty_rows(as.double(fit$condist), fit, nrow(exdat))
}

.npqreg_selected_cdf_cache_atom <- function(x) {
  if (inherits(x, "factor"))
    return(paste0("f:", as.integer(x)))
  if (inherits(x, "Date") || inherits(x, "POSIXt"))
    return(paste0("t:", sprintf("%a", as.double(x))))
  if (is.numeric(x) || is.integer(x) || is.logical(x))
    return(paste0("n:", sprintf("%a", as.double(x))))
  paste0("c:", as.character(x))
}

.npqreg_selected_cdf_cache_row_key <- function(exdat, i) {
  row <- exdat[i, , drop = FALSE]
  parts <- vapply(row, function(x) .npqreg_selected_cdf_cache_atom(x[[1L]]), character(1L))
  paste(parts, collapse = "\r")
}

.npqreg_selected_cdf_cache_row_keys <- function(exdat) {
  exdat <- as.data.frame(exdat)
  if (!nrow(exdat))
    return(character(0L))
  vapply(seq_len(nrow(exdat)), function(i) .npqreg_selected_cdf_cache_row_key(exdat, i), character(1L))
}

.npqreg_selected_cdf_cache_key <- function(row.key, ycand) {
  paste(c(row.key, paste0("y:", sprintf("%a", as.double(ycand)))), collapse = "\r")
}

.npqreg_selected_cdf_cache_new <- function(enabled) {
  cache <- new.env(parent = emptyenv(), hash = FALSE)
  cache$enabled <- isTRUE(enabled)
  cache$store <- new.env(parent = emptyenv(), hash = TRUE)
  cache$visits <- 0L
  cache$hits <- 0L
  cache$misses <- 0L
  cache$unique <- 0L
  cache
}

.npqreg_selected_cdf_cache_clear <- function(cache) {
  if (is.environment(cache)) {
    cache$enabled <- FALSE
    cache$store <- new.env(parent = emptyenv(), hash = TRUE)
    cache$empty.store <- NULL
  }
  invisible(NULL)
}

.npqreg_selected_cdf_cache_should_enable <- function(tau, exdat) {
  isTRUE(npObjectiveCacheEnabled()) &&
    (length(tau) > 1L || anyDuplicated(as.data.frame(exdat)) > 0L)
}

.npqreg_categorical_first_differences <- function(bws,
                                                   xdat,
                                                   ydat,
                                                   exdat,
                                                   tau,
                                                   tol,
                                                   small,
                                                   itmax,
                                                   cdf.cache = NULL,
                                                   allow.external = FALSE) {
  cat.idx <- which(bws$ixuno | bws$ixord)
  out <- matrix(NA_real_, nrow = nrow(exdat), ncol = bws$xndim)
  if (!length(cat.idx))
    return(out)

  eval.quantile <- function(z) {
    row.keys <- if (is.environment(cdf.cache) && isTRUE(cdf.cache$enabled)) {
      .npqreg_selected_cdf_cache_row_keys(z)
    } else {
      NULL
    }
    .npqreg_invert_selected_cdf(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = z,
      tau = tau,
      tol = tol,
      small = small,
      itmax = itmax,
      cdf.cache = cdf.cache,
      cdf.row.keys = row.keys,
      allow.external = allow.external
    )
  }

  flags <- NULL
  for (jj in cat.idx) {
    frames <- npCategoricalFirstDifferenceFrames(
      exdat = exdat,
      index = jj,
      where = "npqreg"
    )
    upper <- eval.quantile(frames$upper)
    lower <- eval.quantile(frames$lower)
    flags <- .npreg_merge_empty_rows(flags, .npreg_merge_empty_rows(
      .npqreg_empty_rows(upper, nrow(exdat)),
      .npqreg_empty_rows(lower, nrow(exdat))))
    out[, jj] <- as.vector(upper) - as.vector(lower)
  }

  if (!is.null(flags)) attr(out, ".np.empty.rows") <- flags
  out
}

.npqreg_selected_cdf_values_cached <- function(bws,
                                               xdat,
                                               ydat,
                                               exdat,
                                               ycand,
                                               cdf.cache = NULL,
                                               row.keys = NULL,
                                               cdf.values = .npqreg_selected_cdf_values) {
  if (!is.environment(cdf.cache) || !isTRUE(cdf.cache$enabled))
    return(cdf.values(bws, xdat, ydat, exdat, ycand))

  exdat <- toFrame(exdat)
  n.eval <- nrow(exdat)
  ycand <- as.double(ycand)
  if (length(ycand) == 1L && n.eval > 1L)
    ycand <- rep.int(ycand, n.eval)
  if (length(ycand) != n.eval)
    return(cdf.values(bws, xdat, ydat, exdat, ycand))
  if (is.null(row.keys) || length(row.keys) != n.eval)
    row.keys <- .npqreg_selected_cdf_cache_row_keys(exdat)

  out <- rep.int(NA_real_, n.eval)
  keys <- character(n.eval)
  miss.first <- integer(0L)
  miss.keys <- character(0L)
  miss.map <- new.env(parent = emptyenv(), hash = TRUE)

  for (i in seq_len(n.eval)) {
    key <- .npqreg_selected_cdf_cache_key(row.keys[[i]], ycand[[i]])
    keys[[i]] <- key
    cdf.cache$visits <- cdf.cache$visits + 1L
    if (exists(key, envir = cdf.cache$store, inherits = FALSE)) {
      cdf.cache$hits <- cdf.cache$hits + 1L
      out[[i]] <- get(key, envir = cdf.cache$store, inherits = FALSE)
    } else {
      cdf.cache$misses <- cdf.cache$misses + 1L
      if (!exists(key, envir = miss.map, inherits = FALSE)) {
        assign(key, length(miss.first) + 1L, envir = miss.map)
        miss.first <- c(miss.first, i)
        miss.keys <- c(miss.keys, key)
      }
    }
  }

  if (length(miss.first)) {
    values <- cdf.values(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat[miss.first, , drop = FALSE],
      ycand = ycand[miss.first]
    )
    for (j in seq_along(miss.first)) {
      assign(miss.keys[[j]], values[[j]], envir = cdf.cache$store)
      cdf.cache$unique <- cdf.cache$unique + 1L
    }
    flags <- .npqreg_empty_rows(values, length(miss.first))
    base.flags <- .npqreg_empty_rows(values, length(miss.first), base = TRUE)
    if (!is.null(flags) || !is.null(base.flags)) {
      if (!is.environment(cdf.cache$empty.store))
        cdf.cache$empty.store <- new.env(parent = emptyenv(), hash = TRUE)
      for (j in seq_along(miss.first))
        assign(miss.keys[[j]], c(
          if (is.null(flags)) 0L else flags[[j]],
          if (is.null(base.flags)) 0L else base.flags[[j]]
        ), envir = cdf.cache$empty.store)
    }
  }

  missing.out <- which(!is.finite(out) & is.na(out))
  if (length(missing.out)) {
    for (i in missing.out)
      out[[i]] <- get(keys[[i]], envir = cdf.cache$store, inherits = FALSE)
  }
  if (is.environment(cdf.cache$empty.store)) {
    flags <- base.flags <- integer(n.eval)
    for (i in seq_len(n.eval)) {
      if (exists(keys[[i]], envir = cdf.cache$empty.store, inherits = FALSE)) {
        row.flags <- get(keys[[i]], envir = cdf.cache$empty.store, inherits = FALSE)
        flags[[i]] <- row.flags[[1L]]
        base.flags[[i]] <- row.flags[[2L]]
      }
    }
    if (any(flags == 1L)) attr(out, ".np.empty.rows") <- flags
    if (any(base.flags == 1L)) attr(out, ".np.empty.base.rows") <- base.flags
  }
  out
}

.npqreg_tau_empty_rows <- function(value, n, ntau, base = FALSE) {
  key <- if (base) ".npqreg.empty.tau.base.rows" else ".npqreg.empty.tau.rows"
  flags <- attr(value, key, exact = TRUE)
  if (is.null(flags)) return(NULL)
  if (!is.integer(flags) || !identical(dim(flags), c(as.integer(n), as.integer(ntau))) ||
      anyNA(flags) || any(flags != 0L & flags != 1L))
    stop("npqreg received malformed tau empty-row status", call. = FALSE)
  if (any(flags == 1L)) flags else NULL
}

.npqreg_copy_tau_empty_rows <- function(value, source, n, ntau) {
  value <- .npqreg_copy_empty_rows(value, source, n)
  for (base in c(FALSE, TRUE)) {
    key <- if (base) ".npqreg.empty.tau.base.rows" else ".npqreg.empty.tau.rows"
    attr(value, key) <- .npqreg_tau_empty_rows(source, n, ntau, base)
  }
  value
}

.npqreg_bind_tau_pieces <- function(pieces) {
  out <- do.call(cbind, pieces)
  n <- nrow(out)
  for (base in c(FALSE, TRUE)) {
    flags <- lapply(pieces, .npqreg_empty_rows, n = n, base = base)
    if (all(vapply(flags, is.null, logical(1L)))) next
    flags <- do.call(cbind, lapply(flags, function(x) if (is.null(x)) integer(n) else x))
    key <- if (base) ".npqreg.empty.tau.base.rows" else ".npqreg.empty.tau.rows"
    attr(out, key) <- flags
    row.key <- if (base) ".np.empty.base.rows" else ".np.empty.rows"
    attr(out, row.key) <- as.integer(rowSums(flags) > 0L)
  }
  base <- attr(out, ".npqreg.empty.tau.base.rows", exact = TRUE)
  if (!is.null(base) && any(base != base[, 1L]))
    stop("npqreg base-support status differs across tau values", call. = FALSE)
  out
}

.npqreg_collect_empty_rows <- function(out, parts, tasks, ntau = NULL) {
  n <- nrow(out)
  starts <- vapply(tasks, function(x) as.integer(x$start), integer(1L))
  sizes <- vapply(tasks, function(x) as.integer(x$bsz), integer(1L))
  expected <- cumsum(c(1L, utils::head(sizes, -1L)))
  if (anyNA(starts) || anyNA(sizes) || any(sizes < 1L) ||
      !identical(starts, as.integer(expected)) || sum(sizes) != n ||
      length(parts) != length(tasks))
    stop("npqreg received malformed chunk row identity", call. = FALSE)
  for (base in c(FALSE, TRUE)) {
    flags <- lapply(seq_along(parts), function(i) {
      if (is.null(ntau)) .npqreg_empty_rows(parts[[i]], sizes[[i]], base)
      else .npqreg_tau_empty_rows(parts[[i]], sizes[[i]], ntau, base)
    })
    if (all(vapply(flags, is.null, logical(1L)))) next
    full <- if (is.null(ntau)) integer(n) else matrix(0L, n, ntau)
    for (i in seq_along(parts)) {
      if (is.null(flags[[i]])) next
      idx <- seq.int(starts[[i]], length.out = sizes[[i]])
      if (is.null(ntau)) full[idx] <- flags[[i]] else full[idx, ] <- flags[[i]]
    }
    key <- if (base) ".np.empty.base.rows" else ".np.empty.rows"
    if (is.null(ntau)) {
      if (base && any(!is.na(out[full == 1L, 1L])))
        stop("npqreg received finite values for empty CDF rows", call. = FALSE)
      attr(out, key) <- full
    } else {
      if (base && any(full != full[, 1L]))
        stop("npqreg base-support status differs across tau values", call. = FALSE)
      if (base) {
        width <- ncol(out)/ntau
        if (!is.finite(width) || width < 1L || width != floor(width))
          stop("npqreg received malformed tau numeric layout", call. = FALSE)
        points <- out[, 1L + (seq_len(ntau) - 1L)*width, drop = FALSE]
        if (any(!is.na(points[full == 1L])))
          stop("npqreg received finite quantiles for empty rows", call. = FALSE)
      }
      tau.key <- if (base) ".npqreg.empty.tau.base.rows" else ".npqreg.empty.tau.rows"
      attr(out, tau.key) <- full
      attr(out, key) <- as.integer(rowSums(full) > 0L)
    }
  }
  out
}

.npRmpi_npqreg_parallel_context <- function(bws, comm = 1L) {
  isTRUE(isa(bws, "condbandwidth")) &&
    isTRUE(.npRmpi_has_active_slave_pool(comm = comm)) &&
    !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
    !isTRUE(.npRmpi_autodispatch_called_from_bcast())
}

.npRmpi_npqreg_chunk_size <- function(n.eval, comm = 1L) {
  n.eval <- as.integer(n.eval)
  if (is.na(n.eval) || n.eval < 1L)
    return(1L)

  opt <- suppressWarnings(as.integer(getOption("npRmpi.npqreg.chunk.size", NA_integer_))[1L])
  if (!is.na(opt) && opt > 0L)
    return(min(n.eval, opt))

  workers <- .npRmpi_bootstrap_worker_count(comm = comm)
  slots <- max(1L, workers + 1L)
  max(1L, as.integer(ceiling(n.eval / slots)))
}

.npRmpi_npqreg_parallel_min_eval <- function(comm = 1L) {
  opt <- suppressWarnings(as.integer(getOption("npRmpi.npqreg.parallel.min.eval", NA_integer_))[1L])
  if (!is.na(opt) && opt > 0L)
    return(opt)

  workers <- .npRmpi_bootstrap_worker_count(comm = comm)
  max(2000L, 8L * (workers + 1L))
}

.npRmpi_npqreg_parallel_ready <- function(bws,
                                          n.eval,
                                          comm = 1L,
                                          what = "npqreg",
                                          force.parallel = FALSE) {
  n.eval <- as.integer(n.eval)
  if (is.na(n.eval))
    return(FALSE)
  if (!isTRUE(force.parallel) &&
      n.eval < .npRmpi_npqreg_parallel_min_eval(comm = comm))
    return(FALSE)
  if (!.npRmpi_npqreg_parallel_context(bws, comm = comm))
    return(FALSE)

  .npRmpi_bootstrap_fanout_enabled(
    comm = comm,
    n = n.eval,
    B = n.eval,
    chunk.size = .npRmpi_npqreg_chunk_size(n.eval = n.eval, comm = comm),
    what = what
  )
  TRUE
}

.npRmpi_npqreg_reset_worker_comm_state <- function(comm = 1L) {
  if (isTRUE(.npRmpi_autodispatch_called_from_bcast()) ||
      !isTRUE(.npRmpi_has_active_slave_pool(comm = comm)))
    return(invisible(FALSE))

  cmd <- quote({
    try(.Call("C_np_set_local_regression_mode", FALSE, PACKAGE = "npRmpi"),
        silent = TRUE)
    try(.Call("C_np_set_active_comm", FALSE, as.integer(1L), PACKAGE = "npRmpi"),
        silent = TRUE)
    invisible(NULL)
  })
  try(.npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE),
      silent = TRUE)
  invisible(TRUE)
}

.npqreg_selected_cdf_values_parallel <- function(bws,
                                                 xdat,
                                                 ydat,
                                                 exdat,
                                                 ycand,
                                                  comm = 1L,
                                                  allow.external = FALSE) {
  exdat <- toFrame(exdat)
  n.eval <- nrow(exdat)
  if (!.npRmpi_npqreg_parallel_ready(
        bws = bws,
        n.eval = n.eval,
        comm = comm,
        what = "npqreg selected CDF"
      )) {
    return(.npRmpi_with_local_cdist_eval(.npqreg_selected_cdf_values(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      ycand = ycand,
      allow.external = allow.external
    )))
  }

  ycand <- as.double(ycand)
  on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = comm), add = TRUE)
  tasks <- .npRmpi_bootstrap_chunk_tasks(
    B = n.eval,
    chunk.size = .npRmpi_npqreg_chunk_size(n.eval = n.eval, comm = comm)
  )
  worker <- function(task, bws, xdat, ydat, exdat, ycand, allow.external) {
    idx <- seq.int(as.integer(task$start),
                   length.out = as.integer(task$bsz))
    .npRmpi_with_local_cdist_eval(.npqreg_selected_cdf_values(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat[idx, , drop = FALSE],
      ycand = ycand[idx],
      allow.external = allow.external
    ))
  }

  out <- .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker,
    ncol.out = 1L,
    what = "npqreg selected CDF",
    progress.label = "npqreg selected CDF",
    profile.where = "npqreg:selected-cdf",
    comm = comm,
    master_local_chunk = TRUE,
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    ycand = ycand,
    allow.external = allow.external,
    metadata.reducer = if (isTRUE(allow.external)) .npqreg_collect_empty_rows else NULL
  )

  .npqreg_copy_empty_rows(as.double(out[, 1L]), out, n.eval)
}

.npqreg_tau_layout <- function(grad.cols = 0L, gradients = FALSE, se = TRUE) {
  grad.cols <- if (isTRUE(gradients)) as.integer(grad.cols) else 0L
  if (length(grad.cols) != 1L || is.na(grad.cols) || grad.cols < 0L)
    stop("internal error: invalid npqreg gradient column count", call. = FALSE)
  level.se <- if (isTRUE(se)) 2L else integer()
  grad <- seq.int(2L + as.integer(se), length.out = grad.cols)
  grad.se <- if (isTRUE(se)) grad + grad.cols else integer()
  list(width = 1L + as.integer(se) + grad.cols * (1L + as.integer(se)),
       point = 1L, error = level.se, gradient = grad, gradient.error = grad.se)
}

.npqreg_quantile_delta_matrix <- function(delta, gradients = FALSE, se = TRUE) {
  pieces <- list()
  if (isTRUE(se)) pieces <- c(pieces, list(as.double(delta$quanterr)))
  if (isTRUE(gradients)) {
    pieces <- c(pieces, list(as.matrix(delta$quantgrad)))
    if (isTRUE(se)) pieces <- c(pieces, list(as.matrix(delta$quantgerr)))
  }
  if (!length(pieces))
    return(NULL)
  do.call(cbind, pieces)
}

.npqreg_tau_piece <- function(yq, delta = NULL, gradients = FALSE, se = TRUE) {
  out <- if (is.null(delta)) matrix(yq, ncol = 1L) else
    cbind(yq, .npqreg_quantile_delta_matrix(delta, gradients = gradients, se = se))
  out <- .npqreg_copy_empty_rows(out, yq, length(yq))
  flags <- .npreg_merge_empty_rows(.npqreg_empty_rows(out, length(yq)),
    .npqreg_empty_rows(delta, length(yq)))
  if (!is.null(flags)) attr(out, ".np.empty.rows") <- flags
  out
}

.npqreg_fit_tau_vector_parallel_matrix <- function(bws,
                                                   xdat,
                                                   ydat,
                                                   exdat,
                                                   tau,
                                                   gradients = FALSE,
                                                   tol,
                                                   small,
                                                   itmax,
                                                   comm = 1L,
                                                   force.parallel = FALSE,
                                                   lp.first.se.demand = NULL,
                                                   cat.se.demand = NULL,
                                                   se = TRUE,
                                                   allow.external = FALSE) {
  exdat <- toFrame(exdat)
  n.eval <- nrow(exdat)
  tau <- .npqreg_validate_tau(tau)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  se <- npValidateScalarLogical(se, "se")
  lp.first.se.demand <- .np_conditional_first_se_demand(
    lp.first.se.demand, bws$xncon)
  cat.se.demand <- .np_conditional_cat_se_demand(
    cat.se.demand, bws$xnuno + bws$xnord)
  grad.cols <- if (isTRUE(gradients)) as.integer(bws$xndim) else 0L
  if (is.na(grad.cols) || grad.cols < 0L)
    grad.cols <- 0L
  cols.per.tau <- .npqreg_tau_layout(grad.cols, gradients, se)$width

  fit_chunk <- function(ex.chunk) {
    pieces <- vector("list", length(tau))
    for (j in seq_along(tau)) {
      yq <- .npRmpi_with_local_cdist_eval(.npqreg_invert_selected_cdf(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = ex.chunk,
        tau = tau[[j]],
        tol = tol,
        small = small,
        itmax = itmax,
        parallel = FALSE,
        allow.external = allow.external
      ))
      qclamp <- .npqreg_quantile_clamp(yq)
      if (!gradients && !se) {
        pieces[[j]] <- .npqreg_tau_piece(yq)
        next
      }
      delta <- .npRmpi_with_local_cdist_eval(.npqreg_quantile_delta_from_conditional(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = ex.chunk,
        quantile = yq,
        gradients = gradients,
        tau = tau[[j]],
        tol = tol,
        small = small,
        itmax = itmax,
        lp.first.se.demand = lp.first.se.demand,
        cat.se.demand = cat.se.demand,
        se = se,
        allow.external = allow.external
      ))
      delta <- .npqreg_mark_clamped_delta(delta, qclamp)
      pieces[[j]] <- .npqreg_tau_piece(yq, delta, gradients, se)
    }
    .npqreg_bind_tau_pieces(pieces)
  }

  if (!.npRmpi_npqreg_parallel_ready(
        bws = bws,
        n.eval = n.eval,
        comm = comm,
        what = "npqreg tau block",
        force.parallel = force.parallel
      )) {
    return(fit_chunk(exdat))
  }

  on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = comm), add = TRUE)
  tasks <- .npRmpi_bootstrap_chunk_tasks(
    B = n.eval,
    chunk.size = .npRmpi_npqreg_chunk_size(n.eval = n.eval, comm = comm)
  )
  worker <- function(task, bws, xdat, ydat, exdat, tau, gradients, tol, small, itmax,
                     lp.first.se.demand, cat.se.demand, se, allow.external) {
    idx <- seq.int(as.integer(task$start),
                   length.out = as.integer(task$bsz))
    ex.chunk <- exdat[idx, , drop = FALSE]
    pieces <- vector("list", length(tau))
    for (j in seq_along(tau)) {
      yq <- .npRmpi_with_local_cdist_eval(.npqreg_invert_selected_cdf(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = ex.chunk,
        tau = tau[[j]],
        tol = tol,
        small = small,
        itmax = itmax,
        parallel = FALSE,
        allow.external = allow.external
      ))
      qclamp <- .npqreg_quantile_clamp(yq)
      if (!gradients && !se) {
        pieces[[j]] <- .npqreg_tau_piece(yq)
        next
      }
      delta <- .npRmpi_with_local_cdist_eval(.npqreg_quantile_delta_from_conditional(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = ex.chunk,
        quantile = yq,
        gradients = gradients,
        tau = tau[[j]],
        tol = tol,
        small = small,
        itmax = itmax,
        lp.first.se.demand = lp.first.se.demand,
        cat.se.demand = cat.se.demand,
        se = se,
        allow.external = allow.external
      ))
      delta <- .npqreg_mark_clamped_delta(delta, qclamp)
      pieces[[j]] <- .npqreg_tau_piece(yq, delta, gradients, se)
    }
    .npqreg_bind_tau_pieces(pieces)
  }

  .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker,
    ncol.out = length(tau) * cols.per.tau,
    what = "npqreg tau block",
    progress.label = "npqreg tau block",
    profile.where = "npqreg:tau-block",
    comm = comm,
    master_local_chunk = TRUE,
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    tau = tau,
    gradients = gradients,
    tol = tol,
    small = small,
    itmax = itmax,
    lp.first.se.demand = lp.first.se.demand,
    cat.se.demand = cat.se.demand,
    se = se,
    allow.external = allow.external,
    metadata.reducer = if (isTRUE(allow.external))
      function(out, parts, tasks)
        .npqreg_collect_empty_rows(out, parts, tasks, ntau = length(tau))
    else NULL
  )
}

.npqreg_fit_tau_vector_from_parallel_matrix <- function(mat,
                                                        tau,
                                                        gradients = FALSE,
                                                        grad.names = NULL,
                                                        se = TRUE,
                                                        expected.grad.cols = NULL) {
  tau <- .npqreg_validate_tau(tau)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  se <- npValidateScalarLogical(se, "se")
  grad.cols <- if (isTRUE(gradients)) {
    per.tau.raw <- ncol(mat) / length(tau)
    (per.tau.raw - 1L - as.integer(se)) / (1L + as.integer(se))
  } else {
    0L
  }
  if (!is.finite(grad.cols) || grad.cols < 0L || grad.cols != floor(grad.cols))
    stop("internal error: malformed npqreg parallel gradient payload", call. = FALSE)
  grad.cols <- as.integer(grad.cols)
  layout <- .npqreg_tau_layout(grad.cols, gradients, se)
  cols.per.tau <- layout$width
  if (ncol(mat) != length(tau) * cols.per.tau ||
      (gradients && !is.null(grad.names) && length(grad.names) != grad.cols) ||
      (gradients && !is.null(expected.grad.cols) &&
       !identical(grad.cols, as.integer(expected.grad.cols))))
    stop("internal error: malformed npqreg parallel tau payload", call. = FALSE)

  n.eval <- nrow(mat)
  yq <- matrix(NA_real_, nrow = n.eval, ncol = length(tau))
  yqerr <- if (se) matrix(NA_real_, nrow = n.eval, ncol = length(tau)) else NULL
  if (isTRUE(gradients)) {
    yqgrad <- array(NA_real_, dim = c(n.eval, grad.cols, length(tau)))
    yqgerr <- if (se) array(NA_real_, dim = c(n.eval, grad.cols, length(tau))) else NULL
  }

  for (j in seq_along(tau)) {
    offset <- (j - 1L) * cols.per.tau
    yq[, j] <- mat[, offset + layout$point]
    if (se) yqerr[, j] <- mat[, offset + layout$error]
    if (isTRUE(gradients) && grad.cols > 0L) {
      grad.idx <- offset + layout$gradient
      gerr.idx <- offset + layout$gradient.error
      yqgrad[, , j] <- mat[, grad.idx, drop = FALSE]
      if (se) yqgerr[, , j] <- mat[, gerr.idx, drop = FALSE]
    }
  }

  tau.labels <- .npqreg_tau_labels(tau)
  if (length(tau) == 1L) {
    return(.npqreg_copy_tau_empty_rows(list(
      yq = as.double(yq[, 1L]),
      yqerr = if (se) as.double(yqerr[, 1L]) else NULL,
      yqgrad = if (isTRUE(gradients)) {
        out <- yqgrad[, , 1L, drop = FALSE]
        dim(out) <- c(n.eval, grad.cols)
        if (!is.null(grad.names) && length(grad.names) == grad.cols)
          colnames(out) <- grad.names
        out
      } else NA,
      yqgerr = if (isTRUE(gradients) && se) {
        out <- yqgerr[, , 1L, drop = FALSE]
        dim(out) <- c(n.eval, grad.cols)
        if (!is.null(grad.names) && length(grad.names) == grad.cols)
          colnames(out) <- grad.names
        out
      } else if (se) NA else NULL
    ), mat, n.eval, length(tau)))
  }

  colnames(yq) <- tau.labels
  if (se) colnames(yqerr) <- tau.labels
  if (isTRUE(gradients)) {
    dimnames(yqgrad) <- list(NULL, NULL, tau.labels)
    if (se) dimnames(yqgerr) <- list(NULL, NULL, tau.labels)
    if (!is.null(grad.names) && length(grad.names) == grad.cols) {
      dimnames(yqgrad)[[2L]] <- grad.names
      if (se) dimnames(yqgerr)[[2L]] <- grad.names
    }
  }
  .npqreg_copy_tau_empty_rows(list(
    yq = yq,
    yqerr = yqerr,
    yqgrad = if (isTRUE(gradients)) yqgrad else NA,
    yqgerr = if (isTRUE(gradients) && se) yqgerr else if (se) NA else NULL
  ), mat, n.eval, length(tau))
}

.npqreg_assert_selected_cdf_metadata <- function(bws) {
  .npqreg_assert_monotone_cdf_kernel(bws)
  npConditionalRegEngineSpec(
    bws,
    where = "selected conditional distribution bandwidth"
  )
  invisible(TRUE)
}

.npqreg_invert_selected_cdf <- function(bws,
                                        xdat,
                                        ydat,
                                        exdat,
                                        tau,
                                        tol,
                                        small,
                                        itmax,
                                        parallel = FALSE,
                                        comm = 1L,
                                        cdf.cache = NULL,
                                        cdf.row.keys = NULL,
                                        cdf.values = NULL,
                                        allow.external = FALSE) {
  .npqreg_assert_selected_cdf_metadata(bws)

  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  itmax <- .npqreg_validate_itmax(itmax)
  y <- as.double(ydat[[1L]])
  y <- y[is.finite(y)]
  if (!length(y))
    stop("npqreg selected-CDF inversion requires finite response support")

  n.eval <- nrow(exdat)
  y.min <- min(y)
  y.max <- max(y)
  if (!is.finite(y.min) || !is.finite(y.max))
    stop("npqreg selected-CDF inversion found non-finite response support")
  if (identical(y.min, y.max)) {
    out <- rep.int(y.min, n.eval)
    attr(out, .npqreg_quantile_clamp_attr) <- rep.int("constant", n.eval)
    return(out)
  }

  lo <- rep.int(y.min, n.eval)
  hi <- rep.int(y.max, n.eval)
  cdf_values <- if (isTRUE(parallel)) {
    function(bws, xdat, ydat, exdat, ycand) {
      .npqreg_selected_cdf_values_parallel(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = exdat,
        ycand = ycand,
        comm = comm,
        allow.external = allow.external
      )
    }
  } else {
    function(bws, xdat, ydat, exdat, ycand) {
      .npRmpi_with_local_cdist_eval(.npqreg_selected_cdf_values(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = exdat,
        ycand = ycand,
        allow.external = allow.external
      ))
    }
  }
  if (is.function(cdf.values))
    cdf_values <- cdf.values
  cdf_values_cached <- function(bws, xdat, ydat, exdat, ycand, row.keys = NULL) {
    .npqreg_selected_cdf_values_cached(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      ycand = ycand,
      cdf.cache = cdf.cache,
      row.keys = row.keys,
      cdf.values = cdf_values
    )
  }

  flo <- cdf_values_cached(bws, xdat, ydat, exdat, lo, row.keys = cdf.row.keys)
  fhi <- cdf_values_cached(bws, xdat, ydat, exdat, hi, row.keys = cdf.row.keys)
  empty <- .npqreg_empty_rows(flo, n.eval, base = TRUE)
  empty.hi <- .npqreg_empty_rows(fhi, n.eval, base = TRUE)
  if (!identical(empty, empty.hi))
    stop("npqreg received inconsistent base-support status across CDF brackets", call. = FALSE)
  if (!is.null(empty) && !isTRUE(allow.external))
    stop("npqreg received empty rows during a required CDF evaluation", call. = FALSE)
  empty.idx <- if (is.null(empty)) integer(0L) else which(empty == 1L)
  if (length(empty.idx) &&
      (any(!is.na(flo[empty.idx])) || any(!is.na(fhi[empty.idx]))))
    stop("npqreg received finite values for empty CDF rows", call. = FALSE)
  flags <- .npreg_merge_empty_rows(
    .npqreg_empty_rows(flo, n.eval), .npqreg_empty_rows(fhi, n.eval))
  flags <- .npreg_merge_empty_rows(flags, empty)
  invalid.bracket <- if (is.null(empty)) {
    any(!is.finite(flo)) || any(!is.finite(fhi))
  } else {
    any(!is.finite(flo[empty == 0L])) || any(!is.finite(fhi[empty == 0L]))
  }
  if (invalid.bracket)
    stop("npqreg selected-CDF inversion encountered non-finite bracket values")

  done.low <- flo >= tau
  done.high <- fhi < tau
  done.low[empty.idx] <- done.high[empty.idx] <- FALSE
  active <- !(done.low | done.high)
  active[empty.idx] <- FALSE

  maxiter <- itmax
  iter <- 0L
  while (any(active) && iter < maxiter) {
    iter <- iter + 1L
    mid <- (lo[active] + hi[active]) / 2.0
    fmid <- cdf_values_cached(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat[active, , drop = FALSE],
      ycand = mid,
      row.keys = if (is.null(cdf.row.keys)) NULL else cdf.row.keys[active]
    )
    if (!is.null(.npqreg_empty_rows(fmid, sum(active), base = TRUE)))
      stop("npqreg base-support status changed during CDF refinement", call. = FALSE)
    mid.flags <- .npqreg_empty_rows(fmid, sum(active))
    if (!is.null(mid.flags)) {
      expanded <- integer(n.eval)
      expanded[active] <- mid.flags
      flags <- .npreg_merge_empty_rows(flags, expanded)
    }
    if (any(!is.finite(fmid)))
      stop("npqreg selected-CDF inversion encountered non-finite refinement values")

    active.idx <- which(active)
    upper <- fmid >= tau
    hi[active.idx[upper]] <- mid[upper]
    lo[active.idx[!upper]] <- mid[!upper]

    width <- hi[active.idx] - lo[active.idx]
    scale <- pmax(abs(hi[active.idx]), abs(lo[active.idx]), 1.0)
    active[active.idx] <- width > (tol * scale + small)
  }

  if (any(active))
    stop("npqreg selected-CDF inversion failed to converge within 'itmax'")

  out <- (lo + hi) / 2.0
  out[done.low] <- y.min
  out[done.high] <- y.max
  out[empty.idx] <- NA_real_
  clamp <- rep.int("none", n.eval)
  clamp[done.low] <- "lower"
  clamp[done.high] <- "upper"
  attr(out, .npqreg_quantile_clamp_attr) <- clamp
  if (!is.null(flags)) attr(out, ".np.empty.rows") <- flags
  if (!is.null(empty)) attr(out, ".np.empty.base.rows") <- empty
  out
}

npqreg.formula <-
  function(bws, data = NULL, newdata = NULL, ..., se = FALSE){
    se <- npValidateScalarLogical(se, "se")

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    if (!is.null(data))
      tmf[["data"]] <- data
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    tydat <- tmf[, bws$variableNames[["response"]], drop = FALSE]
    txdat <- tmf[, bws$variableNames[["terms"]], drop = FALSE]

    has.eval <- !is.null(newdata)
    if (has.eval) {
      .npqreg_validate_newdata_terms(newdata, bws$variableNames[["terms"]])
      tt <- drop.terms(tt, match(bws$variableNames$response, attr(tt, 'term.labels')))
      umf.args <- list(formula = tt, data = newdata)
      umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
      emf <- umf
      exdat <- emf[, bws$variableNames[["terms"]], drop = FALSE]
    }

    q.args <- list(txdat = txdat, tydat = tydat, se = se)
    if (has.eval)
      q.args$exdat <- exdat
    q.args$bws <- bws
    tbw <- do.call(npqreg, c(q.args, .npqreg_fit_dots(list(...))))

    tbw$omit <- attr(umf,"na.action")
    tbw$rows.omit <- as.vector(tbw$omit)
    tbw$nobs.omit <- length(tbw$rows.omit)

    tbw$quantile <- .npqreg_napredict_eval(tbw$omit, tbw$quantile)
    tbw$quanterr <- .npqreg_napredict_eval(tbw$omit, tbw$quanterr)

    if(tbw$gradients){
        tbw$quantgrad <- .npqreg_napredict_eval(tbw$omit, tbw$quantgrad)
        tbw$quantgerr <- .npqreg_napredict_eval(tbw$omit, tbw$quantgerr)
    }

    return(tbw)
  }

npqreg.call <-
  function(bws, ...) {
    npqreg(txdat = .np_eval_bws_call_arg(bws, "xdat"),
           tydat = .np_eval_bws_call_arg(bws, "ydat"),
           bws = bws, ...)
  }

npqreg.conbandwidth <-
  function(bws, ...){
    stop("incorrect bandwidth type: expected conditional distribution bandwidths instead of conditional density bandwidths")
  }

.npRmpi_npqreg_should_localize <- function(bws) {
  isa(bws, "condbandwidth")
}

.npRmpi_npqreg_eval_local_no_dispatch <- function(expr) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  options(npRmpi.autodispatch.context = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
  on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = 1L), add = TRUE)
  force(expr)
}

npqreg.condbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat,
           tau = 0.5,
           gradients = FALSE,
           tol = 1.490116e-04,
           small = 1.490116e-05, itmax = 10000,
           ..., se = FALSE){

    fit.start <- proc.time()[3]
    se <- npValidateScalarLogical(se, "se")
    tau <- .npqreg_validate_tau(tau)
    fit.dots <- list(...)
    cat.se.demand <- .np_conditional_cat_se_demand(
      fit.dots[[".np_conditional_cat_se_demand", exact = TRUE]],
      bws$xnuno + bws$xnord)
    fit.dots[[".np_conditional_cat_se_demand"]] <- NULL
    fit.dots <- .npqreg_fit_dots(fit.dots)
    if (length(fit.dots))
      stop(sprintf("unused npqreg fit argument '%s'", names(fit.dots)[1L]))
    gradients <- npValidateScalarLogical(gradients, "gradients")
    itmax <- .npqreg_validate_itmax(itmax)
    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol <= 0)
      stop("'tol' must be a positive finite numeric scalar")
    if (!is.numeric(small) || length(small) != 1L || is.na(small) ||
        !is.finite(small) || small <= 0)
      stop("'small' must be a positive finite numeric scalar")
    tol <- as.double(tol)
    small <- as.double(small)
    .npRmpi_require_active_slave_pool(where = "npqreg()")
    parallel.cond <- .npRmpi_npqreg_parallel_context(bws, comm = 1L)
    if (isTRUE(parallel.cond)) {
      n.eval.pre <- suppressWarnings(as.integer(tryCatch(
        if (missing(exdat)) NROW(txdat) else NROW(exdat),
        error = function(e) NA_integer_
      ))[1L])
      if (is.na(n.eval.pre) ||
          n.eval.pre < .npRmpi_npqreg_parallel_min_eval(comm = 1L))
        parallel.cond <- FALSE
    }
    if (isTRUE(parallel.cond))
      on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = 1L), add = TRUE)
    dispatch.call <- match.call()
    dispatch.call$se <- se
    dispatch.call$.np_conditional_cat_se_demand <- cat.se.demand
    if (.npRmpi_npqreg_should_localize(bws) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
        !isTRUE(.npRmpi_autodispatch_in_context()) &&
        !isTRUE(parallel.cond))
      return(.npRmpi_npqreg_eval_local_no_dispatch(
        .npRmpi_eval_without_dispatch(dispatch.call, parent.frame())
      ))
    if (.npRmpi_autodispatch_active() && !isTRUE(parallel.cond))
      return(.npRmpi_autodispatch_call(dispatch.call, parent.frame()))

    no.ex = missing(exdat)

    txdat = toFrame(txdat)
    tydat = toFrame(tydat)

    ntau <- length(tau)
    tau.labels <- .npqreg_tau_labels(tau)

    if (dim(tydat)[2] != 1)
      stop("'tydat' has more than one column")

    if (!no.ex){
      exdat = toFrame(exdat)
      
      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")
    }

    if (length(bws$xbw) != length(txdat))
      stop("length of bandwidth vector does not match number of columns of 'txdat'")

    if (length(bws$ybw) != 1)
      stop("length of bandwidth vector does not match number of columns of 'tydat'")

    if (any(bws$iyord) || any(bws$iyuno) || coarseclass(tydat[,1]) != "numeric")
      stop("'tydat' is not continuous")

    if ((any(bws$ixcon) &&
         !all(vapply(txdat[, bws$ixcon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$ixord) &&
         !all(vapply(txdat[, bws$ixord, drop = FALSE], inherits, logical(1), "ordered"))) ||
        (any(bws$ixuno) &&
         !all(vapply(txdat[, bws$ixuno, drop = FALSE], inherits, logical(1), "factor"))))
      stop("supplied bandwidths do not match 'txdat' in type")

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(data.frame(txdat, tydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]
    tydat <- tydat[keep.rows,,drop = FALSE]

    eval.omit <- NULL
    if (!no.ex){
      keep.eval <- rep_len(TRUE, nrow(exdat))
      eval.omit <- attr(na.omit(exdat), "na.action")
      if (length(eval.omit) > 0L)
        keep.eval[as.integer(eval.omit)] <- FALSE
      exdat <- exdat[keep.eval,,drop = FALSE]
    }
    
    tnrow = dim(txdat)[1]
    enrow = (if (no.ex) tnrow else dim(exdat)[1])

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    txdat <- adjustLevels(txdat, bws$xdati)
    tydat <- adjustLevels(tydat, bws$ydati)
    
    if (!no.ex){
      exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    }

    ## grab the evaluation data before it is converted to numeric
    if(no.ex){
      txeval <- txdat
    } else {
      txeval <- exdat
    }
    txdat.df <- txdat
    tydat.df <- tydat
    empty.flags <- NULL
    empty.row.labels <- rownames(txeval)
    if (!no.ex)
      exdat.df <- exdat

    if (isTRUE(parallel.cond)) {
      mat <- .npqreg_fit_tau_vector_parallel_matrix(
        bws = bws,
        xdat = txdat.df,
        ydat = tydat.df,
        exdat = txeval,
        tau = tau,
        gradients = gradients,
        tol = tol,
        small = small,
        itmax = itmax,
        comm = 1L,
        force.parallel = TRUE,
        cat.se.demand = cat.se.demand,
        se = se,
        allow.external = !no.ex
      )
      empty.flags <- .npqreg_empty_rows(mat, nrow(txeval))
      myout <- .npqreg_fit_tau_vector_from_parallel_matrix(
        mat,
        tau = tau,
        gradients = gradients,
        se = se,
        expected.grad.cols = bws$xndim
      )
    } else {
      cdf.cache <- .npqreg_selected_cdf_cache_new(
        .npqreg_selected_cdf_cache_should_enable(tau = tau, exdat = txeval)
      )
      cdf.row.keys <- if (isTRUE(cdf.cache$enabled))
        .npqreg_selected_cdf_cache_row_keys(txeval)
      else NULL
      on.exit(.npqreg_selected_cdf_cache_clear(cdf.cache), add = TRUE)

      fit_one_tau <- function(tau_i) {
        yq <- .npqreg_invert_selected_cdf(
          bws = bws,
          xdat = txdat.df,
          ydat = tydat.df,
          exdat = txeval,
          tau = tau_i,
          tol = tol,
          small = small,
          itmax = itmax,
          cdf.cache = cdf.cache,
          cdf.row.keys = cdf.row.keys,
          allow.external = !no.ex
        )
        empty.flags <<- .npreg_merge_empty_rows(empty.flags,
          .npqreg_empty_rows(yq, nrow(txeval)))
        qclamp <- .npqreg_quantile_clamp(yq)
        if (!gradients && !se)
          return(list(yq = yq, yqerr = NULL, yqgrad = NA, yqgerr = NULL))
        qdelta <- .npRmpi_with_local_cdist_eval(
          .npqreg_quantile_delta_from_conditional(
            bws = bws,
            xdat = txdat.df,
            ydat = tydat.df,
            exdat = txeval,
            quantile = yq,
            gradients = gradients,
            tau = tau_i,
            tol = tol,
            small = small,
            itmax = itmax,
            cdf.cache = cdf.cache,
            cat.se.demand = cat.se.demand,
            se = se,
            allow.external = !no.ex
          )
        )
        empty.flags <<- .npreg_merge_empty_rows(empty.flags,
          .npqreg_empty_rows(qdelta, nrow(txeval)))
        qdelta <- .npqreg_mark_clamped_delta(qdelta, qclamp)
        list(
          yq = yq,
          yqerr = qdelta$quanterr,
          yqgrad = if (gradients) qdelta$quantgrad else NA,
          yqgerr = if (se) { if (gradients) qdelta$quantgerr else NA } else NULL
        )
      }

      tau.out <- lapply(tau, fit_one_tau)

      if (ntau == 1L) {
        myout <- tau.out[[1L]]
      } else {
        myout <- list(
          yq = do.call(cbind, lapply(tau.out, `[[`, "yq")),
          yqerr = if (se) do.call(cbind, lapply(tau.out, `[[`, "yqerr")) else NULL,
          yqgrad = NA,
          yqgerr = if (se) NA else NULL
        )
        colnames(myout$yq) <- tau.labels
        if (se)
          colnames(myout$yqerr) <- tau.labels
        if (gradients) {
          p <- ncol(tau.out[[1L]]$yqgrad)
          grad.names <- colnames(tau.out[[1L]]$yqgrad)
          myout$yqgrad <- array(NA_real_,
                                dim = c(enrow, p, ntau),
                                dimnames = list(NULL, grad.names, tau.labels))
          if (se) myout$yqgerr <- array(NA_real_,
                                dim = c(enrow, p, ntau),
                                dimnames = list(NULL, grad.names, tau.labels))
          for (j in seq_len(ntau)) {
            myout$yqgrad[, , j] <- tau.out[[j]]$yqgrad
            if (se) myout$yqgerr[, , j] <- tau.out[[j]]$yqgerr
          }
        }
      }
    }

    if (!no.ex && length(eval.omit) > 0L) {
      myout$yq <- .npqreg_napredict_eval(eval.omit, myout$yq)
      myout$yqerr <- .npqreg_napredict_eval(eval.omit, myout$yqerr)
      if (gradients) {
        myout$yqgrad <- .npqreg_napredict_eval(eval.omit, myout$yqgrad)
        myout$yqgerr <- .npqreg_napredict_eval(eval.omit, myout$yqgerr)
      }
      txeval <- .npqreg_napredict_eval(eval.omit, txeval)
    }


    fit.elapsed <- proc.time()[3] - fit.start
    optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
    total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)

    fit <- qregression(bws = bws,
                xeval = txeval,
                tau = tau,
                quantile = myout$yq,
                quanterr = myout$yqerr,
                quantgrad = myout$yqgrad,
                quantgerr = myout$yqgerr,
                ntrain = tnrow,
                trainiseval = no.ex,
                gradients = gradients,
                se = se,
                timing = bws$timing, total.time = total.time,
                optim.time = optim.time, fit.time = fit.elapsed)
    .npreg_finish_empty_rows(fit, empty.flags, omitted = eval.omit,
                             row.labels = empty.row.labels, owner = "npqreg")
  }


npqreg.default <- function(bws, txdat, tydat, nomad = FALSE, ..., se = FALSE){
  se <- npValidateScalarLogical(se, "se")
  nomad <- npValidateNomadControl(nomad, "nomad")
  early.dots <- list(...)
  .npqreg_reject_gradient_order_dots(early.dots)
  if ("tau" %in% names(early.dots))
    .npqreg_validate_tau(early.dots$tau)

  if (!missing(bws) && inherits(bws, "formula")) {
    dots <- list(...)
    dot.names <- names(dots)
    if (is.null(dot.names))
      dot.names <- rep("", length(dots))
    fit.names <- c("newdata", "exdat", "tau", "gradients", "tol", "small", "itmax",
                   ".np_conditional_cat_se_demand")
    fit.dots <- .npqreg_fit_dots(dots[nzchar(dot.names) & dot.names %in% fit.names])
    bw.dots <- dots[!(nzchar(dot.names) & dot.names %in% fit.names)]
    bw.args <- c(list(formula = bws, nomad = nomad), bw.dots)
    bw.call <- as.call(c(list(quote(npcdistbw)), bw.args))
    use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
      bw.call,
      caller_env = parent.frame()
    )
    tbw <- if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting conditional distribution bandwidth",
        do.call(npcdistbw, bw.args)
      )
    } else {
      do.call(npcdistbw, bw.args)
    }
    return(do.call(npqreg, c(list(bws = tbw, se = se), fit.dots)))
  }

  if (!missing(txdat) && inherits(txdat, "formula") &&
      !missing(bws) && !isa(bws, "condbandwidth")) {
    dots <- list(...)
    dot.names <- names(dots)
    if (is.null(dot.names))
      dot.names <- rep("", length(dots))
    fit.names <- c("newdata", "exdat", "tau", "gradients", "tol", "small", "itmax",
                   ".np_conditional_cat_se_demand")
    fit.dots <- .npqreg_fit_dots(dots[nzchar(dot.names) & dot.names %in% fit.names])
    bw.dots <- dots[!(nzchar(dot.names) & dot.names %in% fit.names)]
    bw.args <- c(list(formula = txdat, bws = bws, nomad = nomad), bw.dots)
    bw.call <- as.call(c(list(quote(npcdistbw)), bw.args))
    use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
      bw.call,
      caller_env = parent.frame()
    )
    tbw <- if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting conditional distribution bandwidth",
        do.call(npcdistbw, bw.args)
      )
    } else {
      do.call(npcdistbw, bw.args)
    }
    return(do.call(npqreg, c(list(bws = tbw, se = se), fit.dots)))
  }

  .npRmpi_require_active_slave_pool(where = "npqreg()")
  parallel.cond <- (!missing(bws)) &&
    .npRmpi_npqreg_should_localize(bws) &&
    .npRmpi_npqreg_parallel_context(bws, comm = 1L)
  if (isTRUE(parallel.cond)) {
    n.eval.pre <- suppressWarnings(as.integer(tryCatch(
      if (missing(txdat)) NA_integer_ else NROW(txdat),
      error = function(e) NA_integer_
    ))[1L])
    if (is.na(n.eval.pre) ||
        n.eval.pre < .npRmpi_npqreg_parallel_min_eval(comm = 1L))
      parallel.cond <- FALSE
  }
  if (isTRUE(parallel.cond))
    on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = 1L), add = TRUE)
  dispatch.call <- match.call()
  dispatch.call$se <- se
  dispatch.call$.np_conditional_cat_se_demand <-
    early.dots[[".np_conditional_cat_se_demand", exact = TRUE]]
  if (!missing(bws) &&
      .npRmpi_npqreg_should_localize(bws) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
      !isTRUE(.npRmpi_autodispatch_in_context()) &&
      !isTRUE(parallel.cond))
    return(.npRmpi_npqreg_eval_local_no_dispatch(
      .npRmpi_eval_without_dispatch(dispatch.call, parent.frame())
    ))
  if (.npRmpi_autodispatch_active() && !isTRUE(parallel.cond))
    return(.npRmpi_autodispatch_call(dispatch.call, parent.frame()))

  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  txdat.named <- any(sc.names == "txdat")
  tydat.named <- any(sc.names == "tydat")

  no.bws <- missing(bws)
  no.txdat <- missing(txdat)
  no.tydat <- missing(tydat)
  has.explicit.bws <- (!no.bws) && isa(bws, "condbandwidth")
  bws.formula <- (!no.bws) && inherits(bws, "formula")

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  sc.bw <- sc
  if (bws.formula) {
    sc.bw$`bws` <- NULL
    bws.named <- FALSE
  }

  sc.bw[[1]] <- quote(npcdistbw)
  sc.bw <- .npqreg_strip_fit_controls_from_bw_call(sc.bw)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('txdat','tydat')
  nstxy <- c('xdat','ydat')
  
  m.txy <- match(ostxy, names(sc.bw), nomatch = 0)

  if(any(m.txy > 0)) {
    names(sc.bw)[m.txy] <- nstxy[m.txy > 0]
  }
  sc.bw <- .np_public_dots_filter_call(sc.bw, "npcdistbw")

  use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
    sc.bw,
    caller_env = parent.frame()
  )

  tbw <- if (!has.explicit.bws) {
    if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting conditional distribution bandwidth",
        .np_eval_bw_call(sc.bw, caller_env = parent.frame())
      )
    } else {
      .np_eval_bw_call(sc.bw, caller_env = parent.frame())
    }
  } else {
    .np_eval_bw_call(sc.bw, caller_env = parent.frame())
  }

  call.args <- list(bws = tbw, se = se)
  if (no.bws) {
    call.args$txdat <- txdat
    call.args$tydat <- tydat
  } else {
    if (txdat.named) call.args$txdat <- txdat
    if (tydat.named) call.args$tydat <- tydat
    if ((!bws.named) && (!txdat.named) && (!no.tydat) && (!tydat.named)) {
      call.args <- c(call.args, list(tydat))
    }
  }
  dots <- list(...)
  if (has.explicit.bws)
    fit.dots <- .npqreg_fit_dots(dots)
  else
    fit.dots <- .npqreg_fit_dots(dots, allow.bandwidth.controls = TRUE)
  do.call(npqreg, c(call.args, fit.dots))
}
