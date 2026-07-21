.np_iv_pipe_parts <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("|")))
    return(c(.np_iv_pipe_parts(expr[[2L]]), list(expr[[3L]])))
  list(expr)
}

.np_iv_expr_contains_symbol <- function(expr, symbol) {
  if (is.symbol(expr))
    return(identical(as.character(expr), symbol))
  if (!is.call(expr))
    return(FALSE)
  any(vapply(as.list(expr)[-1L], .np_iv_expr_contains_symbol,
             logical(1L), symbol = symbol))
}

.np_iv_add_exprs <- function(exprs) {
  if (length(exprs) == 1L)
    return(exprs[[1L]])
  Reduce(function(lhs, rhs) call("+", lhs, rhs), exprs)
}

.np_iv_make_formula <- function(lhs = NULL, rhs, env) {
  ans <- if (is.null(lhs)) {
    as.call(list(as.name("~"), rhs))
  } else {
    as.call(list(as.name("~"), lhs, rhs))
  }
  class(ans) <- "formula"
  environment(ans) <- env
  ans
}

.np_iv_parse_formula <- function(formula, where) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop(sprintf("%s requires a two-sided formula", where), call. = FALSE)

  env <- environment(formula)
  if (is.null(env))
    env <- parent.frame()
  parts <- .np_iv_pipe_parts(formula[[3L]])
  if (!length(parts) %in% c(2L, 3L)) {
    stop(sprintf(
      "%s formula must be 'y ~ z | w' or 'y ~ z | w | x'",
      where
    ), call. = FALSE)
  }

  role.names <- c("z", "w", "x")[seq_along(parts)]
  roles <- setNames(vector("list", length(parts)), role.names)
  for (i in seq_along(parts)) {
    expr <- parts[[i]]
    if (.np_iv_expr_contains_symbol(expr, ".")) {
      stop(sprintf(
        "%s does not support '.' in IV formula partitions; name each role explicitly",
        where
      ), call. = FALSE)
    }
    role.formula <- .np_iv_make_formula(rhs = expr, env = env)
    role.terms <- terms(role.formula)
    labels <- attr(role.terms, "term.labels")
    if (!length(labels)) {
      stop(sprintf("%s formula partition '%s' is empty", where, role.names[[i]]),
           call. = FALSE)
    }
    if (any(attr(role.terms, "order") > 1L)) {
      stop(sprintf(
        "%s does not support interaction operators in IV formula partitions; use an explicit transformed variable such as I(a * b)",
        where
      ), call. = FALSE)
    }
    roles[[i]] <- list(expr = expr, formula = role.formula,
                       terms = role.terms, labels = labels)
  }

  combined.rhs <- .np_iv_add_exprs(lapply(roles, `[[`, "expr"))
  combined <- .np_iv_make_formula(lhs = formula[[2L]], rhs = combined.rhs,
                                  env = env)
  response.vars <- all.vars(formula[[2L]])
  if (length(response.vars) != 1L)
    stop(sprintf("%s requires one response variable", where), call. = FALSE)

  list(formula = formula, combined = combined, roles = roles,
       response = response.vars[[1L]], environment = env)
}

.np_iv_training_frame <- function(parsed, method.call, data, where, eval.env) {
  if (!is.null(data)) {
    npValidateNewdataColumns(data, all.vars(parsed$combined), argname = "data")
  }

  keep <- match(c("data", "subset", "na.action"), names(method.call),
                nomatch = 0L)
  mf.call <- method.call[c(1L, keep)]
  mf.call[[1L]] <- quote(stats::model.frame)
  mf.call$formula <- parsed$combined
  if (!is.null(data))
    mf.call$data <- data
  eval(mf.call, envir = eval.env)
}

.np_iv_role_frame <- function(frame, role, where) {
  labels <- role$labels
  missing.labels <- setdiff(labels, names(frame))
  if (length(missing.labels)) {
    stop(sprintf(
      "%s could not construct IV formula terms: %s",
      where, paste(shQuote(missing.labels), collapse = ", ")
    ), call. = FALSE)
  }
  ans <- frame[, labels, drop = FALSE]
  for (j in seq_along(ans)) {
    if (is.matrix(ans[[j]])) {
      stop(sprintf(
        "%s formula term %s produces a matrix; supply its columns as explicit variables",
        where, shQuote(names(ans)[[j]])
      ), call. = FALSE)
    }
    if (inherits(ans[[j]], "AsIs")) {
      cls <- setdiff(class(ans[[j]]), "AsIs")
      if (length(cls)) class(ans[[j]]) <- cls else class(ans[[j]]) <- NULL
    }
  }
  ans
}

.np_iv_eval_frame <- function(parsed, newdata, roles, na.action, where) {
  role.parts <- parsed$roles[roles]
  eval.formula <- .np_iv_make_formula(
    rhs = .np_iv_add_exprs(lapply(role.parts, `[[`, "expr")),
    env = parsed$environment
  )
  npValidateNewdataFormula(newdata, eval.formula, include.response = FALSE,
                           argname = "newdata")
  args <- list(formula = eval.formula, data = newdata)
  if (!is.null(na.action))
    args$na.action <- na.action
  do.call(stats::model.frame, args, envir = parsed$environment)
}

.np_iv_formula_metadata <- function(fit, parsed, train.frame, eval.frame,
                                    call, trainiseval) {
  fit$call <- call
  fit$formula <- parsed$formula
  fit$terms <- terms(parsed$combined)
  fit$role.terms <- lapply(parsed$roles, `[[`, "terms")
  fit$omit <- attr(train.frame, "na.action")
  fit$rows.omit <- as.vector(fit$omit)
  fit$nobs.omit <- length(fit$rows.omit)
  fit$eval.omit <- if (is.null(eval.frame)) NULL else attr(eval.frame, "na.action")
  fit$eval.rows.omit <- as.vector(fit$eval.omit)
  fit$eval.nobs.omit <- length(fit$eval.rows.omit)
  fit$trainiseval <- isTRUE(trainiseval)
  fit
}

npregiv.formula <- function(y, data = NULL, subset, na.action,
                            newdata = NULL, zeval = NULL, xeval = NULL,
                            ...) {
  where <- "npregiv()"
  mc <- match.call()
  mc[[1L]] <- quote(npregiv)
  parsed <- .np_iv_parse_formula(y, where = where)
  train.frame <- .np_iv_training_frame(
    parsed = parsed, method.call = mc, data = data,
    where = where, eval.env = parent.frame()
  )
  ydat <- model.response(train.frame)
  zdat <- .np_iv_role_frame(train.frame, parsed$roles$z, where)
  wdat <- .np_iv_role_frame(train.frame, parsed$roles$w, where)
  xdat <- if (is.null(parsed$roles$x)) NULL else
    .np_iv_role_frame(train.frame, parsed$roles$x, where)

  eval.frame <- NULL
  formula.zeval <- NULL
  formula.xeval <- NULL
  if (!is.null(newdata)) {
    eval.roles <- c(if (is.null(zeval)) "z",
                    if (!is.null(parsed$roles$x) && is.null(xeval)) "x")
    if (length(eval.roles)) {
      eval.frame <- .np_iv_eval_frame(
        parsed = parsed, newdata = newdata, roles = eval.roles,
        na.action = if (missing(na.action)) NULL else na.action,
        where = where
      )
      if ("z" %in% eval.roles)
        formula.zeval <- .np_iv_role_frame(eval.frame, parsed$roles$z, where)
      if ("x" %in% eval.roles)
        formula.xeval <- .np_iv_role_frame(eval.frame, parsed$roles$x, where)
    }
  }

  native.args <- list(y = ydat, z = zdat, w = wdat)
  if (!is.null(xdat))
    native.args$x <- xdat
  effective.zeval <- if (!is.null(zeval)) zeval else formula.zeval
  effective.xeval <- if (!is.null(xeval)) xeval else formula.xeval
  if (!is.null(effective.zeval))
    native.args$zeval <- effective.zeval
  if (!is.null(effective.xeval))
    native.args$xeval <- effective.xeval

  fit <- do.call(npregiv, c(native.args, list(...)))
  .np_iv_formula_metadata(
    fit = fit, parsed = parsed, train.frame = train.frame,
    eval.frame = eval.frame, call = mc,
    trainiseval = is.null(effective.zeval) && is.null(effective.xeval)
  )
}

npregivderiv.formula <- function(y, data = NULL, subset, na.action,
                                 newdata = NULL, zeval = NULL, weval = NULL,
                                 xeval = NULL, ...) {
  where <- "npregivderiv()"
  mc <- match.call()
  mc[[1L]] <- quote(npregivderiv)
  parsed <- .np_iv_parse_formula(y, where = where)
  train.frame <- .np_iv_training_frame(
    parsed = parsed, method.call = mc, data = data,
    where = where, eval.env = parent.frame()
  )
  ydat <- model.response(train.frame)
  zdat <- .np_iv_role_frame(train.frame, parsed$roles$z, where)
  wdat <- .np_iv_role_frame(train.frame, parsed$roles$w, where)
  xdat <- if (is.null(parsed$roles$x)) NULL else
    .np_iv_role_frame(train.frame, parsed$roles$x, where)

  eval.frame <- NULL
  formula.zeval <- NULL
  formula.weval <- NULL
  formula.xeval <- NULL
  if (!is.null(newdata)) {
    eval.roles <- c(if (is.null(zeval)) "z", if (is.null(weval)) "w",
                    if (!is.null(parsed$roles$x) && is.null(xeval)) "x")
    if (length(eval.roles)) {
      eval.frame <- .np_iv_eval_frame(
        parsed = parsed, newdata = newdata, roles = eval.roles,
        na.action = if (missing(na.action)) NULL else na.action,
        where = where
      )
      if ("z" %in% eval.roles)
        formula.zeval <- .np_iv_role_frame(eval.frame, parsed$roles$z, where)
      if ("w" %in% eval.roles)
        formula.weval <- .np_iv_role_frame(eval.frame, parsed$roles$w, where)
      if ("x" %in% eval.roles)
        formula.xeval <- .np_iv_role_frame(eval.frame, parsed$roles$x, where)
    }
  }

  native.args <- list(y = ydat, z = zdat, w = wdat)
  if (!is.null(xdat))
    native.args$x <- xdat
  effective.zeval <- if (!is.null(zeval)) zeval else formula.zeval
  effective.weval <- if (!is.null(weval)) weval else formula.weval
  effective.xeval <- if (!is.null(xeval)) xeval else formula.xeval
  if (!is.null(effective.zeval))
    native.args$zeval <- effective.zeval
  if (!is.null(effective.weval))
    native.args$weval <- effective.weval
  if (!is.null(effective.xeval))
    native.args$xeval <- effective.xeval

  fit <- do.call(npregivderiv, c(native.args, list(...)))
  .np_iv_formula_metadata(
    fit = fit, parsed = parsed, train.frame = train.frame,
    eval.frame = eval.frame, call = mc,
    trainiseval = is.null(effective.zeval) && is.null(effective.weval) &&
      is.null(effective.xeval)
  )
}

.np_iv_validate_scalar_degree <- function(degree, argname = "degree") {
  if (!is.numeric(degree) || length(degree) != 1L || anyNA(degree) ||
      !is.finite(degree) || degree < 0 || degree != floor(degree)) {
    stop(sprintf("%s must be one finite non-negative integer", argname),
         call. = FALSE)
  }
  if (degree > .np_glp_degree_hard_max) {
    stop(sprintf("%s must lie in [0,%d]", argname, .np_glp_degree_hard_max),
         call. = FALSE)
  }
  .np_warn_high_glp_degree(degree, argname = argname)
  as.integer(degree)
}

.np_iv_p_regtype <- function(p) {
  if (p == 0L) "lc" else if (p == 1L) "ll" else "lp"
}

.np_iv_resolve_npregiv_smoothing <- function(p, p.missing, regtype, degree,
                                             nomad) {
  nomad.mode <- npValidateNomadControl(nomad, "nomad")
  if (!identical(nomad.mode, "false")) {
    stop(
      "npregiv() automatic-degree NOMAD search is not available: use nomad=FALSE with regtype='lc', regtype='ll', or regtype='lp' and a scalar degree",
      call. = FALSE
    )
  }

  legacy.p <- .np_iv_validate_scalar_degree(p, "p")
  regtype.supplied <- !is.null(regtype)
  degree.supplied <- !is.null(degree)
  if (!regtype.supplied && degree.supplied)
    stop("degree requires an explicit regtype='lp'", call. = FALSE)

  modern.p <- NULL
  requested.regtype <- NULL
  if (regtype.supplied) {
    requested.regtype <- match.arg(regtype, c("lc", "ll", "lp"))
    if (identical(requested.regtype, "lc")) {
      if (degree.supplied)
        stop("regtype='lc' does not accept degree; use regtype='lp', degree=0",
             call. = FALSE)
      modern.p <- 0L
    } else if (identical(requested.regtype, "ll")) {
      if (degree.supplied)
        stop("regtype='ll' does not accept degree; use regtype='lp', degree=1",
             call. = FALSE)
      modern.p <- 1L
    } else {
      if (!degree.supplied)
        stop("degree must be supplied explicitly when regtype='lp'",
             call. = FALSE)
      modern.p <- .np_iv_validate_scalar_degree(degree)
    }
  }

  if (!p.missing && !is.null(modern.p) && legacy.p != modern.p) {
    stop(sprintf(
      "conflicting smoothing controls: p=%d but regtype='%s'%s implies degree %d",
      legacy.p, requested.regtype,
      if (degree.supplied) sprintf(", degree=%d", modern.p) else "",
      modern.p
    ), call. = FALSE)
  }

  effective.p <- if (is.null(modern.p)) legacy.p else modern.p
  engine.p <- if (is.null(modern.p) || !p.missing) p else modern.p
  list(
    p = engine.p,
    requested = list(
      p = if (p.missing) NULL else legacy.p,
      regtype = requested.regtype,
      degree = if (degree.supplied) .np_iv_validate_scalar_degree(degree) else NULL,
      nomad = nomad.mode
    ),
    effective = list(
      regtype = .np_iv_p_regtype(effective.p),
      degree = effective.p,
      p = effective.p,
      engine = "legacy-p"
    )
  )
}

.np_iv_resolve_deriv_smoothing <- function(regtype, degree, nomad,
                                           regtype.missing,
                                           degree.missing,
                                           nomad.missing) {
  nomad.mode <- if (nomad.missing) "false" else npValidateNomadControl(nomad, "nomad")
  if (!identical(nomad.mode, "false")) {
    stop(
      "npregivderiv() NOMAD routing is not yet available for its multiple internal regression stages; use nomad=FALSE and a fixed regtype/degree",
      call. = FALSE
    )
  }

  if (regtype.missing && degree.missing) {
    return(list(
      requested = list(regtype = NULL, degree = NULL, nomad = nomad.mode),
      effective = list(regtype = "ll", degree = 1L,
                       source = "derivative-default"),
      explicit = FALSE
    ))
  }
  if (regtype.missing)
    stop("degree requires an explicit regtype='lp'", call. = FALSE)

  regtype <- match.arg(regtype, c("lc", "ll", "lp"))
  if (identical(regtype, "lc")) {
    if (!degree.missing)
      stop("regtype='lc' does not accept degree; use regtype='lp', degree=0",
           call. = FALSE)
    canonical.degree <- 0L
  } else if (identical(regtype, "ll")) {
    if (!degree.missing)
      stop("regtype='ll' does not accept degree; use regtype='lp', degree=1",
           call. = FALSE)
    canonical.degree <- 1L
  } else {
    if (degree.missing)
      stop("degree must be supplied explicitly when regtype='lp'",
           call. = FALSE)
    canonical.degree <- .np_iv_validate_scalar_degree(degree)
  }

  list(
    requested = list(regtype = regtype,
                     degree = if (degree.missing) NULL else canonical.degree,
                     nomad = nomad.mode),
    effective = list(regtype = regtype, degree = canonical.degree,
                     source = "explicit"),
    explicit = TRUE
  )
}

.np_iv_deriv_stage_args <- function(spec, txdat) {
  txdat <- toFrame(txdat)
  ncon <- sum(vapply(txdat, is.numeric, logical(1L)))
  if (ncon == 0L)
    return(list(regtype = "lc"))
  ans <- list(regtype = spec$effective$regtype)
  if (identical(spec$effective$regtype, "lp")) {
    ans$degree <- rep.int(spec$effective$degree, ncon)
  }
  ans
}

.np_iv_stage_spec <- function(name, spec, txdat) {
  txdat <- toFrame(txdat)
  ncon <- sum(vapply(txdat, is.numeric, logical(1L)))
  stage.regtype <- if (ncon == 0L) "lc" else spec$effective$regtype
  degree <- if (identical(stage.regtype, "lp")) {
    rep.int(spec$effective$degree, ncon)
  } else if (identical(stage.regtype, "ll")) {
    rep.int(1L, ncon)
  } else {
    rep.int(0L, ncon)
  }
  list(name = name, regtype = stage.regtype, degree = degree,
       ncon = ncon)
}

.np_iv_bw_payload <- function(x) {
  if (is.list(x) && !is.data.frame(x) && !is.null(x[["bw"]]))
    return(x[["bw"]])
  x
}

.np_iv_object_trainiseval <- function(object) {
  if (!is.null(object$trainiseval))
    return(isTRUE(object$trainiseval))
  if (NROW(object$y) != length(object$phi))
    return(FALSE)
  if (inherits(object, "npregiv"))
    return(is.null(object$zeval) && is.null(object$xeval))

  nx <- if (is.null(object$x)) 0L else NCOL(object$x)
  z.train <- object$z
  if (nx > 0L)
    z.train <- z.train[, seq_len(NCOL(z.train) - nx), drop = FALSE]
  same.z <- isTRUE(all.equal(as.data.frame(object$zeval),
                             as.data.frame(z.train),
                             check.attributes = FALSE))
  same.w <- isTRUE(all.equal(as.data.frame(object$weval),
                             as.data.frame(object$w),
                             check.attributes = FALSE))
  same.x <- if (is.null(object$x)) {
    is.null(object$xeval)
  } else {
    isTRUE(all.equal(as.data.frame(object$xeval),
                     as.data.frame(object$x),
                     check.attributes = FALSE))
  }
  same.z && same.w && same.x
}

.np_iv_summary_payload <- function(object, derivative = FALSE) {
  is.tikhonov <- !derivative && !is.null(object$alpha)
  nx <- if (is.null(object$x)) 0L else NCOL(object$x)
  nz <- max(0L, NCOL(object$z) - nx)
  nw <- max(0L, NCOL(object$w) - nx)
  selected <- if (derivative) object$num.iterations else object$norm.index
  stopping <- if (!is.null(selected) && length(object$norm.stop) >= selected) {
    object$norm.stop[[selected]]
  } else {
    NULL
  }
  bws <- object$bws
  if (is.null(bws)) {
    bws <- if (is.tikhonov) {
      list(E.y.w = .np_iv_bw_payload(object$bw.E.y.w),
           E.E.y.w.z = .np_iv_bw_payload(object$bw.E.E.y.w.z),
           E.phi.w = .np_iv_bw_payload(object$bw.E.phi.w),
           E.E.phi.w.z = .np_iv_bw_payload(object$bw.E.E.phi.w.z))
    } else if (derivative) {
      list(E.y.w = .np_iv_bw_payload(object$bw.E.y.w),
           E.y.z = .np_iv_bw_payload(object$bw.E.y.z))
    } else {
      list(E.y.w = object$bw.E.y.w, E.y.z = object$bw.E.y.z,
           resid.w = object$bw.resid.w,
           resid.fitted.w.z = object$bw.resid.fitted.w.z)
    }
  }
  list(
    call = object$call,
    title = if (derivative) {
      "Nonparametric Instrumental Kernel Derivative Estimation"
    } else if (is.tikhonov) {
      "Nonparametric Instrumental Kernel Regression (Tikhonov)"
    } else {
      "Nonparametric Instrumental Kernel Regression"
    },
    derivative = derivative,
    tikhonov = is.tikhonov,
    nz = nz,
    nw = nw,
    nx = nx,
    ntrain = NROW(object$y),
    p = object$p,
    smoothing.spec = object$smoothing.spec,
    alpha = object$alpha,
    alpha.iter = object$alpha.iter,
    selected = selected,
    evaluated = length(object$norm.stop),
    stopping = stopping,
    convergence = object$convergence,
    bws = bws,
    bw.names = list(w = colnames(object$w), z = colnames(object$z)),
    nmulti = object$nmulti,
    elapsed = object$ptm[[1L]],
    capabilities = list(
      fitted = !is.null(object$phi),
      gradients = if (derivative) !is.null(object$phi.prime) else !is.null(object$phi.deriv.1),
      residuals = .np_iv_object_trainiseval(object),
      predict = FALSE,
      se = FALSE
    )
  )
}

.np_iv_print_bandwidth <- function(bw, label, names = NULL) {
  bw <- .np_iv_bw_payload(bw)
  if (is.null(bw))
    return(invisible(NULL))
  if (is.matrix(bw))
    bw <- bw[nrow(bw), , drop = TRUE]
  bw <- as.numeric(bw)
  vals <- if (!is.null(names) && length(names) == length(bw)) {
    paste(paste(names, formatC(bw, digits = 8, format = "g"), sep = ": "),
          collapse = ", ")
  } else {
    paste(formatC(bw, digits = 8, format = "g"), collapse = ", ")
  }
  cat("\n", label, " ", vals, sep = "")
  invisible(NULL)
}

.np_iv_print_summary <- function(x) {
  cat("Call:\n")
  print(x$call)
  cat("\n", x$title, "\n", sep = "")
  cat("\nNumber of continuous endogenous predictors: ", format(x$nz), sep = "")
  cat("\nNumber of continuous instruments: ", format(x$nw), sep = "")
  if (x$nx > 0L)
    cat("\nNumber of continuous exogenous predictors: ", format(x$nx), sep = "")
  if (!is.null(x$p))
    cat("\nLocal polynomial order (p): ", format(x$p), sep = "")
  if (!is.null(x$smoothing.spec$effective$regtype)) {
    cat("\nRegression smoothing type: ",
        x$smoothing.spec$effective$regtype, sep = "")
    if (!is.null(x$smoothing.spec$effective$degree))
      cat("\nRegression smoothing degree: ",
          format(x$smoothing.spec$effective$degree), sep = "")
  }
  cat("\nTraining observations: ", format(x$ntrain), sep = "")

  if (isTRUE(x$tikhonov)) {
    cat("\n\nRegularization method: Tikhonov")
    cat("\nTikhonov parameter (alpha): ", format(x$alpha, digits = 8), sep = "")
    if (!is.null(x$alpha.iter))
      cat("\nIterated Tikhonov parameter (alpha.iter): ",
          format(x$alpha.iter, digits = 8), sep = "")
  } else {
    cat("\n\nRegularization method: Landweber-Fridman")
    cat("\nNumber of iterations: ", format(x$selected), sep = "")
    cat("\nStates evaluated: ", format(x$evaluated), sep = "")
    cat("\nStopping rule value: ", format(x$stopping, digits = 8), sep = "")
    if (!is.null(x$convergence))
      cat("\nConvergence status: ", x$convergence, sep = "")
  }

  if (isTRUE(x$tikhonov)) {
    .np_iv_print_bandwidth(x$bws$E.y.w, "Bandwidth for E(y|w):", x$bw.names$w)
    .np_iv_print_bandwidth(x$bws$E.E.y.w.z, "Bandwidth for E(E(y|w)|z):", x$bw.names$z)
    .np_iv_print_bandwidth(x$bws$E.phi.w, "Bandwidth for E(phi(z)|w):", x$bw.names$w)
    .np_iv_print_bandwidth(x$bws$E.E.phi.w.z, "Bandwidth for E(E(phi(z)|w)|z):", x$bw.names$z)
  } else {
    .np_iv_print_bandwidth(x$bws$E.y.w, "Bandwidth for E(y|w):", x$bw.names$w)
    .np_iv_print_bandwidth(x$bws$E.y.z, "Bandwidth for E(y|z):", x$bw.names$z)
    if (!isTRUE(x$derivative)) {
      .np_iv_print_bandwidth(x$bws$resid.w, "Bandwidth for E(y-phi(z)|w):", x$bw.names$w)
      .np_iv_print_bandwidth(x$bws$resid.fitted.w.z,
                             "Bandwidth for E(E(y-phi(z)|w)|z):", x$bw.names$z)
    }
  }
  cat("\nNumber of multistarts: ", format(x$nmulti), sep = "")
  cat("\nEstimation time: ", formatC(x$elapsed, digits = 1, format = "f"),
      " seconds\n\n", sep = "")
  invisible(x)
}

.np_iv_output_omit <- function(object) {
  if (.np_iv_object_trainiseval(object)) object$omit else object$eval.omit
}

.np_iv_restore_omit <- function(object, value) {
  omit <- .np_iv_output_omit(object)
  if (is.null(omit)) value else napredict(omit, value)
}

fitted.npregiv <- function(object, ...) {
  .np_iv_restore_omit(object, object$phi)
}

gradients.npregiv <- function(x, ...) {
  if (is.null(x$phi.deriv.1))
    stop("structural derivatives are not available in this npregiv object",
         call. = FALSE)
  .np_iv_restore_omit(x, x$phi.deriv.1)
}

residuals.npregiv <- function(object, ...) {
  if (!.np_iv_object_trainiseval(object) || NROW(object$y) != length(object$phi))
    stop("residuals are available only when the structural function was evaluated on the training rows",
         call. = FALSE)
  ans <- as.vector(object$y) - as.vector(object$phi)
  if (is.null(object$omit)) ans else naresid(object$omit, ans)
}

fitted.npregivderiv <- function(object, ...) {
  .np_iv_restore_omit(object, object$phi)
}

gradients.npregivderiv <- function(x, ...) {
  if (is.null(x$phi.prime))
    stop("structural derivatives are not available in this npregivderiv object",
         call. = FALSE)
  .np_iv_restore_omit(x, x$phi.prime)
}

residuals.npregivderiv <- function(object, ...) {
  if (!.np_iv_object_trainiseval(object) || NROW(object$y) != length(object$phi))
    stop("residuals are available only when the structural function was evaluated on the training rows",
         call. = FALSE)
  ans <- as.vector(object$y) - as.vector(object$phi)
  if (is.null(object$omit)) ans else naresid(object$omit, ans)
}

print.summary.npregiv <- function(x, ...) .np_iv_print_summary(x)

print.summary.npregivderiv <- function(x, ...) .np_iv_print_summary(x)
