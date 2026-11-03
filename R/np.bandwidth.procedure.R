# Value-only provenance for density/CDF outer bandwidth selection. This records
# the successful selector's rule; it never evaluates a stored call or changes
# the selection already performed.
.np_bw_procedure_signature <- function(bws) {
  fields <- c('bw', 'method', 'pmethod', 'scaling', 'type', 'ckertype',
    'ckerorder', 'ckerbound', 'ckerlb', 'ckerub', 'ukertype', 'okertype',
    'bandwidth', 'sfactor', 'xnames', 'icon', 'iuno', 'iord')
  setNames(lapply(fields, function(k) bws[[k, exact = TRUE]]), fields)
}

.np_bw_procedure_kernel <- function(bws) {
  fields <- c(bwmethod = 'method', bwscaling = 'scaling', bwtype = 'type',
    ckertype = 'ckertype', ckerorder = 'ckerorder', ckerbound = 'ckerbound',
    ckerlb = 'ckerlb', ckerub = 'ckerub', ukertype = 'ukertype', okertype = 'okertype')
  setNames(lapply(fields, function(k) bws[[k, exact = TRUE]]), names(fields))
}

.np_bw_procedure_values <- function(x) {
  if (is.null(x)) return(TRUE)
  if (is.environment(x) || is.function(x) || is.language(x) ||
      is.expression(x) || typeof(x) == 'externalptr') return(FALSE)
  if (!(is.atomic(x) || is.list(x))) return(FALSE)
  if (is.list(x) && !all(vapply(x, .np_bw_procedure_values, logical(1)))) return(FALSE)
  a <- attributes(x)
  is.null(a) || all(vapply(a, .np_bw_procedure_values, logical(1)))
}

.np_bw_procedure_record <- function(value, original, family, frame) {
  compute <- get('bandwidth.compute', frame, inherits = FALSE)
  evaluation <- get('eval.only', frame, inherits = FALSE)
  rule <- if (!compute) 'manual' else if (identical(original[['method', exact = TRUE]],
    'normal-reference')) 'normal-reference' else 'cv'
  args <- .np_bw_procedure_kernel(original)
  args$bws <- original[['bw', exact = TRUE]]
  args$bandwidth.compute <- compute
  if (rule == 'cv' && !evaluation) {
    fields <- c('nmulti', 'powell.remin', 'bwsolver', 'itmax', 'ftol', 'tol',
      'small', 'invalid.penalty', 'penalty.multiplier',
      'scale.factor.init.lower', 'scale.factor.init.upper', 'scale.factor.init',
      'lbd.init', 'hbd.init', 'dfac.init', 'scale.init.categorical.sample',
      'transform.bounds', 'scale.factor.search.lower')
    solver <- get('bwsolver', frame, inherits = FALSE)
    if (!npBwsolverUsesMads(solver)) fields <- c(fields, 'cfac.dir', 'dfac.dir',
      'dfc.dir', 'hbd.dir', 'initc.dir', 'initd.dir', 'lbc.dir', 'lbd.dir')
    args <- c(args, mget(fields, frame, inherits = FALSE))
    if (npBwsolverUsesMads(solver)) {
      dots <- get('dot.args', frame, inherits = FALSE)
      keep <- intersect(c('mads.nmulti', 'nomad.nmulti', 'nomad.remin', 'nomad.opts'), names(dots))
      args <- c(args, dots[keep])
    }
    if (family == 'distribution') args <- c(args,
      mget(c('gdat', 'do.full.integral', 'ngrid', 'memfac'), frame, inherits = FALSE))
  }
  complete <- (!evaluation || !compute) && .np_bw_procedure_values(args)
  value[['selection.procedure']] <- list(version = 1L, family = family,
    rule = rule, complete = complete,
    reason = if (evaluation && compute) 'evaluation-only object' else if (!complete)
      'selection settings are not value-only data' else NULL,
    arguments = if (complete) args else NULL,
    signature = .np_bw_procedure_signature(value))
  value
}

.np_bw_outer_selection <- function(bws, family) {
  if (!family %in% c('density', 'distribution')) stop('invalid outer bandwidth family', call. = FALSE)
  descriptor <- bws[['selection.procedure', exact = TRUE]]
  messages <- character()
  if (!is.null(descriptor)) {
    if (!is.list(descriptor) || !identical(descriptor[['version', exact = TRUE]], 1L) ||
        !identical(descriptor[['family', exact = TRUE]], family) ||
        !is.logical(descriptor[['complete', exact = TRUE]]) || length(descriptor$complete) != 1L ||
        is.na(descriptor$complete) || !.np_bw_procedure_values(descriptor) ||
        !identical(descriptor[['signature', exact = TRUE]], .np_bw_procedure_signature(bws)))
      stop('invalid or stale bandwidth selection.procedure; refit the bandwidth object', call. = FALSE)
  }
  if (!is.null(descriptor) && isTRUE(descriptor$complete)) {
    args <- descriptor[['arguments', exact = TRUE]]
    allowed <- c(names(.np_bw_procedure_kernel(bws)), 'bws', 'bandwidth.compute',
      'nmulti', 'powell.remin', 'bwsolver', 'itmax', 'ftol', 'tol', 'small',
      'invalid.penalty', 'penalty.multiplier', 'scale.factor.init.lower',
      'scale.factor.init.upper', 'scale.factor.init', 'lbd.init', 'hbd.init',
      'dfac.init', 'scale.init.categorical.sample', 'transform.bounds',
      'scale.factor.search.lower', 'cfac.dir', 'dfac.dir', 'dfc.dir', 'hbd.dir',
      'initc.dir', 'initd.dir', 'lbc.dir', 'lbd.dir', 'mads.nmulti',
      'nomad.nmulti', 'nomad.remin', 'nomad.opts',
      if (family == 'distribution') c('gdat', 'do.full.integral', 'ngrid', 'memfac'))
    if (!is.list(args) || is.null(names(args)) || anyDuplicated(names(args)) ||
        anyNA(names(args)) || any(!names(args) %in% allowed) ||
        !is.character(descriptor[['rule', exact = TRUE]]) || length(descriptor$rule) != 1L ||
        !descriptor$rule %in% c('manual', 'normal-reference', 'cv') ||
        !all(c('bws', 'bandwidth.compute', 'bwmethod', 'bwscaling', 'bwtype') %in% names(args)) ||
        !identical(args[['bandwidth.compute', exact = TRUE]], descriptor$rule != 'manual') ||
        !identical(args[['bwmethod', exact = TRUE]], bws[['method', exact = TRUE]]))
      stop('invalid bandwidth selection.procedure arguments', call. = FALSE)
    kernel <- .np_bw_procedure_kernel(bws)
    start <- args[['bws', exact = TRUE]]
    if (!identical(args[names(kernel)], kernel) || !is.numeric(start) ||
        length(start) != length(bws[['bw', exact = TRUE]]) || any(!is.finite(start)) ||
        (descriptor$rule == 'manual' && !identical(start, bws[['bw', exact = TRUE]])))
      stop('invalid bandwidth selection.procedure arguments', call. = FALSE)
    status <- 'retained-rule'; rule <- descriptor$rule
  } else if (is.null(descriptor) &&
             (identical(bws[['pmethod', exact = TRUE]], 'Manual') ||
              identical(bws[['method', exact = TRUE]], 'normal-reference'))) {
    args <- .np_bw_procedure_kernel(bws)
    rule <- if (identical(bws[['pmethod', exact = TRUE]], 'Manual')) 'manual' else 'normal-reference'
    args$bandwidth.compute <- rule != 'manual'
    args$bws <- if (rule == 'manual') bws[['bw', exact = TRUE]] else rep(0, length(bws[['bw', exact = TRUE]]))
    status <- 'legacy-rule-recovered'
  } else {
    args <- NULL; rule <- 'default-cv-compatibility'; status <- 'legacy-unresolved'
    messages <- 'Original bandwidth selection settings are unavailable; smoothed-percentile retains historical default-CV outer calibration. Refit the original bandwidth object for faithful outer selection.'
  }
  force(args); force(family)
  fit <- function(y) {
    if (is.null(args)) {
      # Outer columns already contain the generated observations. Their names
      # may be transformed formula terms; never evaluate those terms again.
      estimator <- if (family == 'density') npudens else npudist
      return(estimator(tdat = y))
    }
    # Keep selection and estimation inside one public entry when that API
    # expresses the retained rule. A named nonzero CV start is different:
    # estimator default methods interpret numeric bws as manual bandwidths.
    # Such starts therefore require the explicit selector below.
    combined <- !args$bandwidth.compute || args$bwmethod == 'normal-reference' ||
      all(args$bws == 0)
    if (combined) {
      fit.args <- args
      if (isTRUE(args$bandwidth.compute)) fit.args$bws <- NULL
      estimator <- if (family == 'density') npudens else npudist
      return(do.call(estimator, c(list(tdat = y), fit.args)))
    }
    # Native data entry avoids re-evaluating the original formula/subset/call.
    selector <- if (family == 'density') npudensbw else npudistbw
    selected <- do.call(selector, c(list(dat = y), args))
    if (family == 'density') npudens(bws = selected, tdat = y) else npudist(bws = selected, tdat = y)
  }
  list(fit = fit, messages = messages,
    diagnostics = list(status = status, rule = rule,
      objective = if (is.null(args)) if (family == 'density') 'cv.ml' else 'cv.cdf' else args$bwmethod,
      bandwidth.units = if (is.null(args)) 'default' else if (isTRUE(args$bwscaling)) 'scale-factor' else 'physical',
      original.rule.replayed = !is.null(args), descriptor.version = if (is.null(descriptor)) NULL else descriptor$version))
}
