# Private requested-only higher conditional LP uncertainty. Point owners and
# the incumbent scalar pseudo-response variance remain authoritative.
.np_conditional_native_call <- function(higher.se.request,
    tyuno, tyord, tycon, txuno, txord, txcon,
    eyuno, eyord, eycon, exuno, exord, excon,
    rbw, ymcv, ypadnum, xmcv, xpadnum, nconfac, ncatfac,
    mysd, myopti, enrow, xndim, ckerlbx, ckerubx, ckerlby, ckeruby,
    regtype, glp_degree, glp_bernstein, glp_basis, first_se, cat_se_request, se_request,
    allow.external = FALSE) {
  # Fixed signatures keep registration/count checks effective for both entries.
  if (is.null(higher.se.request))
    .Call("C_np_density_conditional",
          tyuno, tyord, tycon, txuno, txord, txcon,
          eyuno, eyord, eycon, exuno, exord, excon,
          rbw, ymcv, ypadnum, xmcv, xpadnum, nconfac, ncatfac,
          mysd, myopti, enrow, xndim, ckerlbx, ckerubx, ckerlby, ckeruby,
          regtype, glp_degree, glp_bernstein, glp_basis, first_se, cat_se_request, se_request,
          allow.external,
          PACKAGE = "np")
  else
    .Call("C_np_density_conditional_variance",
          tyuno, tyord, tycon, txuno, txord, txcon,
          eyuno, eyord, eycon, exuno, exord, excon,
          rbw, ymcv, ypadnum, xmcv, xpadnum, nconfac, ncatfac,
          mysd, myopti, enrow, xndim, ckerlbx, ckerubx, ckerlby, ckeruby,
          regtype, glp_degree, glp_bernstein, glp_basis, first_se, cat_se_request, se_request,
          allow.external,
          PACKAGE = "np")
}

.np_conditional_higher_se_request <- function(se, gradients, reg.engine,
                                               order, available, demand = NULL) {
  if (!se || !gradients || !identical(reg.engine, "lp") || !length(order))
    return(NULL)
  if (is.null(demand)) demand <- TRUE
  if (!is.logical(demand) || anyNA(demand) ||
      !(length(demand) %in% c(1L, length(order))))
    stop("invalid conditional higher-SE demand", call. = FALSE)
  selected <- available & order > 1L & rep_len(demand, length(order))
  if (any(selected)) selected else NULL
}

.np_conditional_higher_hat <- function(args, cdf = FALSE, base.rows = NULL,
                                        allow.external = FALSE, return.norm = TRUE) {
  bws <- args[["bws", exact = TRUE]]
  expected <- if (cdf) "condbandwidth" else "conbandwidth"
  where <- if (cdf) "npcdisthat" else "npcdenshat"
  if (!inherits(bws, expected))
    stop(sprintf("argument 'bws' must inherit from class '%s' in %s()",
                 expected, where), call. = FALSE)
  if (!is.null(base.rows)) {
    exdat <- args[["exdat", exact = TRUE]]
    if (!allow.external || is.null(exdat) || !is.integer(base.rows) ||
        length(base.rows) != nrow(exdat) || anyNA(base.rows) || any(!base.rows %in% 0:1))
      stop("internal conditional derivative row status mismatch", call. = FALSE)
    keep <- which(base.rows == 0L)
    value <- rep.int(NA_real_, length(base.rows))
    norm <- if (return.norm) matrix(NA_real_, length(base.rows), 3L) else NULL
    flags <- base.rows
    if (length(keep)) {
      selected <- args
      selected$exdat <- exdat[keep, , drop = FALSE]
      selected$eydat <- args[["eydat", exact = TRUE]][keep, , drop = FALSE]
      fitted <- .np_conditional_higher_hat(selected, cdf = cdf,
        allow.external = allow.external, return.norm = return.norm)
      rows <- if (return.norm) fitted[["value", exact = TRUE]] else fitted
      value[keep] <- rows
      extra <- attr(rows, ".np.empty.rows", exact = TRUE)
      if (!is.null(extra)) flags[keep] <- pmax(flags[keep], extra)
      if (return.norm) norm[keep, ] <- fitted[["norm", exact = TRUE]]
    }
    attr(value, ".np.empty.rows") <- flags
    return(if (return.norm) list(value = value, norm = norm) else value)
  }
  .npcdhat_core(
    bws = bws,
    txdat = args[["txdat", exact = TRUE]],
    tydat = args[["tydat", exact = TRUE]],
    exdat = args[["exdat", exact = TRUE]],
    eydat = args[["eydat", exact = TRUE]],
    y = args[["y", exact = TRUE]],
    output = "apply",
    operator = if (cdf) "integral" else "normal",
    x.s = args[["s", exact = TRUE]],
    class_name = where, where = where, return.norm = return.norm,
    allow.external = allow.external, .np.defer.empty.rows = TRUE)
}

.np_conditional_higher_se <- function(metadata, norm, nrows, base.rows = NULL) {
  if (!is.null(base.rows)) {
    if (!is.integer(base.rows) || length(base.rows) != nrows || anyNA(base.rows) ||
        any(!base.rows %in% 0:1) || !is.matrix(metadata) || !is.double(metadata) ||
        !identical(dim(metadata), c(as.integer(nrows), 4L)) ||
        !is.matrix(norm) || !is.double(norm) ||
        !identical(dim(norm), c(as.integer(nrows), 3L)))
      stop("internal conditional higher-SE row status mismatch", call. = FALSE)
    result <- rep.int(NA_real_, nrows)
    keep <- which(base.rows == 0L)
    if (length(keep)) result[keep] <- .np_conditional_higher_se(
      metadata[keep, , drop = FALSE], norm[keep, , drop = FALSE], length(keep))
    return(result)
  }
  if (!is.matrix(metadata) || !is.double(metadata) ||
      !identical(dim(metadata), c(as.integer(nrows), 4L)) ||
      !is.matrix(norm) || !is.double(norm) ||
      !identical(dim(norm), c(as.integer(nrows), 3L)) ||
      any(!is.finite(metadata)) || any(!is.finite(norm)) ||
      any(metadata[, 4L] != 1) || any(!metadata[, 3L] %in% 0:1) ||
      any(!norm[, 3L] %in% 0:1) || any(metadata[, 1L] < 0) ||
      any(norm[, 1L] < 0) || any(norm[, 2L] < 1))
    stop("internal conditional higher-SE row/ownership mismatch", call. = FALSE)
  result <- rep.int(NA_real_, nrows)
  valid <- metadata[, 3L] == 0 & norm[, 3L] == 0
  zero <- valid & (metadata[, 1L] == 0 | norm[, 1L] == 0)
  result[zero] <- 0
  rows <- which(valid & !zero)
  if (length(rows)) {
    # Do not square the norm or exp(2*a): either could overflow before the
    # final physical SE. Non-finite uncertainty is unavailable, not zero.
    log.se <- metadata[rows, 2L] + .5 * log(metadata[rows, 1L]) +
      log(norm[rows, 1L]) + .5 * log(norm[rows, 2L])
    value <- exp(log.se)
    value[!is.finite(value)] <- NA_real_
    result[rows] <- value
  }
  result
}
