# Single-index coordinates are a parametric design, not mixed-kernel codes.
# Keep the raw training data in the existing retained-data owner. Only new
# bandwidth objects with factors carry this schema; old objects retain their
# original score interpretation. No work here belongs inside an objective.
.np_index_formula_xdat <- function(frame) {
  tt <- attr(frame, "terms")
  xdat <- frame[, .np_formula_term_names(attr(tt, "term.labels")), drop = FALSE]
  if (any(vapply(xdat, is.factor, logical(1L))))
    attr(xdat, ".np.index.terms") <- delete.response(tt)
  xdat
}

.np_index_design_train <- function(xdat, ydat) {
  categorical <- vapply(xdat, is.factor, logical(1L))
  if (!any(categorical)) return(NULL)
  categorical <- categorical | vapply(xdat, is.logical, logical(1L))
  .np_require_paired_rows(xdat, ydat, "xdat", "ydat")
  keep <- complete.cases(xdat, ydat)
  if (!any(keep)) stop("Data has no rows without NAs", call. = FALSE)
  frame <- xdat[keep, , drop = FALSE]
  for (name in names(frame)[categorical]) {
    original <- frame[[name]]
    if (is.logical(original)) original <- factor(original, levels = c(FALSE, TRUE))
    frame[[name]] <- droplevels(original)
    if (nlevels(frame[[name]]) < 2L)
      stop(sprintf("npindex: factor '%s' needs at least two observed training levels", name),
           call. = FALSE)
    specified <- attr(original, "contrasts", exact = TRUE)
    if (!is.null(specified)) {
      if (is.matrix(specified))
        specified <- specified[match(levels(frame[[name]]), levels(original)), , drop = FALSE]
      attr(frame[[name]], "contrasts") <- specified
    }
  }
  tt <- attr(xdat, ".np.index.terms", exact = TRUE)
  if (is.null(tt)) {
    # Construct a formula from symbols, including non-syntactic native names.
    rhs <- Reduce(function(a, b) call("+", a, b), lapply(names(xdat), as.name))
    tt <- terms(as.formula(call("~", rhs), env = baseenv()))
  }
  attr(frame, "terms") <- tt
  contrasts <- lapply(frame[categorical], stats::contrasts)
  mm <- model.matrix(tt, frame, contrasts.arg = contrasts)
  selected <- colnames(mm) != "(Intercept)"
  design <- mm[, selected, drop = FALSE]
  if (!ncol(design) || qr(cbind(1, design))$rank != ncol(design) + 1L)
    stop(paste0("npindex: factor design is aliased or contains an unidentified constant direction. ",
                "Use an ordinary intercept formula (the index absorbs its intercept) ",
                "or supply independent numeric design columns."), call. = FALSE)
  schema <- list(version = 1L, raw.names = names(xdat),
       factor.names = names(frame)[categorical],
       ordered = vapply(frame[categorical], is.ordered, logical(1L)),
       levels = lapply(frame[categorical], levels), contrasts = contrasts,
       terms = tt, columns = colnames(design),
       assign = attr(mm, "assign")[selected],
       factor.columns = attr(mm, "assign")[selected] %in% which(categorical))
  prepared <- matrix(NA_real_, nrow(xdat), ncol(design),
                     dimnames = list(row.names(xdat), colnames(design)))
  prepared[keep, ] <- design
  list(schema = schema, data = as.data.frame(prepared, optional = TRUE))
}

.np_index_design_apply <- function(design, xdat) {
  if (is.null(design)) return(xdat)
  if (!identical(design[["version", exact = TRUE]], 1L))
    stop("npindex: unsupported retained index design version", call. = FALSE)
  missing <- setdiff(design$raw.names, names(xdat))
  if (length(missing))
    stop(sprintf("npindex: missing predictor(s): %s", paste(missing, collapse = ", ")),
         call. = FALSE)
  frame <- xdat[, design$raw.names, drop = FALSE]
  for (name in design$factor.names) {
    values <- as.character(frame[[name]])
    unseen <- setdiff(values[!is.na(values)], design$levels[[name]])
    if (length(unseen))
      stop(sprintf("npindex: new level(s) for '%s': %s", name,
                   paste(unseen, collapse = ", ")), call. = FALSE)
    frame[[name]] <- factor(values, levels = design$levels[[name]],
                            ordered = design$ordered[[name]])
  }
  numeric.names <- setdiff(design$raw.names, design$factor.names)
  if (any(!vapply(frame[numeric.names], is.numeric, logical(1L))))
    stop("npindex: numeric predictor types must match the training design", call. = FALSE)
  attr(frame, "terms") <- design$terms
  mm <- model.matrix(design$terms, frame, contrasts.arg = design$contrasts)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (!identical(colnames(mm), design$columns))
    stop("npindex: evaluation design columns do not match training", call. = FALSE)
  as.data.frame(mm, optional = TRUE)
}

.np_index_predictor_names <- function(bws) {
  design <- bws[["index.design", exact = TRUE]]
  if (is.null(design)) bws$xnames else design$raw.names
}

# Resolve contrasts in the calling session, not independently on each rank.
# The private transport cache is consumed by the constructor and never retained
# in the public training frame. Numeric-only calls keep the original payload.
.np_index_dispatch <- function(mc, xdat, ydat, caller.env, owner.name) {
  original.call <- mc
  xdat <- toFrame(xdat)
  prepared <- attr(xdat, ".np.index.prepared", exact = TRUE)
  if (is.null(prepared)) prepared <- .np_index_design_train(xdat, ydat)
  if (is.null(prepared))
    return(.npRmpi_autodispatch_call(mc, caller.env, owner.name = owner.name))
  attr(xdat, ".np.index.prepared") <- prepared
  mc$xdat <- xdat
  mc$ydat <- ydat
  result <- .npRmpi_autodispatch_call(mc, caller.env, owner.name = owner.name)
  # The dispatcher materializes arguments in returned calls. Restore the
  # caller's expression so the private transport cache cannot escape with it.
  result$call <- original.call
  environment(result$call) <- caller.env
  result
}
