# Metadata belonging to the public producer must precede immutable publication.
# Descriptors contain values already owned by that producer, never its frame or
# an executable method. The original formula wrapper reapplies this pure finish.
.npRmpi_bw_formula_apply <- function(value, plan) {
  if (is.null(plan)) return(value)
  if (identical(plan$producer, "nplsqregbw")) {
    # The vector combiner inherits its first scalar child's metadata. The
    # original formula wrapper appended aggregate terms after the tau fields.
    if (isTRUE(plan[["scalar.children", exact = TRUE]]) &&
        !isTRUE(plan[["scalar.child", exact = TRUE]]))
      value[["terms"]] <- NULL
    if (!isTRUE(plan[["scalar.child", exact = TRUE]]))
      value <- .nplsqreg_record_omit(value, plan$omit)
    value <- .nplsqreg_retain_formula_terms(value, plan$formula, plan$terms)
    value <- .nplsqreg_set_response_name(value, plan$ynames)
    if (!isTRUE(plan[["scalar.child", exact = TRUE]])) value$call <- plan$call
    return(value)
  }
  value$call <- plan$call
  value$formula <- plan$formula
  value$rows.omit <- as.vector(attr(plan$frame, "na.action"))
  value$nobs.omit <- length(value$rows.omit)
  value$terms <- plan$terms
  if (identical(plan$producer, "npplregbw")) value$xterms <- plan$xterms
  value <- .np_bws_retain_formula_training(value, plan$frame, plan$na.action)
  if (plan$producer %in% c("npplregbw", "npscoefbw"))
    value$chromoly <- plan$chromoly
  value <- updateBwNameMetadata(list(ynames = plan$ynames), value)
  if (identical(plan$producer, "npplregbw"))
    value$bw$yzbw <- updateBwNameMetadata(list(ynames = value$ynames), value$bw$yzbw)
  value
}

.npRmpi_bw_formula_plan <- function(producer, call, owner, formula,
                                    frame, na.action = NULL, terms = attr(frame, "terms"),
                                    xterms = NULL, chromoly = NULL,
                                    ynames = names(frame)[attr(terms, "response")],
                                    scalar.children = FALSE) {
  environment(call) <- owner
  if (identical(producer, "nplsqregbw"))
    call <- .nplsqreg_describe_call(call, owner)
  list(producer = producer, call = call, formula = formula, frame = frame,
       na.action = na.action, terms = terms, xterms = xterms,
       chromoly = chromoly, ynames = ynames, omit = attr(frame, "na.action"),
       scalar.children = scalar.children)
}

.npRmpi_bw_formula_finish <- function(value, plan) {
  complete <- FALSE
  on.exit({
    # A failed public finalizer must not strand an acknowledged publication.
    # No handle exists when the inner distributed producer itself failed.
    if (!complete) {
      fields <- .npRmpi_lease_handle_fields(value)
      if (!is.null(fields) && identical(.npRmpi_lease_state$status, "open"))
        .npRmpi_lease_run_lifecycle("autodispatch.lifecycle.rollback",
          list(ids = fields$id), comm = 1L,
          where = "public bandwidth finalization rollback")
    }
  }, add = TRUE)
  value <- .npRmpi_bw_formula_apply(value, plan)
  complete <- TRUE
  value
}

.npRmpi_bw_formula_owner_plan <- function(producer) {
  expected <- get0(paste0(producer, ".formula"), envir = asNamespace("npRmpi"),
                   mode = "function", inherits = FALSE)
  if (!is.function(expected)) return(NULL)
  parents <- sys.parents()
  at <- sys.parent()
  while (at > 0L) {
    if (identical(sys.function(at), expected)) {
      plan <- get0(".npRmpi.bw.metadata", envir = sys.frame(at), inherits = FALSE)
      # Vector-tau wrappers change only formula terms and labels on children;
      # each scalar child retains its own call and omission metadata.
      if (!is.null(plan) && identical(producer, "nplsqregbw") &&
          isTRUE(plan[["scalar.children", exact = TRUE]])) plan$scalar.child <- TRUE
      return(plan)
    }
    at <- parents[[at]]
  }
  NULL
}

.npRmpi_native_name_plan <- function(definition, owner, caller) {
  if (!is.function(definition)) return(list())
  out <- list()
  for (name in c("ynames", "xnames", "znames")) {
    roles <- switch(name, ynames = c("ydat", "tydat"),
      xnames = c("xdat", "txdat", "dat", "tdat"), znames = c("zdat", "tzdat"))
    role <- intersect(roles, names(formals(definition)))
    if (length(role)) {
      expr <- eval(substitute(substitute(ARG), list(ARG = as.name(role[[1L]]))), owner)
      out[name] <- list(.np_formula_dot_expression(expr, caller))
    }
  }
  out
}

.npRmpi_native_names_apply <- function(result, plan) {
  if (!is.list(result) || !length(plan)) return(result)
  bandwidth <- if (is.list(result[["bws", exact = TRUE]])) result[["bws"]] else result
  for (name in names(plan)) {
    label <- bandwidth[[name, exact = TRUE]]
    if (!is.character(label) || length(label) != 1L ||
        !startsWith(label, ".__npRmpi_autod_")) next
    label <- if (name == "ynames") deparse(plan[[name]]) else
      paste(deparse(plan[[name]]), collapse = "")
    bandwidth <- updateBwNameMetadata(setNames(list(label), name), bandwidth)
    if (inherits(bandwidth, "plbandwidth")) {
      if (name == "ynames")
        bandwidth$bw$yzbw <- updateBwNameMetadata(list(ynames = label), bandwidth$bw$yzbw)
      if (name == "xnames")
        for (i in seq_along(label))
          bandwidth$bw[[i + 1L]] <- updateBwNameMetadata(list(ynames = label[[i]]), bandwidth$bw[[i + 1L]])
      if (name == "znames")
        bandwidth$bw <- lapply(bandwidth$bw, function(bw)
          updateBwNameMetadata(list(xnames = label), bw))
    }
    if (inherits(bandwidth, "lsqregressionbandwidth"))
      bandwidth <- .nplsqreg_set_response_name(bandwidth, bandwidth$ynames,
        predictor.names = if (name == "xnames") label else NULL)
    if (is.list(result[["bws", exact = TRUE]])) {
      result[["bws"]] <- bandwidth
      if (!is.null(result[[name, exact = TRUE]])) result[[name]] <- label
    } else result <- bandwidth
  }
  result
}

.npRmpi_bw_metadata_wire <- function(metadata) {
  # Formula evaluation belongs to the coordinator that already made the model
  # frame. Worker leases retain its descriptive AST and data, not that caller's
  # unrelated objects. bandwidth.v1 already treats these environments this way;
  # the coordinator reapplies the original metadata before attaching its handle.
  project <- function(x) {
    if (is.call(x)) {
      attrs <- attributes(x)
      parts <- as.list(x)
      for (i in seq_along(parts))
        if (!.npRmpi_is_missing_call_arg(parts[[i]]))
          parts[i] <- list(project(parts[[i]]))
      x <- as.call(parts)
    } else if (is.pairlist(x)) {
      attrs <- attributes(x)
      parts <- as.list(x)
      for (i in seq_along(parts))
        if (!.npRmpi_is_missing_call_arg(parts[[i]]))
          parts[i] <- list(project(parts[[i]]))
      x <- as.pairlist(parts)
    } else if (is.list(x)) {
      attrs <- attributes(x)
      for (i in seq_along(x)) x[i] <- list(project(x[[i]]))
    } else attrs <- attributes(x)
    if (length(attrs)) {
      if (is.call(x) || inherits(x, "formula") || inherits(x, "terms"))
        if (".Environment" %in% names(attrs)) attrs[[".Environment"]] <- emptyenv()
      attrs <- lapply(attrs, project)
      attributes(x) <- attrs
    }
    x
  }
  project(metadata)
}

.npRmpi_bw_restore_public_owner <- function(result, originals) {
  if (!length(originals) || !is.list(result)) return(result)
  bws <- result[["bws", exact = TRUE]]
  if (!is.list(bws)) return(result)
  for (original in originals) {
    if (is.null(original[["formula", exact = TRUE]])) next
    if (!inherits(original, c("rbandwidth", "sibandwidth", "scbandwidth",
                             "plbandwidth", "lsqregressionbandwidth"))) next
    restore <- function(value, reference) {
      fields <- c("call", "formula", "terms", "xterms", ".np.formula.training")
      # Fitting may legitimately update omission/diagnostic metadata. Verify
      # the constructor-owned formula representation itself, and restore only
      # that representation; never replace fitted state or numerical fields.
      current <- value[intersect(names(value), fields)]
      expected <- reference[intersect(names(reference), fields)]
      if (!identical(.npRmpi_autodispatch_fingerprint(current),
                     .npRmpi_autodispatch_fingerprint(expected)))
        stop("retained formula metadata changed before public owner restoration",
             call. = FALSE)
      for (name in names(expected)) value[name] <- expected[name]
      if (!is.null(value[["tau.bws", exact = TRUE]]))
        for (i in seq_along(value$tau.bws))
          value$tau.bws[[i]] <- restore(value$tau.bws[[i]], reference$tau.bws[[i]])
      value
    }
    result[["bws"]] <- restore(bws, original)
  }
  result
}
