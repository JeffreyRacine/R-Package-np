# An internal secondary fit cannot consume a conditional bandwidth lease.
# Certify only its immediate package-owned selector/formula-dispatch chain;
# user wrappers, callback ancestors and previously returned bandwidths are not
# evidence of temporary ownership. The caller carries this decision once.
.npRmpi_bw_publication_temporary <- function(producer, method.frame) {
  consumer <- switch(producer, npcdensbw = "npconmode.default",
                     npcdistbw = "npqreg.default")
  if (is.null(consumer)) return(FALSE)
  ns <- asNamespace("npRmpi")
  generic <- get0(producer, ns, mode = "function", inherits = FALSE)
  frames <- sys.frames()
  k <- which(vapply(frames, identical, logical(1L), method.frame))
  if (length(k) != 1L || k < 4L ||
      !identical(sys.function(k - 1L), generic)) return(FALSE)
  methods <- paste0(producer, c(".formula", ".NULL", ".default"))
  if (!any(vapply(methods, function(name)
      identical(sys.function(k), get0(name, ns, mode = "function", inherits = FALSE)),
      logical(1L)))) return(FALSE)
  owner <- NULL
  if (k >= 5L && is.primitive(sys.function(k - 2L)) &&
      .np_bw_selector_call_identical(sys.call(k - 2L), sys.call(k - 3L)) &&
      sys.parents()[[k - 2L]] == k - 3L &&
      identical(sys.function(k - 3L), base::eval) &&
      identical(sys.function(k - 4L), .np_eval_bw_call)) {
    holder <- frames[[k - 4L]]
    value.env <- get0("value.env", holder, inherits = FALSE)
    if (!is.environment(value.env) ||
        !identical(frames[[k - 2L]], value.env)) return(FALSE)
    owner <- get0("native.frame", holder, inherits = FALSE)
  } else if (identical(sys.function(k - 2L), base::do.call) &&
             identical(sys.function(k - 3L), .np_formula_dispatch_call)) {
    holder <- frames[[k - 3L]]
    if (!identical(get0("fun", holder, inherits = FALSE), generic)) return(FALSE)
    owner <- get0("envir", holder, inherits = FALSE)
    if (!is.environment(owner) ||
        !identical(frames[[sys.parents()[[k]]]], owner)) return(FALSE)
  }
  if (!is.environment(owner)) return(FALSE)
  at <- which(vapply(frames, identical, logical(1L), owner))
  length(at) == 1L && identical(sys.function(at),
    get0(consumer, ns, mode = "function", inherits = FALSE))
}

# Metadata belonging to the public producer must precede immutable publication.
# Descriptors contain values already owned by that producer, never its frame or
# an executable method. The original formula wrapper reapplies this pure finish.
.npRmpi_bw_formula_apply <- function(value, plan) {
  if (is.null(plan)) return(value)
  if (plan$producer %in% c("npudensbw", "npudistbw")) {
    value$call <- plan$call
    value$formula <- plan$formula
    value$terms <- plan$terms
    value <- .np_bws_retain_formula_training(value, plan$frame, plan$na.action)
    value$rows.omit <- as.vector(attr(plan$frame, "na.action"))
    value$nobs.omit <- length(value$rows.omit)
    return(value)
  }
  if (plan$producer %in% c("npcdensbw", "npcdistbw")) {
    value$call <- plan$call
    if (is.null(plan[["formula", exact = TRUE]])) return(value)
    value$formula <- plan$formula
    value$rows.omit <- as.vector(attr(plan$frame, "na.action"))
    value$nobs.omit <- length(value$rows.omit)
    value$terms <- plan$terms
    value <- .np_bws_retain_formula_training(value, plan$frame, plan$na.action)
    value$variableNames <- plan$variableNames
    return(value)
  }
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
  call <- .np_public_replay_call(call, producer, frame = parent.frame())
  call <- .np_bw_selector_formula_call(call, producer, parent.frame())
  environment(call) <- owner
  if (identical(producer, "nplsqregbw"))
    environment(call) <- .np_call_owner_environment(owner)
  plan <- list(producer = producer, call = call, formula = formula, frame = frame,
       na.action = na.action, terms = terms, xterms = xterms,
       chromoly = chromoly, ynames = ynames, omit = attr(frame, "na.action"),
       scalar.children = scalar.children)
  if (producer %in% c("npcdensbw", "npcdistbw"))
    plan$temporary <- .npRmpi_bw_publication_temporary(producer, parent.frame())
  plan
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
  native <- if (producer %in% c("npcdensbw", "npcdistbw"))
    get0(paste0(producer, ".NULL"), envir = asNamespace("npRmpi"),
         mode = "function", inherits = FALSE) else NULL
  native.plan <- NULL
  parents <- sys.parents()
  at <- sys.parent()
  while (at > 0L) {
    if (is.function(native) && identical(sys.function(at), native))
      native.plan <- get0(".npRmpi.bw.metadata", envir = sys.frame(at),
                         inherits = FALSE)
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
  native.plan
}

.npRmpi_native_label_pack <- function(label) {
  # Literal vector arguments can have long public names. Preserve the exact
  # coordinator text without sending a larger representation unnecessarily.
  bytes <- serialize(label, NULL)
  if (length(bytes) < 512L || length(bytes) > .Machine$integer.max)
    return(label)
  payload <- tryCatch(memCompress(bytes, type = "bzip2"), error = function(e) NULL)
  if (is.null(payload)) return(label)
  packed <- list(codec = "bzip2", payload = payload)
  if (length(serialize(packed, NULL)) < length(bytes)) packed else label
}

.npRmpi_native_label_unpack <- function(label) {
  if (!is.character(label)) {
    if (!is.list(label) || !identical(names(label), c("codec", "payload")) ||
        !identical(label[["codec", exact = TRUE]], "bzip2") || !is.raw(label[["payload", exact = TRUE]]))
      stop("invalid native label metadata", call. = FALSE)
    label <- unserialize(memDecompress(label[["payload", exact = TRUE]], type = "bzip2"))
  }
  if (!is.character(label) || !length(label) ||
      !is.null(attributes(label)) || anyNA(label))
    stop("invalid native label text", call. = FALSE)
  label
}

.npRmpi_native_name_plan <- function(definition, owner, caller, values,
                                    identity.roles = character()) {
  if (!is.function(definition) ||
      inherits(values[["bws", exact = TRUE]], c("bandwidth", "dbandwidth")))
    return(list())
  out <- list()
  for (name in c("ynames", "xnames", "znames")) {
    roles <- switch(name, ynames = c("ydat", "tydat"),
      xnames = c("xdat", "txdat", "dat", "tdat"), znames = c("zdat", "tzdat"))
    role <- intersect(roles, names(formals(definition)))
    if (length(role) && !role[[1L]] %in% identity.roles &&
        role[[1L]] %in% names(values) &&
        is.atomic(values[[role[[1L]]]]) && is.null(dim(values[[role[[1L]]]]))) {
      expr <- eval(substitute(substitute(ARG), list(ARG = as.name(role[[1L]]))), owner)
      expr <- .np_formula_dot_expression(expr, caller)
      # Public names belong to the original caller. Worker formatting options
      # must not affect the object that the publication certificate compares.
      label <- if (name == "ynames") deparse(expr) else
        paste(deparse(expr), collapse = "")
      out[name] <- list(.npRmpi_native_label_pack(label))
    }
  }
  out
}

.npRmpi_native_names_apply <- function(result, plan) {
  if (!is.list(result) || !length(plan)) return(result)
  if (inherits(result, c("bandwidth", "dbandwidth")) && "xnames" %in% names(plan))
    return(.npRmpi_unconditional_names_apply(result,
      list(dat = .npRmpi_native_label_unpack(plan[["xnames"]]))))
  bandwidth <- if (is.list(result[["bws", exact = TRUE]])) result[["bws"]] else result
  for (name in names(plan)) {
    label <- bandwidth[[name, exact = TRUE]]
    if (!is.character(label) || length(label) != 1L ||
        !startsWith(label, ".__npRmpi_autod_")) next
    label <- .npRmpi_native_label_unpack(plan[[name]])
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

# Literal metadata may reuse data already sent in this invocation. Shape and
# original call fingerprints are certified before replacing any literal. These
# private references are never evaluated as R expressions or retained in results.
.npRmpi_bw_literal_plain <- function(x) {
  attrs <- attributes(x)
  if (is.data.frame(x)) attrs[["row.names"]] <- .row_names_info(x, 0L)
  out <- if (is.list(x)) lapply(x, .npRmpi_bw_literal_plain) else {
    out <- vector(typeof(x), length(x))
    out[] <- unclass(x)
    out
  }
  attributes(out) <- NULL
  attributes(out) <- attrs
  out
}

.npRmpi_bw_literal_frame <- function(x) {
  attrs <- attributes(x)
  identical(attr(x, "class", exact = TRUE), "data.frame") &&
    length(attrs) <= 8L &&
    all(vapply(x, .npRmpi_transport_plain_value, logical(1L))) &&
    all(vapply(attrs, function(a) is.atomic(a) && is.null(attributes(a)),
               logical(1L)))
}

.npRmpi_bw_literal_rebuild <- function(source, code) {
  bad <- function() stop("invalid internal literal data reference", call. = FALSE)
  if (!is.data.frame(source) || !is.character(code) || length(code) != 1L ||
      is.na(code) || nchar(code, type = "bytes") > 32L) bad()
  if (grepl("^F[01]:[1-8](,[1-8]){0,7}$", code)) {
    parts <- strsplit(substring(code, 2L), ":", fixed = TRUE)[[1L]]
    order <- as.integer(strsplit(parts[[2L]], ",", fixed = TRUE)[[1L]])
    attrs <- attributes(source)
    if (!.npRmpi_bw_literal_frame(source) || anyDuplicated(order) ||
        length(order) != length(attrs) || any(order > length(attrs))) bad()
    out <- .npRmpi_bw_literal_plain(source)
    attrs[["row.names"]] <- if (identical(parts[[1L]], "1"))
      c(NA_integer_, -nrow(source)) else .row_names_info(source, 0L)
    attributes(out) <- NULL
    attributes(out) <- attrs[order]
    return(out)
  }
  if (!.npRmpi_transport_plain_value(source)) bad()
  if (code %in% c("C0", "C1")) {
    if (length(source) != 1L) bad()
    out <- .npRmpi_bw_literal_plain(.subset2(source, 1L))
    if (identical(code, "C1")) names(out) <- attr(source, "row.names")
    return(out)
  }
  if (!code %in% paste0("M", 0:3) || !length(source)) bad()
  flag <- as.integer(substring(code, 2L))
  out <- unlist(lapply(source, function(x)
    .npRmpi_bw_literal_plain(unclass(x))),
    use.names = FALSE)
  dim(out) <- c(nrow(source), length(source))
  if (flag > 0L) dimnames(out) <- list(
    if (bitwAnd(flag, 1L)) attr(source, "row.names") else NULL,
    if (bitwAnd(flag, 2L)) names(source) else NULL)
  out
}

.npRmpi_bw_literal_code <- function(value, source) {
  if (!is.data.frame(source)) return(NULL)
  if (is.data.frame(value)) {
    if (!.npRmpi_bw_literal_frame(value) ||
        !.npRmpi_bw_literal_frame(source) ||
        !setequal(names(attributes(value)), names(attributes(source)))) return(NULL)
    code <- paste0("F", as.integer(.row_names_info(value, 1L) < 0L), ":",
      paste(match(names(attributes(value)), names(attributes(source))), collapse = ","))
  } else {
    if (!.npRmpi_transport_plain_value(value) ||
        !.npRmpi_transport_plain_value(source)) return(NULL)
    if (is.matrix(value)) {
      if (!length(source)) return(NULL)
      if (any(!names(attributes(value)) %in% c("dim", "dimnames"))) return(NULL)
      dn <- dimnames(value)
      flag <- if (is.null(dn)) 0L else
        as.integer(!is.null(dn[[1L]])) + 2L * as.integer(!is.null(dn[[2L]]))
      code <- paste0("M", flag)
    } else if (is.atomic(value) && is.null(dim(value)) && length(source) == 1L)
      code <- if (is.null(names(value))) "C0" else "C1"
    else return(NULL)
  }
  # Hash the literal in its actual call position, before reconstruction reads.
  wrap <- function(x) as.call(list(as.name(".np.literal"), value = x))
  before <- .npRmpi_autodispatch_fingerprint(wrap(value))
  if (is.na(before)) return(NULL)
  out <- .npRmpi_bw_literal_rebuild(source, code)
  if (!identical(before, .npRmpi_autodispatch_fingerprint(wrap(out))) ||
      !identical(value, out, num.eq = FALSE, single.NA = FALSE,
                 attrib.as.set = FALSE) ||
      !.npRmpi_transport_same_encoding(value, out)) return(NULL)
  code
}

.npRmpi_bw_literals_encode <- function(metadata, original, call, values) {
  refs <- character()
  for (target in c("f", "s")) for (role in c("xdat", "ydat")) {
    value <- if (identical(target, "f")) original$formula$call[[role]] else
      original[["selector.native", exact = TRUE]][[role]]
    if (is.null(value) || is.language(value)) next
    token <- call[[role]]
    if (!is.symbol(token)) next
    token <- as.character(token)
    if (!token %in% names(values)) next
    code <- .npRmpi_bw_literal_code(value, values[[token]])
    if (is.null(code)) next
    refs[paste0(target, ".", role)] <- code
    if (identical(target, "f")) metadata$formula$call[role] <- list(NULL)
    else metadata[["selector.native"]][role] <- list(NULL)
  }
  if (length(refs)) metadata$literal <- refs
  metadata
}

.npRmpi_bw_literals_complete <- function(metadata) {
  for (target in c("f", "s")) for (role in c("xdat", "ydat")) {
    value <- if (identical(target, "f"))
      metadata[["formula", exact = TRUE]][["call", exact = TRUE]][[role]] else
      metadata[["selector.native", exact = TRUE]][[role]]
    if (!is.null(value) && !is.language(value)) return(FALSE)
  }
  TRUE
}

.npRmpi_bw_literals_decode <- function(metadata, call, values) {
  refs <- metadata[["literal", exact = TRUE]]
  if (is.null(refs)) return(metadata)
  if (!is.character(refs) || !length(refs) || length(refs) > 4L ||
      anyNA(refs) || is.null(names(refs)) || anyNA(names(refs)) ||
      anyDuplicated(names(refs)) ||
      any(!names(refs) %in% c("f.xdat", "f.ydat", "s.xdat", "s.ydat")))
    stop("invalid internal literal metadata", call. = FALSE)
  for (name in names(refs)) {
    target <- substr(name, 1L, 1L)
    role <- substring(name, 3L)
    slot <- if (identical(target, "f")) metadata$formula$call else
      metadata[["selector.native", exact = TRUE]]
    if ((identical(target, "f") && !is.call(slot)) ||
        (identical(target, "s") && !is.list(slot)) ||
        !role %in% names(slot) || !is.null(slot[[role]]))
      stop("invalid internal literal target", call. = FALSE)
    token <- call[[role]]
    if (!is.symbol(token) || !as.character(token) %in% names(values))
      stop("internal literal data reference is unavailable", call. = FALSE)
    out <- .npRmpi_bw_literal_rebuild(values[[as.character(token)]], refs[[name]])
    if (identical(target, "f")) metadata$formula$call[role] <- list(out)
    else metadata[["selector.native"]][role] <- list(out)
  }
  metadata$literal <- NULL
  metadata
}

.npRmpi_bw_metadata_wire <- function(metadata, plain.language = FALSE) {
  # The coordinator has already used this local owner decision; it is not
  # descriptive metadata and does not belong in the worker payload.
  if (is.list(metadata$formula)) metadata$formula$temporary <- NULL
  # Formula evaluation belongs to the coordinator that already made the model
  # frame. Worker leases retain its descriptive AST and data, not that caller's
  # unrelated objects. bandwidth.v1 already treats these environments this way;
  # the coordinator reapplies the original metadata before attaching its handle.
  project <- function(x) {
    if (is.call(x)) {
      attrs <- attributes(x)
      parts <- as.list(if (plain.language && is.object(x)) unclass(x) else x)
      for (i in seq_along(parts))
        if (!.npRmpi_is_missing_call_arg(parts[[i]]))
          parts[i] <- list(project(parts[[i]]))
      x <- as.call(parts)
    } else if (is.pairlist(x)) {
      attrs <- attributes(x)
      parts <- as.list(if (plain.language && is.object(x)) unclass(x) else x)
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
  out <- project(metadata)
  # LSQ carries its omission record and formula terms separately. Its metadata
  # finisher does not consume this frame; preserve projection callbacks first.
  if (is.list(out) && !is.object(out)) {
    plan <- .subset2(out, "formula", exact = TRUE)
    if (is.list(plan) && !is.object(plan) &&
        identical(.subset2(plan, "producer", exact = TRUE), "nplsqregbw")) {
      plan[["frame"]] <- NULL
      out[["formula"]] <- plan
    }
  }
  out
}

.npRmpi_bw_restore_public_owner <- function(result, originals) {
  if (!length(originals) || !is.list(result)) return(result)
  bws <- result[["bws", exact = TRUE]]
  if (!is.list(bws)) return(result)
  for (original in originals) {
    if (is.null(original[["formula", exact = TRUE]])) next
    if (!inherits(original, c("rbandwidth", "sibandwidth", "scbandwidth",
                             "plbandwidth", "lsqregressionbandwidth",
                             "conbandwidth", "condbandwidth",
                             "bandwidth", "dbandwidth"))) next
    fields <- c("call", "formula", "terms", "xterms", ".np.formula.training")
    metadata <- function(x) {
      out <- x[intersect(names(x), fields)]
      training <- out[[".np.formula.training", exact = TRUE]]
      # Native fitting may replace the effective sample. Keep its presence
      # and all retained policies in the guard, but not the sample values.
      if (is.list(training) && "frame" %in% names(training)) {
        training["frame"] <- list(NULL)
        out[".np.formula.training"] <- list(training)
      }
      out
    }
    restore <- function(value, reference,
                        value.key = .npRmpi_autodispatch_fingerprint(metadata(value)),
                        reference.key = .npRmpi_autodispatch_fingerprint(metadata(reference))) {
      fields <- c("call", "formula", "terms", "xterms", ".np.formula.training")
      if (!identical(value.key, reference.key))
        stop("retained formula metadata changed before public owner restoration",
             call. = FALSE)
      restored <- reference[intersect(names(reference), fields)]
      training <- restored[[".np.formula.training", exact = TRUE]]
      if (is.list(training) && "frame" %in% names(training)) {
        frame <- value[[".np.formula.training", exact = TRUE]][["frame", exact = TRUE]]
        tt <- attr(frame, "terms", exact = TRUE)
        if (inherits(tt, "terms") && identical(environment(tt), emptyenv())) {
          owner <- environment(attr(training[["frame", exact = TRUE]], "terms", exact = TRUE))
          if (is.null(owner)) owner <- environment(reference[["formula", exact = TRUE]])
          environment(tt) <- owner
          attr(frame, "terms") <- tt
        }
        training["frame"] <- list(frame)
        restored[".np.formula.training"] <- list(training)
      }
      for (name in names(restored)) value[name] <- restored[name]
      if (!is.null(value[["tau.bws", exact = TRUE]]))
        for (i in seq_along(value$tau.bws))
          value$tau.bws[[i]] <- restore(value$tau.bws[[i]], reference$tau.bws[[i]])
      value
    }
    # Keep the guard order, and reuse only the reference's immutable digest.
    value.key <- .npRmpi_autodispatch_fingerprint(metadata(bws))
    reference.key <- .npRmpi_autodispatch_fingerprint(metadata(original))
    result[["bws"]] <- restore(bws, original, value.key, reference.key)
    # Partially linear fits expose a second bandwidth view. Validate it in its
    # own right and preserve its effective sample and non-metadata fields.
    alias <- result[["bw", exact = TRUE]]
    if (inherits(result, "plregression") && inherits(alias, "plbandwidth") &&
        inherits(original, "plbandwidth"))
      result[["bw"]] <- restore(alias, original, reference.key = reference.key)
  }
  result
}

# Observe a native vector only after the existing argument owner has forced it.
# Formula owners and caller-supplied bandwidth objects retain their own identity.
.npRmpi_native_identity_plan <- function(definition, mc, owner, caller) {
  if (!is.function(definition) || !is.call(mc) ||
      is.list(mc[["bws"]]) || is.list(mc[["bw"]])) return(list())
  out <- list()
  for (role in intersect(c("xdat", "txdat", "ydat", "tydat", "zdat", "tzdat"),
                         intersect(names(formals(definition)), names(mc)))) {
    value <- mc[[role]]
    if (!is.atomic(value) || !is.null(dim(value))) next
    expr <- eval(substitute(substitute(ARG), list(ARG = as.name(role))), owner)
    expr <- .np_formula_dot_expression(expr, caller)
    out[[role]] <- paste(deparse(expr), collapse = "")
  }
  out
}

.npRmpi_native_identity_bind <- function(plan, call, temporaries) {
  out <- character()
  for (role in names(plan)) {
    token <- call[[role]]
    if (is.symbol(token) && as.character(token) %in% temporaries)
      out[as.character(token)] <- plan[[role]]
  }
  out
}

# Only descriptive fields in these producer schemas are rewritten, before the
# existing publication. Calls, formulas, values and acknowledged leases are not.
.npRmpi_native_identity_apply <- function(value, identity) {
  if (!is.list(value) || !length(identity)) return(value)
  if (!is.null(.npRmpi_lease_handle_fields(value))) return(value)
  rename <- function(x) {
    if (!is.character(x)) return(x)
    at <- match(x, names(identity))
    ok <- !is.na(at)
    x[ok] <- unname(identity[at[ok]])
    x
  }
  named <- function(x) {
    if (!is.null(names(x))) names(x) <- rename(names(x))
    if (!is.data.frame(x) && !is.null(dimnames(x)))
      dimnames(x) <- lapply(dimnames(x), rename)
    x
  }
  descriptors <- function(x) {
    if (!is.list(x)) return(named(x))
    lapply(x, function(a) if (is.list(a)) descriptors(a) else named(a))
  }
  for (name in c("xnames", "ynames", "znames", "data.xnames", "data.znames"))
    if (!is.null(value[[name, exact = TRUE]])) value[[name]] <- rename(value[[name]])
  if (is.list(value[["varnames", exact = TRUE]]))
    value$varnames <- lapply(value$varnames, rename)
  for (name in c("xcoef", "xcoeferr", "xcoefvcov", "evalx", "evalz", "xeval", "yeval",
                 "icon", "iuno", "iord", "ixcon", "ixuno", "ixord", "iycon", "iyuno", "iyord", "sdev"))
    if (!is.null(value[[name, exact = TRUE]])) value[[name]] <- named(value[[name]])
  for (name in c("xdati", "ydati", "zdati", "dati"))
    if (is.list(value[[name, exact = TRUE]])) value[[name]] <- descriptors(value[[name]])
  if (is.list(value[["eval", exact = TRUE]]))
    for (name in c("exdat", "ezdat"))
      if (!is.null(value$eval[[name, exact = TRUE]])) value$eval[[name]] <- named(value$eval[[name]])
  training <- value[[".np.native.training", exact = TRUE]]
  if (is.list(training)) {
    for (name in c("xdat", "ydat", "zdat"))
      if (!is.null(training[[name, exact = TRUE]])) training[[name]] <- named(training[[name]])
    value[[".np.native.training"]] <- training
  }
  for (name in c("bws", "reg.bws", "fit", "mean.fit", "scale.fit"))
    if (is.list(value[[name, exact = TRUE]]))
      value[[name]] <- .npRmpi_native_identity_apply(value[[name]], identity)
  for (name in c("tau.bws", "tau.fits"))
    if (is.list(value[[name, exact = TRUE]]))
      value[[name]] <- lapply(value[[name]], .npRmpi_native_identity_apply, identity = identity)
  bw <- value[["bw", exact = TRUE]]
  if (is.list(bw)) {
    if (inherits(bw, "plbandwidth") || !is.null(bw[["xnames", exact = TRUE]]))
      value$bw <- .npRmpi_native_identity_apply(bw, identity)
    else value$bw <- lapply(bw, .npRmpi_native_identity_apply, identity = identity)
  }
  value
}

# Certify the immediate native selector before promise binding. Conditional
# producers share regression's ownership check; their labels are carried from
# the coordinator so workers never re-deparse an expression under local options.
.npRmpi_bw_selector_name_plan <- function(expressions, owner) {
  roles <- intersect(names(expressions), c("xdat", "ydat"))
  prepared <- get("bw.args", envir = owner, inherits = FALSE)
  needed <- roles[vapply(roles, function(name) {
    # Both certified default owners have already resolved their public labels.
    # Reuse them without querying a custom frame's names accessor again.
    field <- switch(name, xdat = "xnames", ydat = "ynames")
    identical(prepared[[field, exact = TRUE]], paste0(".np_native_", name))
  }, logical(1L))]
  if (!length(needed)) return(NULL)
  # Conditional producers take labels from prepared column names. Preserve the
  # rare explicit sentinel name, without deparsing ordinary literal samples.
  setNames(lapply(needed, function(name) {
    text <- deparse(expressions[[name]])
    if (identical(name, "ydat")) text else paste(text, collapse = "")
  }), needed)
}

.npRmpi_bw_native_selector_plan <- function(mc, owner.name, definition,
                                            owner, caller) {
  if (!is.character(owner.name) || length(owner.name) != 1L || is.na(owner.name))
    return(NULL)
  producer <- sub("[.]default$", "", owner.name)
  if (!producer %in% c("npregbw", "npcdensbw", "npcdistbw") ||
      !identical(owner.name, paste0(producer, ".default")) ||
      !identical(definition, get0(owner.name, envir = asNamespace("npRmpi"),
                                 mode = "function", inherits = FALSE))) return(NULL)
  generic <- get0(producer, envir = asNamespace("npRmpi"),
                  mode = "function", inherits = FALSE)
  frames <- sys.frames()
  parents <- sys.parents()
  k <- which(vapply(frames, identical, logical(1L), owner))
  if (length(k) != 1L || k < 5L ||
      !identical(sys.function(k - 1L), generic) ||
      !is.primitive(sys.function(k - 2L)) ||
      !.np_bw_selector_call_identical(sys.call(k - 2L), sys.call(k - 3L)) ||
      parents[[k - 2L]] != k - 3L ||
      !identical(sys.function(k - 3L), base::eval) ||
      !identical(sys.function(k - 4L), .np_eval_bw_call)) return(NULL)
  holder <- frames[[k - 4L]]
  value.env <- get0("value.env", holder, inherits = FALSE)
  native <- get0("native", holder, inherits = FALSE)
  if (!is.environment(value.env) || !identical(caller, value.env) ||
      !is.list(native) || !length(native[["values", exact = TRUE]])) return(NULL)
  expressions <- native[["expressions", exact = TRUE]]
  keys <- names(expressions)
  if (!is.list(expressions) || !length(keys) || anyDuplicated(keys) ||
      !all(c("xdat", "ydat") %in% keys) ||
      !all(keys %in% c("xdat", "ydat", "bws")) ||
      !all(vapply(keys, function(name)
        identical(mc[[name]], as.name(paste0(".np_native_", name))),
        logical(1L)))) return(NULL)
  allowed <- vapply(keys, function(name) {
    expr <- expressions[[name]]
    if (!identical(producer, "npregbw")) return(TRUE)
    is.symbol(expr) || (identical(name, "bws") &&
      is.numeric(expr) && is.null(attributes(expr)))
  }, logical(1L))
  if (!all(allowed)) return(NULL)
  expressions
}

# Formula wrappers own these generated descriptive terms. Native arithmetic
# consumes their prepared columns, not the caller's evaluation environment.
.npRmpi_terms_owner_wire <- function(tt) {
  if (inherits(tt, "terms") && is.environment(attr(tt, ".Environment", exact = TRUE)))
    attr(tt, ".Environment") <- emptyenv()
  tt
}

.npRmpi_plreg_role_wire <- function(x) {
  tt <- attr(x, "terms", exact = TRUE)
  if (inherits(tt, "terms") && is.environment(attr(tt, ".Environment", exact = TRUE)))
    attr(x, "terms") <- .npRmpi_terms_owner_wire(tt)
  x
}

.npRmpi_plreg_role_restore <- function(x, original) {
  expected <- attr(.npRmpi_plreg_role_wire(original), "terms", exact = TRUE)
  if (!is.null(expected) && identical(attr(x, "terms", exact = TRUE), expected))
    attr(x, "terms") <- attr(original, "terms", exact = TRUE)
  x
}


# Unconditional formula publication can share columns already sent as dat.
# The descriptor is private to one invocation and carries only a closed header.
# Other metadata keeps the incumbent materialized route. This does not replace
# the existing publication certificate or change user NA/callback semantics.
.npRmpi_bw_frame_exact <- function(x, y) identical(x, y, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)
.npRmpi_bw_frame_plain <- function(x) typeof(x) == "list" && !isS4(x) &&
  identical(attr(x, "class", exact = TRUE), "data.frame")
# Environments that serialize as a reference rather than by value. The same test
# the bandwidth.v1 fingerprint already uses (baseline np.autodispatch.R:619-620).
.npRmpi_bw_frame_reference_environment <- function(e) is.environment(e) &&
  (identical(e, globalenv()) || identical(e, baseenv()) || identical(e, emptyenv()) || isNamespace(e))
.npRmpi_bw_frame_producers <- c("npudensbw", "npudistbw")
.npRmpi_bw_frame_required_attributes <- c("names", "row.names", "class", "terms")
.npRmpi_bw_frame_optional_attributes <- c("na.action", ".np.subset")
.npRmpi_bw_frame_fields <- c("role", "token", "n", "columns", "order", "rows", "subset")


# A work budget prevents traversing or projecting an unbounded user expression.
# The exact projected-byte limit is checked separately after this plain grammar.
.npRmpi_bw_frame_syntax_cap <- function(p) 4096 + 2048 * p
.npRmpi_bw_frame_language <- function(x, budget = NULL, projected = TRUE, allow.formula = TRUE) {
  spend <- function(k) is.null(budget) || isTRUE((budget$left <- budget$left - k) >= 0)
  if (is.symbol(x)) return(spend(12 + nchar(as.character(x), type = "bytes")))
  if (is.null(x)) return(spend(4))
  if (is.atomic(x)) {
    if (!is.null(attributes(x)) || length(x) > 1L) return(FALSE)
    return(spend(if (is.character(x) && length(x))
      16 + nchar(x, type = "bytes", keepNA = FALSE) else 24))
  }
  if (!is.call(x) && !is.pairlist(x)) return(FALSE)
  cls <- attr(x, "class", exact = TRUE)
  if (!is.null(cls)) {
    if (!allow.formula || !is.call(x) || !identical(cls, "formula")) return(FALSE)
    # Do not dispatch length/[[ methods or copy an arbitrarily wide fake formula.
    wide <- tryCatch({ .subset2(x, 4L); TRUE },
                     subscriptOutOfBoundsError = function(e) FALSE)
    if (wide) return(FALSE)
    x <- unclass(x)
  }
  n <- length(x)
  if (!spend(8 + 4 * n)) return(FALSE)
  a <- attributes(x)
  if (!all(names(a) %in% c("names", "class", ".Environment"))) return(FALSE)
  env <- a[[".Environment"]]
  if (!is.null(env) && !(if (projected)
      .npRmpi_bw_frame_reference_environment(env) else is.environment(env))) return(FALSE)
  if ((!is.null(cls) || !is.null(env)) && !spend(96)) return(FALSE)
  tags <- names(x)
  if (!is.null(tags) && !spend(sum(12 + nchar(tags[nzchar(tags)],
      type = "bytes", keepNA = FALSE)))) return(FALSE)
  parts <- as.list(x)
  for (i in seq_len(n)) {
    if (.npRmpi_is_missing_call_arg(parts[[i]])) {
      if (!spend(4)) return(FALSE)
    } else if (!.npRmpi_bw_frame_language(parts[[i]], budget, projected, allow.formula)) return(FALSE)
  }
  TRUE
}
.npRmpi_bw_frame_preflight <- function(plan, frame, budget) {
  subset <- attr(frame, ".np.subset", exact = TRUE)
  .npRmpi_bw_frame_language(plan[["call", exact = TRUE]], budget, FALSE) &&
    .npRmpi_bw_frame_language(plan[["formula", exact = TRUE]], budget, FALSE) &&
    (!is.list(subset) || .npRmpi_bw_frame_language(subset[["expression"]], budget, FALSE, allow.formula = FALSE))
}
.npRmpi_bw_frame_bare <- function(v, type, n = NULL, attrs = character())
  typeof(v) == type && all(names(attributes(v)) %in% attrs) && (is.null(n) || length(v) == n)
.npRmpi_bw_frame_terms_attributes <- c("variables", "factors", "term.labels", "order", "intercept", "response",
                    "class", ".Environment", "predvars", "dataClasses", "specials")
.npRmpi_bw_frame_terms_plain <- function(tt) {
  a <- attributes(tt)
  if (!is.call(tt) || !identical(a[["class"]], c("terms", "formula")) ||
      anyDuplicated(names(a)) || !all(names(a) %in% .npRmpi_bw_frame_terms_attributes) ||
      !.npRmpi_bw_frame_reference_environment(a[[".Environment"]])) return(FALSE)
  body <- tt; attributes(body) <- NULL
  vars <- a[["variables"]]; labels <- a[["term.labels"]]
  if (!.npRmpi_bw_frame_language(body) || !is.call(vars) || !.npRmpi_bw_frame_language(vars) || !.npRmpi_bw_frame_bare(labels, "character")) return(FALSE)
  v <- length(vars) - 1L; t <- length(labels)
  pv <- a[["predvars"]]; fac <- a[["factors"]]; dc <- a[["dataClasses"]]; sp <- a[["specials"]]
  if (!(is.null(pv) || (is.call(pv) && .npRmpi_bw_frame_language(pv) && length(pv) == length(vars)))) return(FALSE)
  if (!.npRmpi_bw_frame_bare(a[["order"]], "integer", t) || !.npRmpi_bw_frame_bare(a[["intercept"]], "integer", 1L) ||
      !.npRmpi_bw_frame_bare(a[["response"]], "integer", 1L)) return(FALSE)
  if (!(.npRmpi_bw_frame_bare(fac, "integer", 0L) ||
        (.npRmpi_bw_frame_bare(fac, "integer", v * t, c("dim", "dimnames")) && identical(dim(fac), c(v, t)) &&
         is.list(dimnames(fac)) && is.null(attributes(dimnames(fac))) &&
         .npRmpi_bw_frame_bare(dimnames(fac)[[1L]], "character", v) && .npRmpi_bw_frame_bare(dimnames(fac)[[2L]], "character", t))))
    return(FALSE)
  if (!(is.null(dc) || (.npRmpi_bw_frame_bare(dc, "character", attrs = "names") && length(dc) <= v + 2L &&
                        .npRmpi_bw_frame_bare(names(dc), "character")))) return(FALSE)
  is.null(sp) || ((is.pairlist(sp) || (is.list(sp) && !is.object(sp))) && length(sp) <= v &&
    all(vapply(as.list(sp), function(s) is.null(s) || (.npRmpi_bw_frame_bare(s, "integer") && length(s) <= v), logical(1L))))
}
# The four policies stats provides. A user closure ships its body, formals and
# possibly its environment by value; this row does not need to carry user code.
.npRmpi_bw_frame_policy_plain <- function(f) is.null(f) || identical(f, na.omit) ||
  identical(f, na.exclude) || identical(f, na.pass) || identical(f, na.fail)
.npRmpi_bw_frame_omit_plain <- function(o) is.null(o) ||
  (.npRmpi_bw_frame_bare(o, "integer", attrs = c("names", "class")) &&
   (identical(attr(o, "class"), "omit") || identical(attr(o, "class"), "exclude")) &&
   (is.null(names(o)) || .npRmpi_bw_frame_bare(names(o), "character")))
.npRmpi_bw_frame_plan_fields <- c("producer", "call", "formula", "frame", "na.action", "terms", "xterms",
                    "chromoly", "ynames", "omit", "scalar.children")
# TRUE when the PROJECTED metadata of an admitted unconditional formula producer
# holds nothing but the reference and plain descriptive fields.
.npRmpi_bw_frame_plan_plain <- function(wire) {
  if (!is.list(wire) || is.object(wire) ||
      anyDuplicated(names(wire)) ||
      !all(names(wire) %in% c("names", "formula", "identity", "frame.reference")) ||
      !identical(wire[["names", exact = TRUE]], list()) ||
      (!is.null(wire[["identity", exact = TRUE]]) &&
       !identical(wire[["identity", exact = TRUE]], character()))) return(FALSE)
  p <- wire[["formula", exact = TRUE]]
  if (!is.list(p) || is.object(p) || anyDuplicated(names(p)) || !all(names(p) %in% .npRmpi_bw_frame_plan_fields)) return(FALSE)
  cl <- p[["call", exact = TRUE]]; fo <- p[["formula", exact = TRUE]]
  is.call(cl) && is.null(attr(cl, "class")) && .npRmpi_bw_frame_language(cl) &&
    is.call(fo) && identical(attr(fo, "class"), "formula") && .npRmpi_bw_frame_language(fo) &&
    .npRmpi_bw_frame_terms_plain(p[["terms", exact = TRUE]]) &&
    .npRmpi_bw_frame_policy_plain(p[["na.action", exact = TRUE]]) && .npRmpi_bw_frame_omit_plain(p[["omit", exact = TRUE]]) &&
    is.null(p[["xterms", exact = TRUE]]) && is.null(p[["chromoly", exact = TRUE]]) &&
    identical(p[["ynames", exact = TRUE]], character()) &&
    identical(p[["scalar.children", exact = TRUE]], FALSE)
}


.npRmpi_bw_frame_decode <- function(metadata, call, values) {
  ref <- metadata[["frame.reference", exact = TRUE]]
  if (is.null(ref)) return(metadata)
  bad <- function() stop("invalid prepared-frame reference", call. = FALSE)
  if (!is.list(ref) || is.object(ref) || !identical(names(ref), .npRmpi_bw_frame_fields) ||
      !identical(ref[["role"]], "dat") || !is.character(ref[["token"]]) ||
      length(ref[["token"]]) != 1L || is.na(ref[["token"]]) ||
      !is.null(attributes(ref[["token"]]))) bad()
  p <- metadata[["formula", exact = TRUE]]
  if (!is.list(p) || is.object(p) || !isTRUE(p[["producer", exact = TRUE]] %in% .npRmpi_bw_frame_producers) ||
      !"frame" %in% names(p) || !is.null(p[["frame", exact = TRUE]]) ||
      !is.symbol(call[["dat"]]) || !identical(as.character(call[["dat"]]), ref[["token"]]) ||
      !ref[["token"]] %in% names(values)) bad()
  d <- values[[ref[["token"]]]]
  if (!.npRmpi_bw_frame_plain(d)) bad()
  dn <- attr(d, "names", exact = TRUE)
  n <- ref[["n"]]
  if (!identical(sort(names(attributes(d))), c("class", "names", "row.names")) ||
      !is.integer(n) || length(n) != 1L || is.na(n) || !identical(.row_names_info(d, 2L), n) ||
      !.npRmpi_bw_frame_exact(dn, ref[["columns"]])) bad()
  order <- ref[["order"]]; subset <- ref[["subset"]]
  terms <- p[["terms", exact = TRUE]]; omit <- p[["omit", exact = TRUE]]
  if (!is.character(order) || !is.null(attributes(order)) || anyNA(order) || anyDuplicated(order) ||
      !all(.npRmpi_bw_frame_required_attributes %in% order) || !all(order %in% c(.npRmpi_bw_frame_required_attributes, .npRmpi_bw_frame_optional_attributes)) ||
      !inherits(terms, "terms") ||
      ("na.action" %in% order) != !is.null(omit) ||
      (".np.subset" %in% order) != !is.null(subset)) bad()
  rows <- ref[["rows"]]
  rn <- if (identical(rows, "source")) .row_names_info(d, 0L) else
    if (is.integer(rows) && is.null(attributes(rows)) && length(rows) == 2L &&
        is.na(rows[1L]) && !is.na(rows[2L]) && abs(rows[2L]) == n) rows else bad()
  if (!is.null(subset) &&
      (!is.list(subset) || is.object(subset) || !identical(names(subset), c("expression", "envir")) ||
       !is.language(subset[["expression"]]) || !.npRmpi_bw_frame_language(subset[["expression"]], allow.formula = FALSE) ||
       !.npRmpi_bw_frame_reference_environment(subset[["envir"]]))) bad()
  attrs <- setNames(vector("list", length(order)), order)
  for (k in order) attrs[k] <- list(switch(k, names = dn, row.names = rn, class = "data.frame",
    terms = terms, na.action = omit, .np.subset = subset))
  out <- unclass(d)                                   # shares every column with the source
  attributes(out) <- attrs
  if (!identical(.row_names_info(out, 2L), n)) bad()
  p["frame"] <- list(out)
  metadata["formula"] <- list(p)
  metadata[["frame.reference"]] <- NULL
  metadata
}

.npRmpi_bw_frame_encode <- function(metadata, call, values) {
  none <- list(status = "none", reason = "", metadata = metadata)
  decline <- function(why) list(status = "decline", reason = why, metadata = NULL)
  p <- if (is.list(metadata)) metadata[["formula", exact = TRUE]] else NULL
  if (!is.list(p) || is.object(p)) return(none)
  producer <- p[["producer", exact = TRUE]]
  if (!is.character(producer) || length(producer) != 1L || is.na(producer) ||
      !producer %in% .npRmpi_bw_frame_producers) return(none)
  f <- p[["frame", exact = TRUE]]
  if (is.null(f)) return(none)
  token <- call[["dat"]]
  if (!is.symbol(token) || !as.character(token) %in% names(values)) return(decline("no dat token"))
  token <- as.character(token)
  d <- values[[token]]
  if (!.npRmpi_bw_frame_plain(f) || !.npRmpi_bw_frame_plain(d)) return(decline("not a plain data frame"))
  fn <- attr(f, "names", exact = TRUE); dn <- attr(d, "names", exact = TRUE)
  uf <- unclass(f); ud <- unclass(d)
  if (!is.character(fn) || length(uf) != length(ud) || !.npRmpi_bw_frame_exact(fn, dn) ||
      !.npRmpi_transport_same_encoding(fn, dn)) return(decline("column names"))
  for (i in seq_along(uf)) {
    x <- uf[[i]]; y <- ud[[i]]
    if (!(.npRmpi_bw_frame_exact(x, y) && .npRmpi_transport_plain_value(x) &&
          .npRmpi_transport_same_encoding(x, y))) return(decline("column values or attributes"))
  }
  fa <- attributes(f); fa[["row.names"]] <- .row_names_info(f, 0L)
  da <- attributes(d); da[["row.names"]] <- .row_names_info(d, 0L)
  if (anyDuplicated(names(fa)) || !all(.npRmpi_bw_frame_required_attributes %in% names(fa)) ||
      !all(names(fa) %in% c(.npRmpi_bw_frame_required_attributes, .npRmpi_bw_frame_optional_attributes)))
    return(decline("frame attribute outside the admitted set"))
  if (!identical(sort(names(da)), c("class", "names", "row.names"))) return(decline("source attributes"))
  terms <- p[["terms", exact = TRUE]]; omit <- p[["omit", exact = TRUE]]
  if (!inherits(terms, "terms") || !.npRmpi_bw_frame_exact(fa[["terms"]], terms)) return(decline("terms are not the plan's"))
  if ("na.action" %in% names(fa)) {
    if (is.null(omit) || !.npRmpi_bw_frame_exact(fa[["na.action"]], omit)) return(decline("omission record is not the plan's"))
  } else if (!is.null(omit)) return(decline("plan omission without a frame record"))
  rows <- if (.npRmpi_bw_frame_exact(fa[["row.names"]], da[["row.names"]]) &&
              .npRmpi_transport_same_encoding(fa[["row.names"]], da[["row.names"]])) "source" else {
    r <- fa[["row.names"]]
    if (!(is.integer(r) && length(r) == 2L && is.na(r[1L]) && !is.na(r[2L])))
      return(decline("row names differ from the source"))
    r
  }
  subset <- NULL
  if (".np.subset" %in% names(fa)) {
    subset <- fa[[".np.subset"]]
    if (!is.list(subset) || is.object(subset) || !identical(names(subset), c("expression", "envir")) ||
        !is.language(subset[["expression"]]) || !.npRmpi_bw_frame_language(subset[["expression"]], allow.formula = FALSE) ||
        !.npRmpi_bw_frame_reference_environment(subset[["envir"]]))
      return(decline("subset metadata would travel by value"))
  }
  p["frame"] <- list(NULL)
  metadata["formula"] <- list(p)
  metadata[["frame.reference"]] <- list(role = "dat", token = token, n = .row_names_info(d, 2L),
    columns = dn, order = names(fa), rows = rows, subset = subset)
  back <- tryCatch(.npRmpi_bw_frame_decode(metadata, call, values), error = function(e) NULL)
  g <- if (is.null(back)) NULL else back[["formula", exact = TRUE]][["frame", exact = TRUE]]
  if (is.null(g) || !.npRmpi_bw_frame_exact(g, f) || !.npRmpi_transport_same_encoding(g, f) ||
      !.npRmpi_bw_frame_exact(.row_names_info(g, 0L), .row_names_info(f, 0L)))
    return(decline("round trip is not exact"))
  list(status = "admit", reason = "", metadata = metadata)
}
# The whole admission decision for the two new rows, in the order a dispatcher
# would run it. "none" leaves the caller to project as it does today.
.npRmpi_bw_frame_admit <- function(metadata, call, values) {
  decline <- function(why) list(status = "decline", reason = why, metadata = NULL)
  tryCatch({
    plan <- if (is.list(metadata)) metadata[["formula", exact = TRUE]] else NULL
    frame <- if (is.list(plan) && !is.object(plan)) plan[["frame", exact = TRUE]] else NULL
    producer <- if (is.null(frame)) NULL else plan[["producer", exact = TRUE]]
    cap <- NULL
    if (typeof(frame) == "list" && is.character(producer) && length(producer) == 1L &&
        !is.na(producer) && producer %in% .npRmpi_bw_frame_producers) {
      # Older supported R versions lack the precise inspection error classes.
      # Keep their incumbent formula transport; native admission is unaffected.
      if (getRversion() < "4.2.0")
        return(decline("prepared-frame reuse requires R 4.2.0 condition classes"))
      cap <- .npRmpi_bw_frame_syntax_cap(length(attr(frame, "names", exact = TRUE)))
      budget <- new.env(parent = emptyenv())
      budget$left <- cap
      if (!.npRmpi_bw_frame_preflight(plan, frame, budget))
        return(decline(if (isTRUE(budget$left < 0)) "call syntax exceeds the lease's metadata bound"
                       else "call syntax is outside the plain grammar"))
    }
    r <- .npRmpi_bw_frame_encode(metadata, call, values)
    if (!identical(r$status, "admit")) return(r)
    omit <- r$metadata[["formula", exact = TRUE]][["omit", exact = TRUE]]
    if (!is.null(omit) && length(unclass(omit)) > 0L)
      return(decline("nonempty omission record exceeds new transport scope"))
    wire <- .npRmpi_bw_metadata_wire(r$metadata, plain.language = TRUE)
    if (!.npRmpi_bw_frame_plan_plain(wire))
      return(decline("plan holds a value outside the plain grammar"))
    # Only validated, projected syntax reaches the serializer: no caller frames
    # or closures. Use this R's transport format, including its list overhead.
    syntax <- list(call = wire$formula$call, formula = wire$formula$formula,
                   subset = wire$frame.reference$subset[["expression"]])
    if (length(serialize(syntax, NULL)) > cap)
      return(decline("projected call syntax exceeds the lease's metadata bound"))
    list(status = "admit", reason = "", metadata = wire)
  }, stackOverflowError = function(e) decline("call syntax is too deeply nested to inspect"))
}
