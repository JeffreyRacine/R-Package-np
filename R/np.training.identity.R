# Hash only explicitly retained training values. All other fields continue
# through the complete metadata fingerprint and its existing ownership rules.
.npRmpi_training_value_eligible <- function(x) {
  if (isS4(x) || !typeof(x) %in%
      c("NULL", "logical", "integer", "double", "complex", "character", "raw", "list"))
    return(FALSE)
  cl <- attr(x, "class", exact = TRUE)
  classes <- list(NULL, "data.frame", "factor", c("ordered", "factor"),
                  "Date", c("POSIXct", "POSIXt"), "difftime", "omit", "exclude")
  if (!any(vapply(classes, identical, logical(1L), y = cl)))
    return(FALSE)
  view <- if (is.object(x)) unclass(x) else x
  if (is.list(view)) for (i in seq_along(view))
    if (!.npRmpi_training_value_eligible(.subset2(view, i))) return(FALSE)
  attrs <- attributes(x)
  for (i in seq_along(attrs)) {
    a <- .subset2(attrs, i)
    if (identical(names(attrs)[i], "terms") && inherits(a, "terms")) next
    if (!.npRmpi_training_value_eligible(a)) return(FALSE)
  }
  TRUE
}

.npRmpi_training_value_certificate <- function(x) {
  if (is.null(x) || !.npRmpi_training_value_eligible(x)) return(NULL)
  path <- tempfile("npRmpi-training-value-", fileext = ".bin")
  con <- NULL
  on.exit({
    if (!is.null(con)) try(close(con), silent = TRUE)
    unlink(path)
  }, add = TRUE)
  con <- file(path, "wb")
  put <- function(value) serialize(value, con, xdr = FALSE, version = 2L)
  emit <- function(value) {
    view <- if (is.object(value)) unclass(value) else value
    n <- length(view)
    put(list(type = typeof(value), length = n))
    attrs <- attributes(value)
    put(names(attrs))
    for (i in seq_along(attrs)) {
      a <- .subset2(attrs, i)
      if (identical(names(attrs)[i], "terms") && inherits(a, "terms")) {
        digest <- .npRmpi_autodispatch_fingerprint(a)
        if (is.na(digest)) stop("training terms fingerprint failed", call. = FALSE)
        put(list(protected.metadata = digest))
      } else emit(a)
    }
    if (is.list(view)) {
      for (i in seq_along(view)) emit(.subset2(view, i))
    } else if (n) for (start in seq.int(1, n, by = 32768L)) {
      piece <- .subset(view, seq.int(start, min(n, start + 32767L)))
      attributes(piece) <- NULL
      put(piece)
    }
    invisible(NULL)
  }
  domain <- "npRmpi.training.value.v1.chunk32768"
  put(domain)
  emit(x)
  close(con)
  con <- NULL
  digest <- unname(as.character(tools::md5sum(path)))
  bytes <- unname(file.info(path)$size)
  if (length(digest) != 1L || is.na(digest) || !nzchar(digest) ||
      length(bytes) != 1L || !is.finite(bytes))
    stop("training value fingerprint failed", call. = FALSE)
  list(version = domain, bytes = bytes, md5 = digest)
}

.npRmpi_training_value_projection <- function(x) {
  containers <- c("rbandwidth", "lsqregressionbandwidth", "plbandwidth",
                  "sibandwidth", "scbandwidth", "conbandwidth", "condbandwidth",
                  "bandwidth", "dbandwidth")
  roles <- c("xdat", "ydat", "zdat", "dat")
  call.roles <- c(roles, "txdat", "tydat", "tzdat", "tdat", "data")
  project <- function(value) {
    cl <- attr(value, "class", exact = TRUE)
    if (!identical(typeof(value), "list") || isS4(value) ||
        !(is.null(cl) || any(vapply(as.list(containers), identical,
                                  logical(1L), y = cl))))
      return(list(value = value, changed = FALSE))
    view <- if (is.object(value)) unclass(value) else value
    nms <- attr(view, "names", exact = TRUE)
    out <- vector("list", length(view))
    changed <- FALSE
    for (i in seq_along(view)) {
      z <- .subset2(view, i)
      key <- if (is.null(nms)) "" else nms[[i]]
      if (key %in% c(".np.native.training", ".np.formula.training") &&
          identical(typeof(z), "list") && !is.object(z) && !isS4(z)) {
        zn <- attr(z, "names", exact = TRUE)
        fields <- if (identical(key, ".np.native.training")) roles else "frame"
        replacement <- vector("list", length(z))
        local.change <- FALSE
        for (j in seq_along(z)) {
          a <- .subset2(z, j)
          if (!is.null(zn) && zn[[j]] %in% fields && !is.null(a)) {
            h <- .npRmpi_training_value_certificate(a)
            if (!is.null(h)) {
              a <- list(training.role = zn[[j]], certificate = h)
              local.change <- TRUE
            }
          }
          replacement[j] <- list(a)
        }
        if (local.change) {
          attributes(replacement) <- attributes(z)
          z <- replacement
          changed <- TRUE
        }
      } else if (identical(key, "call") && is.call(z)) {
        # A call is projected only in this named field, never at the root.
        parts <- as.list.default(z)
        zn <- attr(parts, "names", exact = TRUE)
        local.change <- FALSE
        for (j in which(zn %in% call.roles)) {
          if (identical(.subset2(parts, j), quote(expr = ))) next
          a <- .subset2(parts, j)
          if (is.null(a) || is.language(a)) next
          h <- .npRmpi_training_value_certificate(a)
          if (!is.null(h)) {
            parts[j] <- list(list(call.role = zn[[j]], certificate = h))
            local.change <- TRUE
          }
        }
        if (local.change) {
          replacement <- as.call(parts)
          attributes(replacement) <- attributes(z)
          z <- replacement
          changed <- TRUE
        }
      } else {
        child <- project(z)
        if (child$changed) {
          z <- child$value
          changed <- TRUE
        }
      }
      out[i] <- list(z)
    }
    if (!changed) return(list(value = value, changed = FALSE))
    attributes(out) <- attributes(value)
    list(value = out, changed = TRUE)
  }
  project(x)
}
