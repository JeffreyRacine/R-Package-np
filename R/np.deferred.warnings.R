# Only explicitly wired package warnings may enter this coordinator-local scope.
# No state is retained by the namespace or sent in an MPI message.
.npRmpi_warning_scope <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (identical(sys.function(i), .npRmpi_with_deferred_package_warnings))
      return(get("pending", envir = frames[[i]], inherits = FALSE))
  }
  NULL
}

.npRmpi_warning_scope_activate <- function() {
  pending <- .npRmpi_warning_scope()
  if (!is.null(pending)) pending$active <- TRUE
  invisible(NULL)
}

.npRmpi_defer_named_warning <- function(code) {
  # warn=2 keeps its established error/coordination path. This policy does not
  # repair rank-asymmetric warning-to-error conversion or manual SPMD handlers.
  if (getOption("warn", 0L) >= 2L) return(force(code))
  pending <- .npRmpi_warning_scope()
  if (is.null(pending) || !isTRUE(pending$active)) return(force(code))
  withCallingHandlers(force(code), warning = function(w) {
    pending$items[[length(pending$items) + 1L]] <- w
    invokeRestart("muffleWarning")
  })
}

.npRmpi_with_deferred_package_warnings <- function(code) {
  pending <- new.env(parent = emptyenv())
  pending$active <- FALSE
  pending$items <- list()
  value <- withVisible(tryCatch(force(code), error = function(e) {
    pending$active <- FALSE
    if (length(pending$items)) {
      attr(e, ".npRmpi.deferred.warnings") <- c(
        pending$items,
        attr(e, ".npRmpi.deferred.warnings", exact = TRUE))
    }
    stop(e)
  }))
  # Distributed execution and its exit guards, plus metadata finalization,
  # have returned. Exiting/reentrant user handlers can no longer strand peers.
  pending$active <- FALSE
  for (item in pending$items)
    warning(item)
  if (value$visible) value$value else invisible(value$value)
}
