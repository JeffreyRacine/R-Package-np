.np_nomad_interrupt_condition <- function(message = "NOMAD search interrupted by user") {
  structure(
    list(message = as.character(message)[1L], call = NULL),
    class = c("interrupt", "condition")
  )
}

.np_nomad_native_status <- function(native,
                                    route,
                                    default.message = "NOMAD rejected the supplied parameters") {
  if (!is.list(native))
    stop(sprintf("%s returned a malformed result", route), call. = FALSE)

  status.raw <- if (!is.null(native$native_status)) {
    native$native_status[1L]
  } else if (!is.null(native$status)) {
    native$status[1L]
  } else {
    0L
  }
  status.integer <- suppressWarnings(as.integer(status.raw))
  status.failed <- if (!is.na(status.integer)) {
    !identical(status.integer, 0L)
  } else {
    !(tolower(as.character(status.raw)) %in% c("", "ok", "success"))
  }

  result.status <- if (!is.null(native$result_status)) {
    suppressWarnings(as.integer(native$result_status[1L]))
  } else {
    0L
  }
  result.failed <- is.na(result.status) || !identical(result.status, 0L)
  if (!status.failed && !result.failed)
    return(invisible(native))

  message <- if (!is.null(native$message) &&
                 nzchar(as.character(native$message[1L]))) {
    as.character(native$message[1L])
  } else {
    default.message
  }

  if (identical(status.integer, 4L) || identical(result.status, 4L))
    stop(.np_nomad_interrupt_condition(message))

  stop(
    sprintf(
      "%s failed (status=%s, result_status=%s): %s",
      route,
      as.character(status.raw),
      result.status,
      message
    ),
    call. = FALSE
  )
}

.npRmpi_mpi_interrupt_state <- function(action = 0L) {
  isTRUE(.Call(
    "C_np_mpi_interrupt_state",
    as.integer(action)[1L],
    PACKAGE = "npRmpi"
  ))
}

.npRmpi_mpi_interrupt_scope <- function(action) {
  isTRUE(.Call(
    "C_np_mpi_interrupt_scope",
    as.integer(action)[1L],
    PACKAGE = "npRmpi"
  ))
}

.npRmpi_mpi_interrupt_set_pending <- function() {
  .npRmpi_mpi_interrupt_state(2L)
  invisible(TRUE)
}

.npRmpi_mpi_interrupt_defer_on_master <- function() {
  if (!.npRmpi_mpi_interrupt_state(3L))
    return(FALSE)
  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) NA_integer_)
  !is.na(rank) && identical(rank, 0L)
}

.npRmpi_signal_pending_interrupt <- function(
    message = "MPI-backed computation interrupted by user") {
  if (.npRmpi_mpi_interrupt_state(4L))
    stop(.np_nomad_interrupt_condition(message))
  invisible(FALSE)
}

.npRmpi_with_command_interrupt_scope <- function(code) {
  scope.started <- .npRmpi_mpi_interrupt_scope(1L)
  if (!scope.started) {
    value <- force(code)
    .npRmpi_signal_pending_interrupt()
    return(value)
  }

  scope.closed <- FALSE
  on.exit({
    if (!scope.closed)
      .npRmpi_mpi_interrupt_scope(2L)
  }, add = TRUE)

  caught <- NULL
  value <- tryCatch(
    force(code),
    error = function(e) {
      caught <<- e
      NULL
    },
    interrupt = function(e) {
      caught <<- e
      NULL
    }
  )
  scope.interrupted <- .npRmpi_mpi_interrupt_scope(2L)
  scope.closed <- TRUE
  if (scope.interrupted)
    stop(.np_nomad_interrupt_condition(
      "MPI-backed computation interrupted by user"
    ))
  if (!is.null(caught))
    stop(caught)
  value
}
