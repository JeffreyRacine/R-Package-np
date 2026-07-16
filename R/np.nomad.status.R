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
