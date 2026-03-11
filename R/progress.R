.np_progress_is_interactive <- function() {
  interactive()
}

.np_progress_now <- function() {
  as.numeric(proc.time()[["elapsed"]])
}

.np_progress_enabled <- function(domain = "general") {
  isTRUE(getOption("np.messages", TRUE)) && isTRUE(.np_progress_is_interactive())
}

.np_io_pkg_prefix <- function() {
  "[np]"
}

.np_io_prefix_text <- function(text) {
  prefix <- .np_io_pkg_prefix()
  if (!is.character(text)) {
    return(text)
  }

  vapply(
    text,
    function(line) {
      if (is.na(line) || startsWith(line, prefix)) {
        return(line)
      }

      paste(prefix, line)
    },
    character(1L)
  )
}

.np_message <- function(...) {
  message(.np_io_prefix_text(paste(..., sep = "", collapse = "")))
  invisible(NULL)
}

.np_warning <- function(..., call. = TRUE, immediate. = FALSE, noBreaks. = FALSE) {
  warning(
    .np_io_prefix_text(paste(..., sep = "", collapse = "")),
    call. = call.,
    immediate. = immediate.,
    noBreaks. = noBreaks.
  )
}

.np_progress_pkg_prefix <- function() {
  .np_io_pkg_prefix()
}

.np_progress_emit <- function(line) {
  .np_message(line)
  invisible(NULL)
}

.np_progress_fmt_num <- function(x) {
  formatC(x, digits = 1L, format = "f")
}

.np_progress_format_known_total <- function(state, done, detail = NULL, now = .np_progress_now()) {
  total <- state$total
  pct <- if (isTRUE(total > 0)) 100 * done / total else 0
  elapsed <- max(0, now - state$started)
  eta <- if (isTRUE(done > 0) && isTRUE(total >= done)) elapsed * (total - done) / done else 0

  line <- sprintf(
    "%s %s %s/%s (%s%%, elapsed %ss, eta %ss)",
    state$pkg_prefix,
    state$label,
    format(done),
    format(total),
    .np_progress_fmt_num(pct),
    .np_progress_fmt_num(elapsed),
    .np_progress_fmt_num(eta)
  )

  if (!is.null(detail)) {
    line <- paste0(line, ": ", detail)
  }

  line
}

.np_progress_format_unknown_total <- function(state, done = NULL, detail = NULL, now = .np_progress_now()) {
  elapsed <- max(0, now - state$started)
  line <- sprintf("%s %s... ", state$pkg_prefix, state$label)

  if (!is.null(done)) {
    line <- paste0(line, "iteration ", format(done), ", ")
  }

  line <- paste0(line, "elapsed ", .np_progress_fmt_num(elapsed), "s")

  if (!is.null(detail)) {
    line <- paste0(line, ": ", detail)
  }

  line
}

.np_progress_begin <- function(label, total = NULL, domain = "general") {
  known_total <- !is.null(total)
  throttle_sec <- if (known_total) 0.5 else 2.0
  started <- .np_progress_now()

  list(
    enabled = .np_progress_enabled(domain = domain),
    pkg_prefix = .np_progress_pkg_prefix(),
    label = label,
    total = total,
    known_total = known_total,
    started = started,
    last_emit = started - throttle_sec,
    throttle_sec = throttle_sec,
    last_done = if (known_total) 0 else NULL,
    domain = domain
  )
}

.np_progress_step <- function(state, done = NULL, detail = NULL) {
  if (!is.null(done)) {
    state$last_done <- done
  }

  if (!isTRUE(state$enabled)) {
    return(state)
  }

  now <- .np_progress_now()
  if ((now - state$last_emit) < state$throttle_sec) {
    return(state)
  }

  line <- if (isTRUE(state$known_total)) {
    .np_progress_format_known_total(
      state = state,
      done = if (is.null(done)) state$last_done else done,
      detail = detail,
      now = now
    )
  } else {
    .np_progress_format_unknown_total(
      state = state,
      done = if (is.null(done)) state$last_done else done,
      detail = detail,
      now = now
    )
  }

  .np_progress_emit(line)
  state$last_emit <- now
  state
}

.np_progress_end <- function(state, detail = NULL) {
  if (!isTRUE(state$enabled)) {
    return(invisible(state))
  }

  now <- .np_progress_now()

  if (isTRUE(state$known_total)) {
    done <- if (is.null(state$total)) state$last_done else state$total
    line <- .np_progress_format_known_total(state = state, done = done, detail = detail, now = now)
    if (!identical(done, state$last_done) || (now - state$last_emit) > 0) {
      .np_progress_emit(line)
      state$last_emit <- now
      state$last_done <- done
    }
    return(invisible(state))
  }

  if (!is.null(state$last_done)) {
    line <- .np_progress_format_unknown_total(state = state, done = state$last_done, detail = detail, now = now)
    if ((now - state$last_emit) > 0) {
      .np_progress_emit(line)
      state$last_emit <- now
    }
  }

  invisible(state)
}

.np_progress_note <- function(label) {
  if (.np_progress_enabled()) {
    .np_progress_emit(sprintf("%s %s", .np_progress_pkg_prefix(), label))
  }

  invisible(NULL)
}

.np_progress_select_bandwidth <- function(label, expr) {
  .np_progress_note(label)
  .np_progress_with_legacy_suppressed(expr)
}

.np_progress_with_legacy_suppressed <- function(expr) {
  old_np_messages <- getOption("np.messages")
  on.exit(options(np.messages = old_np_messages), add = TRUE)
  options(np.messages = FALSE)
  force(expr)
}
