.np_progress_is_interactive <- function() {
  interactive()
}

.np_progress_now <- function() {
  as.numeric(proc.time()[["elapsed"]])
}

.np_progress_enabled <- function(domain = "general") {
  (isTRUE(.np_progress_runtime$force_enabled) || isTRUE(getOption("np.messages", TRUE))) &&
    isTRUE(.np_progress_is_interactive())
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

.np_progress_runtime <- local({
  env <- new.env(parent = emptyenv())
  env$bandwidth_depth <- 0L
  env$bandwidth_label <- NULL
  env$bandwidth_state <- NULL
  env$bandwidth_old_messages <- NULL
  env$force_enabled <- FALSE
  env
})

.np_progress_start_grace_sec <- function(known_total = FALSE, domain = "general") {
  default <- if (identical(domain, "plot")) {
    0.75
  } else if (isTRUE(known_total)) {
    0.75
  } else {
    1.0
  }

  option_name <- if (identical(domain, "plot")) {
    "np.plot.progress.start.grace.sec"
  } else if (isTRUE(known_total)) {
    "np.progress.start.grace.known.sec"
  } else {
    "np.progress.start.grace.unknown.sec"
  }

  val <- suppressWarnings(as.numeric(getOption(option_name, default))[1L])
  if (!is.finite(val) || is.na(val) || val < 0) {
    val <- default
  }

  val
}

.np_progress_format_known_total <- function(state, done, detail = NULL, now = .np_progress_now()) {
  total <- state$total
  pct <- if (isTRUE(total > 0)) 100 * done / total else 0
  elapsed <- max(0, now - state$started)
  eta <- if (isTRUE(done > 0) && isTRUE(total >= done)) elapsed * (total - done) / done else 0

  progress_token <- if (isTRUE(state$bandwidth_multistart)) {
    sprintf("multistart %s/%s", format(done), format(total))
  } else {
    sprintf("%s/%s", format(done), format(total))
  }

  line <- sprintf(
    "%s %s %s (%s%%, elapsed %ss, eta %ss)",
    state$pkg_prefix,
    state$label,
    progress_token,
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
    domain = domain,
    start_note = paste0(label, "..."),
    start_note_pending = TRUE,
    start_note_grace_sec = .np_progress_start_grace_sec(known_total = known_total, domain = domain),
    last_line = NULL,
    last_emitted_done = NULL,
    last_emitted_detail = NULL
  )
}

.np_progress_maybe_emit_start_note <- function(state, now = .np_progress_now()) {
  if (!isTRUE(state$enabled) || !isTRUE(state$start_note_pending)) {
    return(state)
  }

  if ((now - state$started) < state$start_note_grace_sec) {
    return(state)
  }

  .np_progress_emit(state$start_note)
  state$start_note_pending <- FALSE
  state$last_line <- state$start_note

  if (isTRUE(state$start_note_consumes_throttle)) {
    state$last_emit <- now
  }

  state
}

.np_progress_step <- function(state, done = NULL, detail = NULL) {
  if (!is.null(done)) {
    state$last_done <- done
  }

  if (!isTRUE(state$enabled)) {
    return(state)
  }

  now <- .np_progress_now()
  state <- .np_progress_maybe_emit_start_note(state = state, now = now)
  if (isTRUE(state$start_note_pending)) {
    return(state)
  }

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

  if (!identical(line, state$last_line)) {
    .np_progress_emit(line)
    state$last_emit <- now
    state$last_line <- line
    state$last_emitted_done <- if (is.null(done)) state$last_done else done
    state$last_emitted_detail <- detail
  }

  state
}

.np_progress_end <- function(state, detail = NULL) {
  if (!isTRUE(state$enabled)) {
    return(invisible(state))
  }

  now <- .np_progress_now()
  if (is.null(state$last_line) &&
      (now - state$started) < state$start_note_grace_sec) {
    return(invisible(state))
  }

  if (isTRUE(state$known_total)) {
    done <- if (is.null(state$total)) state$last_done else state$total
    line <- .np_progress_format_known_total(state = state, done = done, detail = detail, now = now)
    if (!(identical(done, state$last_emitted_done) && identical(detail, state$last_emitted_detail))) {
      .np_progress_emit(line)
      state$last_emit <- now
      state$last_done <- done
      state$last_line <- line
      state$last_emitted_done <- done
      state$last_emitted_detail <- detail
    }
    return(invisible(state))
  }

  if (!is.null(state$last_done)) {
    line <- .np_progress_format_unknown_total(state = state, done = state$last_done, detail = detail, now = now)
    if (!(identical(state$last_done, state$last_emitted_done) && identical(detail, state$last_emitted_detail))) {
      .np_progress_emit(line)
      state$last_emit <- now
      state$last_line <- line
      state$last_emitted_done <- state$last_done
      state$last_emitted_detail <- detail
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

.np_progress_bandwidth_active <- function() {
  isTRUE(.np_progress_runtime$bandwidth_depth > 0L)
}

.np_progress_bandwidth_detail <- function(done, total) {
  NULL
}

.np_progress_bandwidth_set_total <- function(total) {
  state <- .np_progress_runtime$bandwidth_state
  total <- suppressWarnings(as.integer(total)[1L])

  if (is.null(state) || is.na(total) || total <= 1L) {
    return(invisible(NULL))
  }

  if (isTRUE(state$known_total) && identical(as.integer(state$total), total)) {
    return(invisible(NULL))
  }

  label <- .np_progress_runtime$bandwidth_label
  if (is.null(label))
    label <- state$label

  upgraded <- .np_progress_begin(label = label, total = total, domain = "general")
  upgraded$started <- state$started
  upgraded$last_emit <- state$started - upgraded$throttle_sec
  upgraded$start_note_pending <- state$start_note_pending
  upgraded$last_line <- state$last_line
  upgraded$bandwidth_multistart <- TRUE
  if (!isTRUE(state$start_note_pending) && identical(state$last_line, state$start_note)) {
    upgraded$start_note_pending <- FALSE
  }

  .np_progress_runtime$bandwidth_state <- upgraded
  invisible(NULL)
}

.np_progress_bandwidth_multistart_step <- function(done, total = NULL) {
  state <- .np_progress_runtime$bandwidth_state

  if (is.null(state)) {
    return(invisible(NULL))
  }

  if (!is.null(total)) {
    .np_progress_bandwidth_set_total(total = total)
    state <- .np_progress_runtime$bandwidth_state
  }

  if (is.null(state) || !isTRUE(state$known_total)) {
    return(invisible(NULL))
  }

  done <- suppressWarnings(as.integer(done)[1L])
  if (is.na(done) || done < 1L) {
    return(invisible(NULL))
  }

  total <- as.integer(state$total)
  done <- min(done, total)
  .np_progress_runtime$bandwidth_state <- .np_progress_step(
    state = state,
    done = done,
    detail = .np_progress_bandwidth_detail(done = done, total = total)
  )

  invisible(NULL)
}

.np_progress_bandwidth_finish <- function() {
  state <- .np_progress_runtime$bandwidth_state

  if (is.null(state)) {
    return(invisible(NULL))
  }

  if (isTRUE(state$known_total) && isTRUE(state$total > 1L)) {
    total <- as.integer(state$total)
    .np_progress_runtime$bandwidth_state <- .np_progress_end(
      state,
      detail = .np_progress_bandwidth_detail(done = total, total = total)
    )
  } else {
    .np_progress_runtime$bandwidth_state <- .np_progress_end(state)
  }

  invisible(NULL)
}

.np_progress_bandwidth_clear <- function() {
  .np_progress_runtime$bandwidth_depth <- 0L
  .np_progress_runtime$bandwidth_label <- NULL
  .np_progress_runtime$bandwidth_state <- NULL
  .np_progress_runtime$bandwidth_old_messages <- NULL
  .np_progress_runtime$force_enabled <- FALSE
  invisible(NULL)
}

.np_progress_select_bandwidth <- function(label, expr) {
  starting <- !.np_progress_bandwidth_active()
  if (starting) {
    .np_progress_runtime$bandwidth_old_messages <- getOption("np.messages", TRUE)
    .np_progress_runtime$force_enabled <- isTRUE(.np_progress_runtime$bandwidth_old_messages)
    options(np.messages = FALSE)
    .np_progress_runtime$bandwidth_label <- as.character(label)[1L]
    .np_progress_runtime$bandwidth_state <- .np_progress_begin(
      label = .np_progress_runtime$bandwidth_label,
      domain = "general"
    )
  }

  .np_progress_runtime$bandwidth_depth <- as.integer(.np_progress_runtime$bandwidth_depth) + 1L
  on.exit({
    .np_progress_runtime$bandwidth_depth <- max(
      0L,
      as.integer(.np_progress_runtime$bandwidth_depth) - 1L
    )
    if (starting) {
      .np_progress_bandwidth_finish()
      options(np.messages = .np_progress_runtime$bandwidth_old_messages)
      .np_progress_bandwidth_clear()
    }
  }, add = TRUE)

  force(expr)
}

.np_progress_with_legacy_suppressed <- function(expr) {
  old_np_messages <- getOption("np.messages")
  on.exit(options(np.messages = old_np_messages), add = TRUE)
  options(np.messages = FALSE)
  force(expr)
}
