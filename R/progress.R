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

.np_progress_call_name <- function(call) {
  head <- call[[1L]]

  if (is.symbol(head)) {
    return(as.character(head))
  }

  if (is.call(head) &&
      identical(head[[1L]], as.name("::")) &&
      length(head) >= 3L &&
      is.symbol(head[[2L]]) &&
      is.symbol(head[[3L]])) {
    return(paste(as.character(head[[2L]]), as.character(head[[3L]]), sep = "::"))
  }

  ""
}

.np_progress_is_message_muffled <- function() {
  calls <- sys.calls()
  if (!length(calls)) {
    return(FALSE)
  }

  any(vapply(
    calls,
    function(call) {
      .np_progress_call_name(call) %in% c(
        "suppressMessages",
        "base::suppressMessages",
        "suppressPackageStartupMessages",
        "base::suppressPackageStartupMessages"
      )
    },
    logical(1L)
  ))
}

.np_progress_resolve_message_muffling <- function(state) {
  if (!identical(state$renderer, "single_line")) {
    return(state)
  }

  if (isTRUE(state$message_muffled_checked)) {
    return(state)
  }

  state$message_muffled <- isTRUE(.np_progress_is_message_muffled())
  state$message_muffled_checked <- TRUE
  state
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

.np_progress_make_registry <- function() {
  env <- new.env(parent = emptyenv())
  env$next_id <- 0L
  env$active_id <- NULL
  env
}

.np_progress_registry <- .np_progress_make_registry()

.np_progress_is_rstudio_console <- function() {
  identical(.Platform$GUI, "RStudio") || identical(Sys.getenv("RSTUDIO"), "1")
}

.np_progress_reset_registry <- function() {
  .np_progress_registry$next_id <- 0L
  .np_progress_registry$active_id <- NULL
  invisible(NULL)
}

.np_progress_next_id <- function() {
  .np_progress_registry$next_id <- as.integer(.np_progress_registry$next_id) + 1L
  sprintf("progress-%d", .np_progress_registry$next_id)
}

.np_progress_capability <- function(domain = "general") {
  interactive_console <- isTRUE(.np_progress_is_interactive())
  rstudio_console <- isTRUE(.np_progress_is_rstudio_console())

  list(
    domain = domain,
    interactive = interactive_console,
    rstudio = rstudio_console,
    single_line_viable = interactive_console
  )
}

.np_progress_single_line_surfaces <- function() {
  c("bandwidth", "plot_activity", "plot_bounded", "bootstrap", "lag", "iv_solve")
}

.np_progress_renderer_for_surface <- function(surface, capability) {
  if (surface %in% .np_progress_single_line_surfaces() &&
      isTRUE(capability$single_line_viable)) {
    return("single_line")
  }

  "legacy"
}

.np_progress_claim_owner <- function(session_id) {
  active_id <- .np_progress_registry$active_id
  if (is.null(active_id) || identical(active_id, session_id)) {
    .np_progress_registry$active_id <- session_id
    return(TRUE)
  }

  FALSE
}

.np_progress_release_owner <- function(session_id) {
  if (identical(.np_progress_registry$active_id, session_id)) {
    .np_progress_registry$active_id <- NULL
  }

  invisible(NULL)
}

.np_progress_make_snapshot <- function(state, line, event, now, done = NULL, detail = NULL, render_line = line) {
  list(
    id = state$id,
    surface = state$surface,
    renderer = state$renderer,
    label = state$label,
    kind = if (isTRUE(state$known_total)) "known" else "unknown",
    current = done,
    total = state$total,
    detail = detail,
    line = line,
    render_line = render_line,
    event = event,
    started_at = state$started,
    now = now,
    last_width = state$last_render_width
  )
}

.np_progress_render_legacy <- function(snapshot, event = c("render", "finish", "abort")) {
  event <- match.arg(event)
  .np_progress_emit(snapshot$line)
  invisible(snapshot)
}

.np_progress_single_line_connection <- function() {
  stderr()
}

.np_progress_has_tty <- function() {
  any(vapply(
    list(stdin(), stdout(), stderr()),
    function(con) {
      tryCatch(isTRUE(isatty(con)), error = function(...) FALSE)
    },
    logical(1L)
  ))
}

.np_progress_parse_width <- function(value) {
  value <- suppressWarnings(as.integer(value)[1L])
  if (!is.finite(value) || is.na(value) || value < 20L) {
    return(NA_integer_)
  }

  value
}

.np_progress_terminal_width_probe <- function() {
  if (!isTRUE(.np_progress_is_interactive()) || isTRUE(.np_progress_is_rstudio_console())) {
    return(NA_integer_)
  }
  if (!isTRUE(.np_progress_has_tty())) {
    return(NA_integer_)
  }

  if (identical(.Platform$OS.type, "windows")) {
    return(.np_progress_parse_width(Sys.getenv("COLUMNS", "")))
  }

  probe <- tryCatch(
    suppressWarnings(system2(
      "sh",
      c("-c", "stty size < /dev/tty 2>/dev/null"),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(...) character()
  )
  if (length(probe)) {
    fields <- strsplit(trimws(probe[[1L]]), "[[:space:]]+")[[1L]]
    width <- .np_progress_parse_width(fields[[length(fields)]])
    if (!is.na(width)) {
      return(width)
    }
  }

  probe <- tryCatch(
    suppressWarnings(system2("tput", "cols", stdout = TRUE, stderr = TRUE)),
    error = function(...) character()
  )
  if (length(probe)) {
    width <- .np_progress_parse_width(probe[[1L]])
    if (!is.na(width)) {
      return(width)
    }
  }

  .np_progress_parse_width(Sys.getenv("COLUMNS", ""))
}

.np_progress_output_width <- function() {
  width <- if (isTRUE(.np_progress_is_rstudio_console())) {
    .np_progress_parse_width(getOption("width", 80))
  } else {
    .np_progress_terminal_width_probe()
  }
  if (is.na(width)) {
    width <- .np_progress_parse_width(getOption("width", 80))
  }
  if (is.na(width)) {
    width <- 80L
  }

  reserve <- if (isTRUE(.np_progress_is_rstudio_console())) 4L else 0L
  max(20L, width - reserve)
}

.np_progress_ellipsize_middle <- function(text, max_width) {
  max_width <- suppressWarnings(as.integer(max_width)[1L])
  if (is.na(max_width) || max_width < 1L) {
    return("")
  }

  if (nchar(text, type = "width") <= max_width) {
    return(text)
  }

  if (max_width <= 3L) {
    return(substr("...", 1L, max_width))
  }

  keep_right <- max(1L, floor((max_width - 3L) / 2L))
  keep_left <- max_width - 3L - keep_right
  paste0(
    substr(text, 1L, keep_left),
    "...",
    substr(text, nchar(text, type = "chars") - keep_right + 1L, nchar(text, type = "chars"))
  )
}

.np_progress_compact_single_line <- function(line, max_width) {
  if (!is.character(line) || length(line) != 1L || is.na(line)) {
    return(line)
  }

  max_width <- suppressWarnings(as.integer(max_width)[1L])
  if (is.na(max_width) || max_width < 1L) {
    return(line)
  }

  if (nchar(line, type = "width") <= max_width) {
    return(line)
  }

  compacted <- line
  replacements <- list(
    c(" smooth coefficient ", " smooth coef "),
    c(" coefficient ", " coef "),
    c(" single-index ", " SI "),
    c(" regression ", " reg "),
    c(" density ", " dens "),
    c(" distribution ", " dist "),
    c(" bandwidth", " bw")
  )

  for (rule in replacements) {
    if (nchar(compacted, type = "width") <= max_width) {
      break
    }
    compacted <- gsub(rule[[1L]], rule[[2L]], compacted, fixed = TRUE)
  }

  compacted
}

.np_progress_fit_single_line <- function(line, max_width = .np_progress_output_width()) {
  if (!is.character(line) || length(line) != 1L || is.na(line)) {
    return(line)
  }

  max_width <- suppressWarnings(as.integer(max_width)[1L])
  if (is.na(max_width) || max_width < 1L) {
    return(line)
  }

  if (nchar(line, type = "width") <= max_width) {
    return(line)
  }

  line <- .np_progress_compact_single_line(line, max_width = max_width)
  if (nchar(line, type = "width") <= max_width) {
    return(line)
  }

  detail_pos <- regexpr(": ", line, fixed = TRUE)[1L]
  if (detail_pos > 0L) {
    without_detail <- substr(line, 1L, detail_pos - 1L)
    if (nchar(without_detail, type = "width") <= max_width) {
      return(without_detail)
    }
    line <- without_detail
  }

  .np_progress_ellipsize_middle(line, max_width = max_width)
}

.np_progress_render_single_line <- function(snapshot, event = c("render", "finish", "abort")) {
  event <- match.arg(event)
  render_line <- snapshot$render_line
  con <- .np_progress_single_line_connection()
  width <- nchar(render_line, type = "width")
  if (identical(event, "finish")) {
    clear_width <- max(snapshot$last_width, width)
    clear_line <- if (clear_width > 0L) strrep(" ", clear_width) else ""
    base::cat("\r", clear_line, "\r", file = con, sep = "")
    flush(con)
    flush.console()
    return(invisible(snapshot))
  }

  pad <- max(0L, snapshot$last_width - width)
  suffix <- if (pad > 0L) paste(rep(" ", pad), collapse = "") else ""

  base::cat("\r", render_line, suffix, file = con, sep = "")
  if (identical(event, "abort")) {
    base::cat("\n", file = con, sep = "")
  }
  flush(con)
  flush.console()
  invisible(snapshot)
}

.np_progress_render <- function(state, line, event, now, done = NULL, detail = NULL) {
  if (!isTRUE(state$enabled) || !isTRUE(state$visible)) {
    return(state)
  }

  state <- .np_progress_resolve_message_muffling(state)
  if (isTRUE(state$message_muffled)) {
    return(state)
  }

  render_line <- if (identical(state$renderer, "single_line")) {
    .np_progress_fit_single_line(line)
  } else {
    line
  }

  snapshot <- .np_progress_make_snapshot(
    state = state,
    line = line,
    render_line = render_line,
    event = event,
    now = now,
    done = done,
    detail = detail
  )

  if (identical(state$renderer, "single_line")) {
    .np_progress_render_single_line(
      snapshot = snapshot,
      event = if (identical(event, "finish")) "finish" else if (identical(event, "abort")) "abort" else "render"
    )
  } else {
    .np_progress_render_legacy(
      snapshot = snapshot,
      event = if (identical(event, "finish")) "finish" else if (identical(event, "abort")) "abort" else "render"
    )
  }

  state$rendered <- TRUE
  state$last_render_width <- nchar(render_line, type = "width")
  state$last_line <- line
  state
}

.np_progress_show_now <- function(state, done = NULL, detail = NULL) {
  if (is.null(state) || !isTRUE(state$enabled) || !isTRUE(state$visible)) {
    return(state)
  }

  if (!is.null(done)) {
    state$last_done <- done
  }

  now <- state$started
  line <- .np_progress_format_line(
    state = state,
    done = state$last_done,
    detail = detail,
    now = now
  )
  state <- .np_progress_render(
    state = state,
    line = line,
    event = "start",
    now = now,
    done = state$last_done,
    detail = detail
  )
  state$start_note_pending <- FALSE
  state$last_emit <- now
  state$last_emitted_done <- state$last_done
  state$last_emitted_detail <- detail
  state
}

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

.np_progress_format_line <- function(state, done = NULL, detail = NULL, now = .np_progress_now()) {
  if (isTRUE(state$known_total)) {
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
}

.np_progress_begin <- function(label, total = NULL, domain = "general", surface = domain) {
  known_total <- !is.null(total)
  throttle_sec <- if (known_total) 0.5 else 2.0
  started <- .np_progress_now()
  capability <- .np_progress_capability(domain = domain)
  session_id <- .np_progress_next_id()
  enabled <- .np_progress_enabled(domain = domain)
  renderer <- .np_progress_renderer_for_surface(surface = surface, capability = capability)
  visible <- if (!isTRUE(enabled)) {
    FALSE
  } else {
    isTRUE(.np_progress_claim_owner(session_id))
  }

  list(
    id = session_id,
    enabled = enabled,
    visible = visible,
    pkg_prefix = .np_progress_pkg_prefix(),
    label = label,
    surface = surface,
    total = total,
    known_total = known_total,
    started = started,
    last_emit = started - throttle_sec,
    throttle_sec = throttle_sec,
    last_done = if (known_total) 0 else NULL,
    domain = domain,
    renderer = renderer,
    capability = capability,
    rendered = FALSE,
    start_note = sprintf("%s %s...", .np_progress_pkg_prefix(), label),
    start_note_pending = TRUE,
    start_note_grace_sec = .np_progress_start_grace_sec(known_total = known_total, domain = domain),
    last_line = NULL,
    last_render_width = 0L,
    last_emitted_done = NULL,
    last_emitted_detail = NULL,
    message_muffled = FALSE,
    message_muffled_checked = FALSE
  )
}

.np_progress_maybe_emit_start_note <- function(state, now = .np_progress_now()) {
  if (!isTRUE(state$enabled) || !isTRUE(state$visible) || !isTRUE(state$start_note_pending)) {
    return(state)
  }

  if ((now - state$started) < state$start_note_grace_sec) {
    return(state)
  }

  state <- .np_progress_render(
    state = state,
    line = state$start_note,
    event = "start",
    now = now
  )
  state$start_note_pending <- FALSE

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

  line <- .np_progress_format_line(
    state = state,
    done = if (is.null(done)) state$last_done else done,
    detail = detail,
    now = now
  )

  if (!identical(line, state$last_line)) {
    state <- .np_progress_render(
      state = state,
      line = line,
      event = "update",
      now = now,
      done = if (is.null(done)) state$last_done else done,
      detail = detail
    )
    state$last_emit <- now
    state$last_emitted_done <- if (is.null(done)) state$last_done else done
    state$last_emitted_detail <- detail
  }

  state
}

.np_progress_end <- function(state, detail = NULL) {
  if (!isTRUE(state$enabled)) {
    .np_progress_release_owner(state$id)
    return(invisible(state))
  }

  now <- .np_progress_now()
  if (is.null(state$last_line) &&
      (now - state$started) < state$start_note_grace_sec) {
    .np_progress_release_owner(state$id)
    return(invisible(state))
  }

  must_clear <- identical(state$renderer, "single_line") && isTRUE(state$rendered)

  if (isTRUE(state$known_total)) {
    done <- if (is.null(state$total)) state$last_done else state$total
    line <- .np_progress_format_line(state = state, done = done, detail = detail, now = now)
    if (isTRUE(must_clear) ||
        !(identical(done, state$last_emitted_done) && identical(detail, state$last_emitted_detail))) {
      state <- .np_progress_render(
        state = state,
        line = line,
        event = "finish",
        now = now,
        done = done,
        detail = detail
      )
      state$last_emit <- now
      state$last_done <- done
      state$last_emitted_done <- done
      state$last_emitted_detail <- detail
    }
    .np_progress_release_owner(state$id)
    return(invisible(state))
  }

  if (!is.null(state$last_done)) {
    line <- .np_progress_format_line(state = state, done = state$last_done, detail = detail, now = now)
    if (isTRUE(must_clear) ||
        !(identical(state$last_done, state$last_emitted_done) && identical(detail, state$last_emitted_detail))) {
      state <- .np_progress_render(
        state = state,
        line = line,
        event = "finish",
        now = now,
        done = state$last_done,
        detail = detail
      )
      state$last_emit <- now
      state$last_emitted_done <- state$last_done
      state$last_emitted_detail <- detail
    }
  }

  if (isTRUE(must_clear) && is.null(state$last_done) && !is.null(state$last_line)) {
    state <- .np_progress_render(
      state = state,
      line = state$last_line,
      event = "finish",
      now = now,
      detail = detail
    )
  }

  .np_progress_release_owner(state$id)
  invisible(state)
}

.np_progress_abort <- function(state, detail = NULL) {
  if (!isTRUE(state$enabled)) {
    .np_progress_release_owner(state$id)
    return(invisible(state))
  }

  now <- .np_progress_now()
  if (!isTRUE(state$rendered)) {
    .np_progress_release_owner(state$id)
    return(invisible(state))
  }

  line <- if (!is.null(detail)) paste0(state$pkg_prefix, " ", detail) else state$last_line
  state <- .np_progress_render(
    state = state,
    line = line,
    event = "abort",
    now = now,
    done = state$last_done,
    detail = detail
  )
  .np_progress_release_owner(state$id)
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

  upgraded <- .np_progress_begin(
    label = label,
    total = total,
    domain = state$domain,
    surface = state$surface
  )
  upgraded$id <- state$id
  upgraded$visible <- state$visible
  upgraded$started <- state$started
  upgraded$last_emit <- state$started - upgraded$throttle_sec
  upgraded$renderer <- state$renderer
  upgraded$capability <- state$capability
  upgraded$rendered <- state$rendered
  upgraded$start_note_pending <- state$start_note_pending
  upgraded$last_line <- state$last_line
  upgraded$last_render_width <- state$last_render_width
  upgraded$last_emitted_done <- state$last_emitted_done
  upgraded$last_emitted_detail <- state$last_emitted_detail
  upgraded$message_muffled <- state$message_muffled
  upgraded$message_muffled_checked <- state$message_muffled_checked
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

.np_progress_signal_from_c <- function(event, surface, current = NULL, total = NULL) {
  event <- as.character(event)[1L]
  surface <- as.character(surface)[1L]

  if (is.na(event) || is.na(surface)) {
    return(invisible(FALSE))
  }

  if (identical(event, "bandwidth_multistart_step") && identical(surface, "bandwidth")) {
    .np_progress_bandwidth_multistart_step(done = current, total = total)
    return(invisible(TRUE))
  }

  invisible(FALSE)
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
      domain = "general",
      surface = "bandwidth"
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
