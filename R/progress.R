.np_progress_is_interactive <- function() {
  interactive()
}

.np_progress_is_master <- function() {
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE))) {
    return(TRUE)
  }

  tryCatch(isTRUE(mpi.is.master()), error = function(e) TRUE)
}

.np_progress_now <- function() {
  as.numeric(proc.time()[["elapsed"]])
}

.np_progress_enabled <- function(domain = "general") {
  (isTRUE(.np_progress_runtime$force_enabled) || isTRUE(getOption("np.messages", TRUE))) &&
    isTRUE(.np_progress_is_interactive()) &&
    isTRUE(.np_progress_is_master())
}

.np_io_pkg_prefix <- function() {
  "[npRmpi]"
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
  env$bandwidth_context_label <- NULL
  env$fit_state <- NULL
  env$iv_depth <- 0L
  env$iv_label <- NULL
  env$iv_state <- NULL
  env$iv_old_messages <- NULL
  env$force_enabled <- FALSE
  env
})

.np_progress_bandwidth_set_context <- function(label = NULL) {
  if (is.null(label)) {
    .np_progress_runtime$bandwidth_context_label <- NULL
    return(invisible(NULL))
  }

  label <- as.character(label)[1L]
  if (is.na(label) || !nzchar(label)) {
    .np_progress_runtime$bandwidth_context_label <- NULL
    return(invisible(NULL))
  }

  .np_progress_runtime$bandwidth_context_label <- label
  invisible(NULL)
}

.np_progress_make_registry <- function() {
  env <- new.env(parent = emptyenv())
  env$next_id <- 0L
  env$active_id <- NULL
  env$last_single_line_width <- 0L
  env
}

.np_progress_registry <- .np_progress_make_registry()

.np_progress_is_rstudio_console <- function() {
  identical(.Platform$GUI, "RStudio") || identical(Sys.getenv("RSTUDIO"), "1")
}

.np_progress_reset_registry <- function() {
  .np_progress_registry$next_id <- 0L
  .np_progress_registry$active_id <- NULL
  .np_progress_registry$last_single_line_width <- 0L
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

.np_progress_single_line_supports_ansi <- function(con = .np_progress_single_line_connection()) {
  if (identical(.Platform$OS.type, "windows")) {
    return(FALSE)
  }

  tryCatch(isTRUE(isatty(con)), error = function(...) FALSE)
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
    c(" Preparing plot bootstrap ", " Prep plot "),
    c(" Computing regression plot fit", " Regression plot fit"),
    c(" Computing bootstrap pmzsd errors", " Bootstrap pmzsd"),
    c(" Constructing bootstrap all bands", " Bootstrap all bands"),
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

.np_progress_compact_bandwidth_line <- function(line, max_width) {
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

  matches <- regexec("^(.+Bandwidth selection) \\((.+)\\)$", line)
  capture <- regmatches(line, matches)[[1L]]
  if (length(capture) != 3L) {
    return(line)
  }

  prefix <- capture[[2L]]
  fields <- strsplit(capture[[3L]], ", ", fixed = TRUE)[[1L]]
  if (!length(fields)) {
    return(line)
  }

  render <- function(parts) {
    sprintf("%s (%s)", prefix, paste(parts, collapse = ", "))
  }

  abbreviate_fields <- function(parts) {
    compacted <- parts
    compacted <- sub("^multistart ([0-9]+/[0-9]+)$", "\\1", compacted)
    compacted <- sub("^iteration ", "iter ", compacted)
    compacted
  }

  drop_last_field <- function(parts) {
    if (length(parts) <= 1L) {
      parts
    } else {
      parts[seq_len(length(parts) - 1L)]
    }
  }

  fields_abbrev <- abbreviate_fields(fields)

  candidates <- list(
    render(fields),
    render(fields_abbrev)
  )

  for (candidate in candidates) {
    if (nchar(candidate, type = "width") <= max_width) {
      return(candidate)
    }
  }

  line
}

.np_progress_compact_plot_bootstrap_line <- function(line, max_width) {
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

  matches <- regexec(
    "^(.+(?:Plot bootstrap(?: \\([^)]+\\))?|Bootstrap (?:all bands|pmzsd)(?: \\([^)]+\\))?)) ([0-9]+/[0-9]+) \\((.+)\\)$",
    line
  )
  capture <- regmatches(line, matches)[[1L]]
  if (length(capture) != 4L) {
    return(line)
  }

  prefix <- capture[[2L]]
  counter <- capture[[3L]]
  fields <- strsplit(capture[[4L]], ", ", fixed = TRUE)[[1L]]
  if (!length(fields)) {
    return(line)
  }

  render <- function(prefix_text, parts) {
    sprintf("%s %s (%s)", prefix_text, counter, paste(parts, collapse = ", "))
  }

  compact_prefix <- prefix
  compact_prefix <- sub("\\(grad index ", "(grad idx ", compact_prefix)
  compact_prefix <- sub("\\(index ", "(idx ", compact_prefix)
  short_prefix <- compact_prefix
  short_prefix <- sub("Plot bootstrap", "Plot boot", short_prefix, fixed = TRUE)
  short_prefix <- sub("Bootstrap all bands", "Boot all bands", short_prefix, fixed = TRUE)
  short_prefix <- sub("Bootstrap pmzsd", "Boot pmzsd", short_prefix, fixed = TRUE)
  short_prefix <- sub("\\(surf ", "(s ", short_prefix)
  very_short_prefix <- short_prefix
  very_short_prefix <- sub("Boot all bands", "Bands", very_short_prefix, fixed = TRUE)
  very_short_prefix <- sub("Boot pmzsd", "Pmzsd", very_short_prefix, fixed = TRUE)

  compact_fields <- function(parts) {
    compacted <- parts
    compacted <- sub("^elapsed ", "elap ", compacted)
    compacted
  }

  drop_last_field <- function(parts) {
    if (length(parts) <= 1L) {
      parts
    } else {
      parts[seq_len(length(parts) - 1L)]
    }
  }

  fields_compact <- compact_fields(fields)
  fields_drop_eta <- drop_last_field(fields_compact)
  fields_drop_eta_pct <- drop_last_field(fields_drop_eta)

  candidates <- list(
    render(prefix, fields),
    render(compact_prefix, fields),
    render(compact_prefix, fields_compact),
    render(short_prefix, fields_compact),
    render(short_prefix, fields_drop_eta),
    render(short_prefix, fields_drop_eta_pct),
    render(very_short_prefix, fields_drop_eta),
    render(very_short_prefix, fields_drop_eta_pct),
    sprintf("%s %s", very_short_prefix, counter)
  )

  for (candidate in candidates) {
    if (nchar(candidate, type = "width") <= max_width) {
      return(candidate)
    }
  }

  line
}

.np_progress_compact_degree_line <- function(line, max_width) {
  if (!is.character(line) || length(line) != 1L || is.na(line)) {
    return(line)
  }

  max_width <- suppressWarnings(as.integer(max_width)[1L])
  if (is.na(max_width) || max_width < 1L) {
    return(line)
  }

  if (!grepl("Selecting polynomial degree and bandwidth", line, fixed = TRUE) &&
      !grepl("Selecting polynomial degree and bw", line, fixed = TRUE) &&
      !grepl("Selecting degree and bandwidth", line, fixed = TRUE) &&
      !grepl("Refining bandwidth", line, fixed = TRUE)) {
    return(line)
  }

  compact <- function(text) {
    out <- text
    out <- sub("Selecting polynomial degree and bandwidth", "Degree/bw search", out, fixed = TRUE)
    out <- sub("Selecting polynomial degree and bw", "Degree/bw search", out, fixed = TRUE)
    out <- sub("Selecting degree and bandwidth", "Degree/bw search", out, fixed = TRUE)
    out <- sub("Refining bandwidth", "Refining bw", out, fixed = TRUE)
    out <- gsub("elapsed ", "elap ", out, fixed = TRUE)
    out <- gsub("restart ", "r ", out, fixed = TRUE)
    out <- gsub("cycle ", "cy ", out, fixed = TRUE)
    out <- gsub("coord ", "c ", out, fixed = TRUE)
    out <- gsub("degree ", "deg ", out, fixed = TRUE)
    out
  }

  candidates <- list(
    line,
    compact(line)
  )

  for (candidate in candidates) {
    if (nchar(candidate, type = "width") <= max_width) {
      return(candidate)
    }
  }

  line
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

  line <- .np_progress_compact_bandwidth_line(line, max_width = max_width)
  if (nchar(line, type = "width") <= max_width) {
    return(line)
  }

  line <- .np_progress_compact_plot_bootstrap_line(line, max_width = max_width)
  if (nchar(line, type = "width") <= max_width) {
    return(line)
  }

  line <- .np_progress_compact_degree_line(line, max_width = max_width)
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
  use_ansi <- isTRUE(.np_progress_single_line_supports_ansi(con))
  width <- nchar(render_line, type = "width")
  global_width <- suppressWarnings(as.integer(.np_progress_registry$last_single_line_width)[1L])
  if (is.na(global_width) || global_width < 0L) {
    global_width <- 0L
  }
  if (identical(event, "finish")) {
    if (use_ansi) {
      .np_progress_registry$last_single_line_width <- 0L
      base::cat("\r\033[2K\r", file = con, sep = "")
      flush(con)
      flush.console()
      return(invisible(snapshot))
    }
    clear_width <- max(snapshot$last_width, width, global_width)
    clear_line <- if (clear_width > 0L) strrep(" ", clear_width) else ""
    .np_progress_registry$last_single_line_width <- 0L
    base::cat("\r", clear_line, "\r", file = con, sep = "")
    flush(con)
    flush.console()
    return(invisible(snapshot))
  }

  pad <- max(0L, max(snapshot$last_width, global_width) - width)
  suffix <- if (pad > 0L) paste(rep(" ", pad), collapse = "") else ""

  if (use_ansi) {
    base::cat("\r\033[2K", render_line, file = con, sep = "")
  } else {
    base::cat("\r", render_line, suffix, file = con, sep = "")
  }
  if (identical(event, "abort")) {
    .np_progress_registry$last_single_line_width <- 0L
    base::cat("\n", file = con, sep = "")
  } else {
    .np_progress_registry$last_single_line_width <- width
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

.np_progress_bandwidth_enhanced_enabled <- function() {
  isTRUE(getOption("np.progress.bandwidth.enhanced", FALSE))
}

.np_progress_iv_enhanced_state <- function(state) {
  isTRUE(state$iv_progress_common)
}

.np_progress_bandwidth_enhanced_state <- function(state) {
  isTRUE(state$bandwidth_progress_common)
}

.np_progress_bandwidth_coordinator_active <- function(state) {
  isTRUE(state$bandwidth_coordinator_active)
}

.np_progress_bandwidth_title <- function() {
  "Bandwidth selection"
}

.np_progress_iv_title <- function() {
  "IV regression"
}

.np_progress_iv_initialize_state <- function(state) {
  state$iv_progress_common <- TRUE
  state$iv_object_label <- NULL
  state$iv_iteration <- NULL
  state$start_note <- sprintf(
    "%s %s...",
    state$pkg_prefix,
    .np_progress_iv_title()
  )
  state
}

.np_progress_iv_format_line <- function(state, now = .np_progress_now()) {
  elapsed <- max(0, now - state$started)
  fields <- character()

  if (!is.null(state$iv_object_label) && nzchar(state$iv_object_label)) {
    fields <- c(fields, state$iv_object_label)
  }

  if (!is.null(state$iv_iteration)) {
    fields <- c(fields, sprintf("iteration %s", format(state$iv_iteration)))
  }

  fields <- c(fields, sprintf("elapsed %ss", .np_progress_fmt_num(elapsed)))

  sprintf(
    "%s %s (%s)",
    state$pkg_prefix,
    .np_progress_iv_title(),
    paste(fields, collapse = ", ")
  )
}

.np_progress_bandwidth_initialize_state <- function(state) {
  state$bandwidth_progress_common <- TRUE
  state$bandwidth_nmulti_total <- 1L
  state$bandwidth_multistart_current <- 1L
  state$bandwidth_multistart_completed <- 0L
  state$bandwidth_current_started <- state$started
  state$bandwidth_multistart_durations <- numeric()
  state$bandwidth_mode <- "iteration"
  state$bandwidth_last_iteration <- NULL
  state$bandwidth_estimate_last_emit <- NULL
  state$bandwidth_coordinator_active <- FALSE
  state$bandwidth_coordinator_groups <- NULL
  state$bandwidth_coordinator_local_total <- NULL
  state$bandwidth_coordinator_group_index <- NULL
  state$bandwidth_coordinator_group_label <- NULL
  state$bandwidth_coordinator_offset <- 0L
  state$bandwidth_coordinator_local_current <- 1L
  state$start_note <- sprintf(
    "%s %s...",
    state$pkg_prefix,
    .np_progress_bandwidth_title()
  )
  state
}

.np_progress_bandwidth_set_coordinator <- function(total_groups, local_total) {
  state <- .np_progress_runtime$bandwidth_state
  total_groups <- suppressWarnings(as.integer(total_groups)[1L])
  local_total <- suppressWarnings(as.integer(local_total)[1L])

  if (is.null(state) ||
      !.np_progress_bandwidth_enhanced_state(state) ||
      is.na(total_groups) || total_groups < 1L ||
      is.na(local_total) || local_total < 1L) {
    return(invisible(NULL))
  }

  state$bandwidth_coordinator_active <- TRUE
  state$bandwidth_coordinator_groups <- total_groups
  state$bandwidth_coordinator_local_total <- local_total
  state$bandwidth_coordinator_group_index <- 1L
  state$bandwidth_coordinator_offset <- 0L
  state$bandwidth_nmulti_total <- total_groups * local_total
  state$bandwidth_multistart_current <- 1L
  state$bandwidth_multistart_completed <- 0L
  .np_progress_runtime$bandwidth_state <- state
  invisible(NULL)
}

.np_progress_bandwidth_set_coordinator_group <- function(group_index, group_label = NULL) {
  state <- .np_progress_runtime$bandwidth_state
  group_index <- suppressWarnings(as.integer(group_index)[1L])

  if (is.null(state) ||
      !.np_progress_bandwidth_coordinator_active(state) ||
      is.na(group_index) || group_index < 1L) {
    return(invisible(NULL))
  }

  local_total <- suppressWarnings(as.integer(state$bandwidth_coordinator_local_total)[1L])
  total <- suppressWarnings(as.integer(state$bandwidth_nmulti_total)[1L])
  if (is.na(local_total) || local_total < 1L || is.na(total) || total < 1L) {
    return(invisible(NULL))
  }

  offset <- (group_index - 1L) * local_total
  state$bandwidth_coordinator_group_index <- group_index
  state$bandwidth_coordinator_group_label <- if (is.null(group_label)) NULL else as.character(group_label)[1L]
  state$bandwidth_coordinator_offset <- offset
  state$bandwidth_coordinator_local_current <- 1L
  state$bandwidth_multistart_completed <- min(offset, total)
  state$bandwidth_multistart_current <- min(total, offset + 1L)
  state$bandwidth_current_started <- .np_progress_now()
  state$last_done <- NULL
  if (offset > 0L) {
    state$bandwidth_mode <- "complete_estimate"
  }
  if (!is.na(local_total) && local_total > 1L) {
    state$start_note <- if (!is.null(state$bandwidth_coordinator_group_label)) {
      sprintf(
        "%s %s (%s, multistart 1/%s)",
        state$pkg_prefix,
        .np_progress_bandwidth_title(),
        state$bandwidth_coordinator_group_label,
        format(local_total)
      )
    } else {
      sprintf(
        "%s %s (multistart 1/%s)",
        state$pkg_prefix,
        .np_progress_bandwidth_title(),
        format(local_total)
      )
    }
  } else {
    state$start_note <- if (!is.null(state$bandwidth_coordinator_group_label)) {
      sprintf(
        "%s %s (%s)",
        state$pkg_prefix,
        .np_progress_bandwidth_title(),
        state$bandwidth_coordinator_group_label
      )
    } else {
      sprintf(
        "%s %s...",
        state$pkg_prefix,
        .np_progress_bandwidth_title()
      )
    }
  }
  .np_progress_runtime$bandwidth_state <- state
  invisible(NULL)
}

.np_progress_bandwidth_register_nmulti <- function(state, total) {
  total <- suppressWarnings(as.integer(total)[1L])
  if (is.null(state) || is.na(total) || total < 1L) {
    return(state)
  }

  if (.np_progress_bandwidth_coordinator_active(state)) {
    state$bandwidth_coordinator_local_total <- total
    groups <- suppressWarnings(as.integer(state$bandwidth_coordinator_groups)[1L])
    if (!is.na(groups) && groups >= 1L) {
      state$bandwidth_nmulti_total <- groups * total
    }
    return(state)
  }

  state$bandwidth_nmulti_total <- total
  state$bandwidth_multistart_current <- min(
    max(1L, suppressWarnings(as.integer(state$bandwidth_multistart_current)[1L])),
    total
  )
  if (total > 1L) {
    state$start_note <- sprintf(
      "%s %s (multistart 1/%s)",
      state$pkg_prefix,
      .np_progress_bandwidth_title(),
      format(total)
    )
  } else {
    state$start_note <- sprintf(
      "%s %s...",
      state$pkg_prefix,
      .np_progress_bandwidth_title()
    )
  }
  state
}

.np_progress_bandwidth_estimated_total_time <- function(state) {
  durations <- state$bandwidth_multistart_durations
  total <- suppressWarnings(as.integer(state$bandwidth_nmulti_total)[1L])
  if (!length(durations) || is.na(total) || total < 1L) {
    return(NA_real_)
  }

  mean(durations) * total
}

.np_progress_bandwidth_format_iteration <- function(state, done = NULL, now = .np_progress_now()) {
  elapsed <- max(0, now - state$started)
  total <- suppressWarnings(as.integer(state$bandwidth_nmulti_total)[1L])
  current <- suppressWarnings(as.integer(state$bandwidth_multistart_current)[1L])
  fields <- character()

  if (.np_progress_bandwidth_coordinator_active(state)) {
    label <- state$bandwidth_coordinator_group_label
    local_total <- suppressWarnings(as.integer(state$bandwidth_coordinator_local_total)[1L])
    local_current <- suppressWarnings(as.integer(state$bandwidth_coordinator_local_current)[1L])
    if (!is.null(label) && nzchar(label)) {
      fields <- c(fields, label)
    }
    if (!is.na(local_total) && local_total > 1L && !is.na(local_current) && local_current >= 1L) {
      fields <- c(fields, sprintf("multistart %s/%s", format(local_current), format(local_total)))
    }
  } else if (!is.na(total) && total > 1L && !is.na(current) && current >= 1L) {
    fields <- c(fields, sprintf("multistart %s/%s", format(current), format(total)))
  }

  if (!is.null(done)) {
    fields <- c(fields, sprintf("iteration %s", format(done)))
  }

  fields <- c(fields, sprintf("elapsed %ss", .np_progress_fmt_num(elapsed)))

  sprintf(
    "%s %s (%s)",
    state$pkg_prefix,
    .np_progress_bandwidth_title(),
    paste(fields, collapse = ", ")
  )
}

.np_progress_bandwidth_format_estimate <- function(state, now = .np_progress_now()) {
  elapsed <- max(0, now - state$started)
  total <- suppressWarnings(as.integer(state$bandwidth_nmulti_total)[1L])
  completed <- suppressWarnings(as.integer(state$bandwidth_multistart_completed)[1L])
  current <- suppressWarnings(as.integer(state$bandwidth_multistart_current)[1L])
  if (is.na(total) || total < 1L || completed < 1L) {
    return(.np_progress_bandwidth_format_iteration(state = state, done = state$last_done, now = now))
  }

  est.total <- .np_progress_bandwidth_estimated_total_time(state)
  base.share <- 100 * completed / total
  pct <- base.share
  eta <- 0

  if (is.finite(est.total) && !is.na(est.total) && est.total > 0) {
    pct.calc <- 100 * elapsed / est.total
    pct.upper <- if (completed < total) 99.9 else 100.0
    pct <- max(base.share, min(pct.upper, pct.calc))
    if (completed < total) {
      eta <- max(0, est.total - elapsed)
    }
  }

  fields <- character()

  if (.np_progress_bandwidth_coordinator_active(state)) {
    label <- state$bandwidth_coordinator_group_label
    local_total <- suppressWarnings(as.integer(state$bandwidth_coordinator_local_total)[1L])
    local_current <- suppressWarnings(as.integer(state$bandwidth_coordinator_local_current)[1L])
    if (!is.null(label) && nzchar(label)) {
      fields <- c(fields, label)
    }
    if (!is.na(local_current) && local_current >= 1L && !is.na(local_total) && local_total >= 1L) {
      fields <- c(fields, sprintf("multistart %s/%s", format(local_current), format(local_total)))
    }
  } else if (!is.na(current) && current >= 1L) {
    fields <- c(fields, sprintf("multistart %s/%s", format(current), format(total)))
  }

  if (completed < total && !is.null(state$last_done)) {
    fields <- c(fields, sprintf("iteration %s", format(state$last_done)))
  }

  fields <- c(
    fields,
    sprintf("elapsed %ss", .np_progress_fmt_num(elapsed)),
    sprintf("%s%%", .np_progress_fmt_num(pct)),
    sprintf("eta %ss", .np_progress_fmt_num(eta))
  )

  sprintf(
    "%s %s (%s)",
    state$pkg_prefix,
    .np_progress_bandwidth_title(),
    paste(fields, collapse = ", ")
  )
}

.np_progress_bandwidth_format_line <- function(state, done = NULL, now = .np_progress_now()) {
  mode <- state$bandwidth_mode
  if (identical(mode, "complete_estimate")) {
    .np_progress_bandwidth_format_estimate(state = state, now = now)
  } else {
    .np_progress_bandwidth_format_iteration(
      state = state,
      done = if (is.null(done)) state$last_done else done,
      now = now
    )
  }
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
  if (is.function(state$unknown_total_fields)) {
    fields <- state$unknown_total_fields(
      state = state,
      done = done,
      detail = detail,
      now = now
    )
    if (length(fields)) {
      return(sprintf(
        "%s %s (%s)",
        state$pkg_prefix,
        state$label,
        paste(as.character(fields), collapse = ", ")
      ))
    }
  }

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
  if (.np_progress_iv_enhanced_state(state)) {
    return(.np_progress_iv_format_line(
      state = state,
      now = now
    ))
  }

  if (.np_progress_bandwidth_enhanced_state(state)) {
    return(.np_progress_bandwidth_format_line(
      state = state,
      done = done,
      now = now
    ))
  }

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
    unknown_total_fields = NULL,
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

.np_progress_step_at <- function(state, now, done = NULL, detail = NULL, force = FALSE) {
  if (!is.null(done)) {
    state$last_done <- done
  }

  if (!isTRUE(state$enabled)) {
    return(state)
  }

  state <- .np_progress_maybe_emit_start_note(state = state, now = now)
  if (isTRUE(state$start_note_pending)) {
    return(state)
  }

  if (!isTRUE(force) && (now - state$last_emit) < state$throttle_sec) {
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

.np_progress_step <- function(state, done = NULL, detail = NULL) {
  .np_progress_step_at(
    state = state,
    now = .np_progress_now(),
    done = done,
    detail = detail
  )
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

  if (.np_progress_bandwidth_enhanced_state(state)) {
    total <- suppressWarnings(as.integer(state$bandwidth_nmulti_total)[1L])
    if (!is.na(total) && total > 1L) {
      completed <- suppressWarnings(as.integer(state$bandwidth_multistart_completed)[1L])
      if (completed < total) {
        gap <- total - completed
        dur <- max(0, now - state$bandwidth_current_started)
        if (gap > 0L) {
          state$bandwidth_multistart_durations <- c(
            state$bandwidth_multistart_durations,
            rep.int(dur, gap)
          )
        }
        state$bandwidth_multistart_completed <- total
        state$bandwidth_multistart_current <- total
      }
      state$bandwidth_mode <- "complete_estimate"
      line <- .np_progress_bandwidth_format_estimate(state = state, now = now)
      if (isTRUE(must_clear) || !identical(line, state$last_line)) {
        state <- .np_progress_render(
          state = state,
          line = line,
          event = "finish",
          now = now,
          done = state$last_done,
          detail = detail
        )
      }
      .np_progress_release_owner(state$id)
      return(invisible(state))
    }
  }

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

.np_progress_iv_active <- function() {
  isTRUE(.np_progress_runtime$iv_depth > 0L)
}

.np_progress_iv_set_object <- function(label = NULL, iteration = NULL) {
  state <- .np_progress_runtime$iv_state

  if (is.null(state) || !.np_progress_iv_enhanced_state(state)) {
    return(invisible(NULL))
  }

  if (!is.null(label)) {
    label <- as.character(label)[1L]
    if (is.na(label) || !nzchar(label)) {
      label <- NULL
    }
  }

  if (!is.null(iteration)) {
    iteration <- suppressWarnings(as.integer(iteration)[1L])
    if (is.na(iteration) || iteration < 1L) {
      iteration <- NULL
    }
  }

  state$iv_object_label <- label
  state$iv_iteration <- iteration
  .np_progress_runtime$iv_state <- .np_progress_step_at(
    state = state,
    now = .np_progress_now(),
    done = iteration,
    force = TRUE
  )

  invisible(NULL)
}

.np_progress_iv_activity_step <- function(iteration = NULL) {
  state <- .np_progress_runtime$iv_state

  if (is.null(state) || !.np_progress_iv_enhanced_state(state)) {
    return(invisible(NULL))
  }

  if (!is.null(iteration)) {
    iteration <- suppressWarnings(as.integer(iteration)[1L])
    if (is.na(iteration) || iteration < 1L) {
      iteration <- NULL
    }
    state$iv_iteration <- iteration
  }

  .np_progress_runtime$iv_state <- .np_progress_step(
    state = state,
    done = state$iv_iteration
  )

  invisible(NULL)
}

.np_progress_iv_finish <- function() {
  state <- .np_progress_runtime$iv_state

  if (is.null(state)) {
    return(invisible(NULL))
  }

  .np_progress_runtime$iv_state <- .np_progress_end(state)
  invisible(NULL)
}

.np_progress_iv_clear <- function() {
  .np_progress_runtime$iv_depth <- 0L
  .np_progress_runtime$iv_label <- NULL
  .np_progress_runtime$iv_state <- NULL
  .np_progress_runtime$iv_old_messages <- NULL
  .np_progress_runtime$force_enabled <- FALSE
  invisible(NULL)
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

  # NOMAD restart search is unknown-bound even when inner helpers expose a
  # local multistart count; do not synthesize percent/ETA from that count.
  if (!is.null(state$nomad_nmulti) && is.function(state$unknown_total_fields)) {
    return(invisible(NULL))
  }

  if (.np_progress_bandwidth_enhanced_state(state)) {
    .np_progress_runtime$bandwidth_state <- .np_progress_bandwidth_register_nmulti(
      state = state,
      total = total
    )
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

  if (.np_progress_bandwidth_enhanced_state(state)) {
    done <- suppressWarnings(as.integer(done)[1L])
    if (is.null(state) || is.na(done) || done < 1L) {
      return(invisible(NULL))
    }

    total <- suppressWarnings(as.integer(state$bandwidth_nmulti_total)[1L])
    local_total <- if (.np_progress_bandwidth_coordinator_active(state)) {
      suppressWarnings(as.integer(state$bandwidth_coordinator_local_total)[1L])
    } else {
      total
    }
    if (is.na(local_total) || local_total < 2L || is.na(total) || total < 2L) {
      return(invisible(NULL))
    }

    done <- min(done, local_total)
    now <- .np_progress_now()
    completed <- suppressWarnings(as.integer(state$bandwidth_multistart_completed)[1L])
    completed <- max(0L, completed)
    global_done <- if (.np_progress_bandwidth_coordinator_active(state)) {
      offset <- suppressWarnings(as.integer(state$bandwidth_coordinator_offset)[1L])
      offset + done
    } else {
      done
    }
    gap <- global_done - completed
    if (gap > 0L) {
      dur <- max(0, now - state$bandwidth_current_started)
      state$bandwidth_multistart_durations <- c(
        state$bandwidth_multistart_durations,
        rep.int(dur, gap)
      )
    }
    state$bandwidth_multistart_completed <- global_done
    state$bandwidth_multistart_current <- if (done < local_total) global_done + 1L else min(total, global_done)
    if (.np_progress_bandwidth_coordinator_active(state)) {
      state$bandwidth_coordinator_local_current <- if (done < local_total) done + 1L else local_total
    }
    if (done < local_total) {
      state$bandwidth_current_started <- now
      state$last_done <- NULL
    }
    state$bandwidth_mode <- "complete_estimate"
    state$last_emit <- -Inf
    .np_progress_runtime$bandwidth_state <- .np_progress_step_at(
      state = state,
      now = now,
      done = state$last_done,
      force = TRUE
    )
    return(invisible(NULL))
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

.np_progress_bandwidth_activity_step <- function(done = NULL) {
  state <- .np_progress_runtime$bandwidth_state

  if (.np_progress_bandwidth_enhanced_state(state)) {
    if (!is.null(done)) {
      done <- suppressWarnings(as.integer(done)[1L])
      if (is.na(done) || done < 1L) {
        done <- NULL
      }
    }
    .np_progress_runtime$bandwidth_state <- .np_progress_step(
      state = state,
      done = done
    )
    return(invisible(NULL))
  }

  if (is.null(state) || isTRUE(state$known_total)) {
    return(invisible(NULL))
  }

  if (!is.null(done)) {
    done <- suppressWarnings(as.integer(done)[1L])
    if (is.na(done) || done < 1L) {
      done <- NULL
    }
  }

  .np_progress_runtime$bandwidth_state <- .np_progress_step(
    state = state,
    done = done
  )

  invisible(NULL)
}

.np_fit_progress_begin <- function(label, total, handoff = FALSE, detail = NULL) {
  total <- suppressWarnings(as.integer(total)[1L])
  if (is.na(total) || total < 1L) {
    .np_progress_runtime$fit_state <- NULL
    return(invisible(NULL))
  }

  if (isTRUE(getOption("npRmpi.mpi.initialized", FALSE))) {
    attach.size <- tryCatch(as.integer(mpi.comm.size(1L)), error = function(e) NA_integer_)
    attach.rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) NA_integer_)

    if (!is.na(attach.size) && attach.size > 1L &&
        !is.na(attach.rank) && attach.rank != 0L) {
      .np_progress_runtime$fit_state <- NULL
      return(invisible(NULL))
    }
  }

  state <- .np_progress_begin(
    label = label,
    total = total,
    domain = "bandwidth",
    surface = "bandwidth"
  )

  if (isTRUE(handoff)) {
    state <- .np_progress_show_now(
      state = state,
      done = 0L,
      detail = detail
    )
  }

  .np_progress_runtime$fit_state <- state
  invisible(state)
}

.np_fit_progress_step <- function(done = NULL, detail = NULL) {
  state <- .np_progress_runtime$fit_state

  if (is.null(state)) {
    return(invisible(NULL))
  }

  if (!is.null(done)) {
    done <- suppressWarnings(as.integer(done)[1L])
    if (is.na(done) || done < 1L) {
      done <- NULL
    }
  }

  .np_progress_runtime$fit_state <- .np_progress_step(
    state = state,
    done = done,
    detail = detail
  )

  invisible(NULL)
}

.np_fit_progress_finish <- function(detail = NULL) {
  state <- .np_progress_runtime$fit_state

  if (is.null(state)) {
    return(invisible(NULL))
  }

  .np_progress_end(state, detail = detail)
  .np_progress_runtime$fit_state <- NULL
  invisible(NULL)
}

.np_fit_progress_abort <- function(detail = NULL) {
  state <- .np_progress_runtime$fit_state

  if (is.null(state)) {
    return(invisible(NULL))
  }

  .np_progress_abort(state, detail = detail)
  .np_progress_runtime$fit_state <- NULL
  invisible(NULL)
}

.np_with_compiled_fit_progress <- function(label, total, handoff = FALSE, handoff.detail = NULL, expr) {
  total <- suppressWarnings(as.integer(total)[1L])
  fit.progress.active <- FALSE
  c.progress.active <- FALSE

  on.exit({
    if (isTRUE(c.progress.active)) {
      .Call("C_np_progress_fit_end", PACKAGE = "npRmpi")
      c.progress.active <- FALSE
    }

    if (isTRUE(fit.progress.active)) {
      .np_fit_progress_abort()
      fit.progress.active <- FALSE
    }
  }, add = TRUE)

  if (isTRUE(.np_progress_enabled(domain = "bandwidth")) &&
      !is.na(total) && total >= 1L) {
    .np_fit_progress_begin(
      label = label,
      total = total,
      handoff = handoff,
      detail = handoff.detail
    )
    fit.progress.active <- TRUE
    .Call("C_np_progress_fit_begin", total, PACKAGE = "npRmpi")
    c.progress.active <- TRUE
  }

  value <- force(expr)

  if (isTRUE(c.progress.active)) {
    .Call("C_np_progress_fit_end", PACKAGE = "npRmpi")
    c.progress.active <- FALSE
  }

  if (isTRUE(fit.progress.active)) {
    .np_fit_progress_finish()
    fit.progress.active <- FALSE
  }

  value
}

.np_densdist_fit_total <- function(bws, tnrow, enrow) {
  if (identical(as.character(bws$type)[1L], "adaptive_nn")) {
    as.integer(tnrow)
  } else {
    as.integer(enrow)
  }
}

.np_condensdist_fit_total <- function(bws, tnrow, enrow) {
  reg.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  base.total <- .np_densdist_fit_total(bws = bws, tnrow = tnrow, enrow = enrow)
  lp.route <- identical(reg.engine, "lp") && isTRUE(bws$xncon > 0L)

  if (isTRUE(lp.route)) {
    as.integer(enrow)
  } else {
    as.integer(2L * base.total)
  }
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

  if (identical(event, "bandwidth_activity_step") && identical(surface, "bandwidth")) {
    .np_progress_bandwidth_activity_step(done = current)
    return(invisible(TRUE))
  }

  if (identical(event, "fit_step") && identical(surface, "bandwidth")) {
    .np_fit_progress_step(done = current)
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
  .np_progress_runtime$bandwidth_context_label <- NULL
  .np_progress_runtime$force_enabled <- FALSE
  invisible(NULL)
}

.np_progress_with_bandwidth_enhanced <- function(expr) {
  old <- getOption("np.progress.bandwidth.enhanced", FALSE)
  on.exit(options(np.progress.bandwidth.enhanced = old), add = TRUE)
  options(np.progress.bandwidth.enhanced = TRUE)
  force(expr)
}

.np_progress_select_bandwidth <- function(label, expr) {
  starting <- !.np_progress_bandwidth_active()
  worker_silent <- FALSE
  if (starting) {
    .np_progress_runtime$bandwidth_old_messages <- getOption("np.messages", TRUE)
    .np_progress_runtime$force_enabled <- isTRUE(.np_progress_runtime$bandwidth_old_messages)
    options(np.messages = FALSE)
    .np_progress_runtime$bandwidth_label <- as.character(label)[1L]

    if (isTRUE(getOption("npRmpi.mpi.initialized", FALSE))) {
      attach.size <- tryCatch(as.integer(mpi.comm.size(1L)), error = function(e) NA_integer_)
      attach.rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) NA_integer_)

      worker_silent <- !is.na(attach.size) && attach.size > 1L &&
        !is.na(attach.rank) && attach.rank != 0L
    }

    if (!worker_silent) {
      .np_progress_runtime$bandwidth_state <- .np_progress_begin(
        label = .np_progress_runtime$bandwidth_label,
        domain = "general",
        surface = "bandwidth"
      )
      if (.np_progress_bandwidth_enhanced_enabled()) {
        .np_progress_runtime$bandwidth_state <- .np_progress_bandwidth_initialize_state(
          .np_progress_runtime$bandwidth_state
        )
        context_label <- .np_progress_runtime$bandwidth_context_label
        if (!is.null(context_label) && nzchar(context_label)) {
          .np_progress_bandwidth_set_coordinator(total_groups = 1L, local_total = 1L)
          .np_progress_bandwidth_set_coordinator_group(1L, context_label)
          .np_progress_runtime$bandwidth_state$start_note_grace_sec <- 0
          .np_progress_runtime$bandwidth_state <- .np_progress_maybe_emit_start_note(
            .np_progress_runtime$bandwidth_state,
            now = .np_progress_runtime$bandwidth_state$started
          )
        }
      }
    }
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

.np_progress_select_bandwidth_enhanced <- function(label, expr) {
  .np_progress_with_bandwidth_enhanced(
    .np_progress_select_bandwidth(label = label, expr = expr)
  )
}

.np_progress_select_iv <- function(label = .np_progress_iv_title(), expr) {
  starting <- !.np_progress_iv_active()
  if (starting) {
    .np_progress_runtime$iv_old_messages <- getOption("np.messages", TRUE)
    .np_progress_runtime$force_enabled <- isTRUE(.np_progress_runtime$iv_old_messages)
    options(np.messages = FALSE)
    .np_progress_runtime$iv_label <- as.character(label)[1L]
    .np_progress_runtime$iv_state <- .np_progress_iv_initialize_state(
      .np_progress_begin(
        label = .np_progress_runtime$iv_label,
        domain = "general",
        surface = "iv_solve"
      )
    )
  }

  .np_progress_runtime$iv_depth <- as.integer(.np_progress_runtime$iv_depth) + 1L
  on.exit({
    .np_progress_runtime$iv_depth <- max(
      0L,
      as.integer(.np_progress_runtime$iv_depth) - 1L
    )
    if (starting) {
      .np_progress_iv_finish()
      options(np.messages = .np_progress_runtime$iv_old_messages)
      .np_progress_iv_clear()
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
