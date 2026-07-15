.np_degree_eval_key <- function(degree) {
  paste(as.integer(degree), collapse = ",")
}

.np_degree_better <- function(candidate, incumbent, direction = c("min", "max")) {
  direction <- match.arg(direction)
  if (!is.finite(candidate))
    return(FALSE)
  if (!is.finite(incumbent))
    return(TRUE)
  if (identical(direction, "min")) candidate < incumbent else candidate > incumbent
}

.np_degree_normalize_bounds <- function(ncon,
                                        degree.min = NULL,
                                        degree.max = NULL,
                                        default.max = 3L) {
  if (!is.numeric(ncon) || length(ncon) != 1L || is.na(ncon) || ncon < 1L)
    stop("automatic degree search requires at least one continuous regressor")

  normalize_one <- function(x, nm, default) {
    if (is.null(x))
      x <- default
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x != floor(x)))
      stop(sprintf("%s must contain finite integers", nm))
    if (length(x) == 1L)
      x <- rep.int(as.integer(x), ncon)
    if (length(x) != ncon)
      stop(sprintf("%s must have length 1 or %d", nm, ncon))
    as.integer(x)
  }

  lower <- normalize_one(degree.min, "degree.min", 0L)
  upper <- normalize_one(degree.max, "degree.max", default.max)

  if (any(lower < 0L))
    stop("degree.min must contain non-negative integers")
  if (any(upper < lower))
    stop("degree.max must be greater than or equal to degree.min coordinate-wise")
  if (any(upper > .np_glp_degree_hard_max))
    stop(sprintf("degree.max must contain integers in [0,%d]", .np_glp_degree_hard_max))

  .np_warn_high_glp_degree(upper, argname = "degree.max")

  candidates <- lapply(seq_len(ncon), function(j) seq.int(lower[j], upper[j]))
  grid.size <- .np_degree_grid_size(candidates)
  singleton <- identical(as.integer(grid.size), 1L)

  list(
    lower = lower,
    upper = upper,
    candidates = candidates,
    grid.size = grid.size,
    singleton = singleton,
    fixed.degree = if (singleton) vapply(candidates, `[`, integer(1L), 1L) else NULL
  )
}

.np_degree_search_engine_controls <- function(search.engine) {
  match.arg(search.engine, c("nomad+powell", "cell", "nomad"))
}

.np_degree_resolve_auto_engine <- function(search.engine,
                                           degree.select,
                                           ncon,
                                           source = "explicit",
                                           auto.filled = character()) {
  search.engine <- .np_degree_search_engine_controls(search.engine)
  degree.select <- match.arg(degree.select, c("manual", "coordinate", "exhaustive"))
  source <- as.character(source)[1L]
  if (is.na(source) || !nzchar(source))
    source <- "explicit"
  auto.filled <- if (is.null(auto.filled)) character() else as.character(auto.filled)

  reason <- "explicit degree-search controls"
  if (identical(source, "auto") &&
      all(c("search.engine", "degree.select") %in% auto.filled)) {
    if (as.integer(ncon[1L]) == 1L) {
      search.engine <- "cell"
      degree.select <- "exhaustive"
      reason <- "auto policy: p=1 LP degree lattice uses exhaustive/cell search"
    } else {
      search.engine <- "nomad+powell"
      degree.select <- "coordinate"
      reason <- "auto policy: p>=2 LP degree search uses NOMAD"
    }
  } else if (identical(source, "auto")) {
    reason <- "auto requested with explicit degree-search controls"
  }

  list(
    search.engine = search.engine,
    degree.select = degree.select,
    source = source,
    reason = reason
  )
}

.np_degree_search_label <- function(method, source = "explicit") {
  method <- as.character(method)[1L]
  source <- as.character(source)[1L]
  nomad.method <- method %in% c("nomad", "nomad+powell")
  base <- if (isTRUE(nomad.method)) {
    "NOMAD degree/bw"
  } else if (identical(method, "exhaustive")) {
    "Exhaustive degree/bw"
  } else {
    "Degree/bw"
  }

  if (identical(source, "auto")) {
    return(if (isTRUE(nomad.method)) {
      "Auto:NOMAD degree/bw"
    } else if (identical(method, "exhaustive")) {
      "Auto:exhaustive degree/bw"
    } else {
      paste0("Auto:", base)
    })
  }

  base
}

.np_degree_search_summary_str <- function(x) {
  ds <- if (is.list(x) && !is.null(x$degree.search)) x$degree.search else NULL
  if (!is.list(ds))
    return("")

  mode <- if (!is.null(ds$mode) && length(ds$mode)) as.character(ds$mode[1L]) else NULL
  engine <- if (!is.null(ds$engine) && length(ds$engine)) as.character(ds$engine[1L]) else mode
  if (is.null(engine) || is.na(engine) || !nzchar(engine))
    return("")

  label <- if (engine %in% c("nomad", "nomad+powell")) {
    "NOMAD"
  } else if (identical(engine, "cell") || identical(mode, "exhaustive")) {
    "Exhaustive"
  } else {
    tools::toTitleCase(engine)
  }

  source <- if (!is.null(ds$source) && length(ds$source)) as.character(ds$source[1L]) else ""
  if (identical(source, "auto"))
    label <- paste0(label, " (auto)")

  paste("\nDegree Search Method:", label)
}

.np_degree_trace_to_frame <- function(records, objective_name = "objective") {
  if (!length(records)) {
    out <- data.frame(
      trace_id = integer(0),
      eval_id = integer(0),
      degree = character(0),
      objective = numeric(0),
      status = character(0),
      cached = logical(0),
      message = character(0),
      elapsed = numeric(0),
      num.feval = numeric(0),
      num.feval.fast = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      trace_id = vapply(records, `[[`, integer(1), "trace_id"),
      eval_id = vapply(records, `[[`, integer(1), "eval_id"),
      degree = vapply(records, function(x) paste(x$degree, collapse = ","), character(1)),
      objective = vapply(records, `[[`, numeric(1), "objective"),
      status = vapply(records, `[[`, character(1), "status"),
      cached = vapply(records, `[[`, logical(1), "cached"),
      message = vapply(records, function(x) {
        if (is.null(x$message)) "" else as.character(x$message)
      }, character(1)),
      elapsed = vapply(records, `[[`, numeric(1), "elapsed"),
      num.feval = vapply(records, function(x) {
        if (is.null(x$num.feval) || length(x$num.feval) != 1L || is.na(x$num.feval))
          NA_real_
        else as.numeric(x$num.feval)
      }, numeric(1)),
      num.feval.fast = vapply(records, function(x) {
        if (is.null(x$num.feval.fast) || length(x$num.feval.fast) != 1L || is.na(x$num.feval.fast))
          NA_real_
        else as.numeric(x$num.feval.fast)
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }

  if (!identical(objective_name, "objective"))
    names(out)[names(out) == "objective"] <- objective_name

  out
}

.np_degree_clip_to_grid <- function(degree, candidates) {
  out <- as.integer(degree)
  for (j in seq_along(candidates)) {
    cand <- as.integer(candidates[[j]])
    out[j] <- cand[which.min(abs(cand - out[j]))]
  }
  out
}

.np_degree_grid_size <- function(candidates) {
  prod(vapply(candidates, length, integer(1)))
}

.np_degree_singleton_search_result <- function(degree.search,
                                               eval_result,
                                               direction = c("min", "max"),
                                               objective_name = "objective") {
  direction <- match.arg(direction)
  fixed.degree <- as.integer(degree.search$fixed.degree)
  objective <- as.numeric(eval_result$objective[1L])
  num.feval <- if (is.null(eval_result$num.feval)) NA_real_ else as.numeric(eval_result$num.feval[1L])
  num.feval.fast <- if (is.null(eval_result$num.feval.fast)) NA_real_ else as.numeric(eval_result$num.feval.fast[1L])
  nn.cache <- if (is.null(eval_result$nn.cache)) NULL else eval_result$nn.cache
  rec <- list(
    degree = fixed.degree,
    objective = objective,
    status = "ok",
    cached = FALSE,
    message = "singleton degree grid; fixed-degree bandwidth search",
    elapsed = NA_real_,
    num.feval = num.feval,
    num.feval.fast = num.feval.fast
  )
  trace <- data.frame(
    trace_id = 1L,
    eval_id = 1L,
    degree = paste(fixed.degree, collapse = ","),
    objective = objective,
    status = "ok",
    cached = FALSE,
    message = rec$message,
    elapsed = NA_real_,
    num.feval = num.feval,
    num.feval.fast = num.feval.fast,
    stringsAsFactors = FALSE
  )
  if (!identical(objective_name, "objective"))
    names(trace)[names(trace) == "objective"] <- objective_name

  list(
    method = degree.search$engine,
    source = if (!is.null(degree.search$source)) degree.search$source else "explicit",
    reason = if (!is.null(degree.search$reason)) degree.search$reason else NULL,
    direction = direction,
    verify = FALSE,
    completed = TRUE,
    certified = TRUE,
    interrupted = FALSE,
    singleton = TRUE,
    fixed.degree = fixed.degree,
    baseline = rec,
    best = rec,
    best_payload = eval_result$payload,
    nn.cache = nn.cache,
    nomad.time = NA_real_,
    powell.time = if (!is.null(eval_result$payload$powell.time)) as.numeric(eval_result$payload$powell.time[1L]) else NA_real_,
    optim.time = if (!is.null(eval_result$payload$total.time)) as.numeric(eval_result$payload$total.time[1L]) else NA_real_,
    n.unique = 1L,
    n.visits = 1L,
    n.cached = 0L,
    grid.size = 1L,
    best.restart = NA_integer_,
    restart.starts = NULL,
    restart.degree.starts = fixed.degree,
    restart.bandwidth.starts = NULL,
    restart.start.info = NULL,
    restart.results = NULL,
    trace = trace
  )
}

.np_degree_search_metadata <- function(search_result,
                                       default_direction = c("min", "max")) {
  default_direction <- match.arg(default_direction)
  scalar_character <- function(value, default = NA_character_) {
    if (is.null(value) || !length(value))
      return(default)
    value <- as.character(value[1L])
    if (is.na(value) || !nzchar(value)) default else value
  }
  scalar_logical <- function(value, default = FALSE) {
    if (is.null(value) || !length(value))
      return(default)
    isTRUE(value[1L])
  }
  scalar_real <- function(value) {
    if (is.null(value) || !length(value))
      return(NA_real_)
    as.numeric(value[1L])
  }
  scalar_integer <- function(value) {
    if (is.null(value) || !length(value))
      return(NA_integer_)
    as.integer(value[1L])
  }
  degree_value <- function(value) {
    if (is.null(value) || !length(value))
      return(integer(0L))
    as.integer(value)
  }

  baseline <- search_result$baseline
  best <- search_result$best
  direction <- scalar_character(search_result$direction, default_direction)

  list(
    mode = scalar_character(search_result$method, "explicit"),
    source = scalar_character(search_result$source, "explicit"),
    reason = if (!is.null(search_result$reason)) search_result$reason else NULL,
    engine = scalar_character(search_result$engine,
                              scalar_character(search_result$method, "explicit")),
    direction = direction,
    verify = scalar_logical(search_result$verify),
    completed = scalar_logical(search_result$completed),
    certified = scalar_logical(search_result$certified),
    interrupted = scalar_logical(search_result$interrupted),
    baseline.degree = degree_value(baseline$degree),
    baseline.fval = scalar_real(baseline$objective),
    best.degree = degree_value(best$degree),
    best.fval = scalar_real(best$objective),
    nomad.time = scalar_real(search_result$nomad.time),
    powell.time = scalar_real(search_result$powell.time),
    optim.time = scalar_real(search_result$optim.time),
    n.unique = scalar_integer(search_result$n.unique),
    n.visits = scalar_integer(search_result$n.visits),
    n.cached = scalar_integer(search_result$n.cached),
    grid.size = scalar_integer(search_result$grid.size),
    singleton = scalar_logical(search_result$singleton),
    fixed.degree = degree_value(search_result$fixed.degree),
    best.restart = scalar_integer(search_result$best.restart),
    restart.starts = search_result$restart.starts,
    restart.degree.starts = search_result$restart.degree.starts,
    restart.bandwidth.starts = search_result$restart.bandwidth.starts,
    restart.start.info = search_result$restart.start.info,
    restart.results = search_result$restart.results,
    nn.cache = search_result$nn.cache,
    trace = search_result$trace
  )
}

.np_degree_format_degree <- function(degree) {
  sprintf("(%s)", paste(as.integer(degree), collapse = ","))
}

.np_degree_format_objective <- function(value) {
  value <- as.numeric(value)[1L]
  if (!is.finite(value))
    return("NA")

  formatC(value, digits = 6L, format = "fg", flag = "#")
}

.np_degree_format_candidate_set <- function(cand) {
  cand <- as.integer(cand)
  if (!length(cand))
    return("c()")

  if (length(cand) >= 2L && all(diff(cand) == 1L))
    return(sprintf("%s:%s", cand[1L], cand[length(cand)]))

  sprintf("c(%s)", paste(cand, collapse = ","))
}

.np_degree_format_search_space <- function(candidates) {
  paste(
    vapply(candidates, .np_degree_format_candidate_set, character(1L)),
    collapse = " x "
  )
}

.np_degree_coordinate_visit_cap <- function(candidates,
                                            max_cycles,
                                            restart_total) {
  restart_total <- npValidatePositiveInteger(restart_total, "restart_total")
  max_cycles <- npValidatePositiveInteger(max_cycles, "degree.max.cycles")
  per_cycle <- sum(vapply(candidates, function(cand) {
    max(0L, length(cand) - 1L)
  }, integer(1L)))
  restart_total * (1L + max_cycles * per_cycle)
}

.np_degree_progress_context <- function() {
  label <- .np_progress_runtime$bandwidth_context_label
  if (is.null(label))
    return(NULL)

  label <- as.character(label)[1L]
  if (is.na(label) || !nzchar(label))
    return(NULL)

  label
}

.np_degree_progress_label <- function(label = NULL) {
  if (!is.null(label) && length(label)) {
    label <- as.character(label)[1L]
    if (!is.na(label) && nzchar(label))
      return(label)
  }
  "Selecting degree and bandwidth"
}

.np_degree_progress_context_fields <- function() {
  context <- .np_degree_progress_context()
  if (is.null(context)) character(0L) else context
}

.np_degree_progress_best_detail <- function(best_record,
                                            objective_name = "objective") {
  if (is.null(best_record))
    return("best pending")

  sprintf(
    "best %s, %s=%s",
    .np_degree_format_degree(best_record$degree),
    objective_name,
    .np_degree_format_objective(best_record$objective)
  )
}

.np_degree_progress_best_degree_detail <- function(best_record) {
  if (is.null(best_record))
    return("best pending")

  sprintf("best %s", .np_degree_format_degree(best_record$degree))
}

.np_degree_progress_detail <- function(phase,
                                       best_record,
                                       objective_name = "objective",
                                       step = NULL,
                                       step_max = NULL,
                                       cycle = NULL,
                                       max_cycles = NULL,
                                       coordinate = NULL,
                                       ncon = NULL,
                                       restart_index = NULL,
                                       restart_total = NULL,
                                       evaluating = NULL) {
  phase <- switch(
    as.character(phase)[1L],
    coordinate = "coord",
    exhaustive = "exhaustive",
    verify = "verify",
    as.character(phase)[1L]
  )
  fields <- c(.np_degree_progress_context_fields(), phase)

  if (!is.null(step) && !is.na(step) && step >= 1L) {
    if (!is.null(step_max) && !is.na(step_max) && step_max >= step) {
      fields <- c(fields, sprintf("step %s/%s", format(step), format(step_max)))
    } else {
      fields <- c(fields, sprintf("step %s", format(step)))
    }
  }

  if (!is.null(restart_total) && !is.na(restart_total) && restart_total > 1L &&
      !is.null(restart_index) && !is.na(restart_index) && restart_index >= 1L) {
    fields <- c(fields, sprintf("restart %s/%s", format(restart_index), format(restart_total)))
  }

  if (!is.null(cycle) && !is.na(cycle) && cycle >= 1L &&
      !is.null(max_cycles) && !is.na(max_cycles) && max_cycles >= 1L) {
    fields <- c(fields, sprintf("cycle %s/%s", format(cycle), format(max_cycles)))
  }

  if (!is.null(coordinate) && !is.na(coordinate) && coordinate >= 1L &&
      !is.null(ncon) && !is.na(ncon) && ncon >= 1L) {
    fields <- c(fields, sprintf("coord %s/%s", format(coordinate), format(ncon)))
  }

  if (!is.null(evaluating)) {
    fields <- c(fields, sprintf("deg %s", .np_degree_format_degree(evaluating)))
  }

  fields <- c(fields, .np_degree_progress_best_detail(
    best_record = best_record,
    objective_name = objective_name
  ))

  paste(fields, collapse = ", ")
}

.np_degree_progress_begin <- function(total = NULL, detail = NULL, label = NULL) {
  state <- .np_progress_begin(
    label = .np_degree_progress_label(label),
    total = total,
    domain = "general",
    surface = "bandwidth"
  )

  .np_progress_show_now(
    state = state,
    done = if (is.null(total)) NULL else 0L,
    detail = detail
  )
}

.np_degree_progress_step <- function(state,
                                     done = NULL,
                                     detail = NULL,
                                     force = FALSE) {
  if (is.null(state))
    return(state)

  if (isTRUE(force)) {
    return(.np_progress_step_at(
      state = state,
      now = .np_progress_now(),
      done = done,
      detail = detail,
      force = TRUE
    ))
  }

  .np_progress_step(
    state = state,
    done = done,
    detail = detail
  )
}

.np_degree_progress_end <- function(state,
                                    detail = NULL,
                                    interrupted = FALSE) {
  if (is.null(state))
    return(invisible(state))

  if (isTRUE(interrupted)) {
    return(invisible(.np_progress_abort(
      state = state,
      detail = paste("Polynomial degree search interrupted;", detail)
    )))
  }

  invisible(.np_progress_end(state = state, detail = detail))
}

.np_degree_check_grid_budget <- function(candidates,
                                         method,
                                         verify = FALSE,
                                         warn_threshold = getOption("np.degree.search.warn.grid", 10000L),
                                         max_threshold = getOption("np.degree.search.max.grid", Inf)) {
  grid.size <- .np_degree_grid_size(candidates)

  if (!is.finite(grid.size) || grid.size < 1)
    stop("searched degree grid is invalid")

  if ((identical(method, "exhaustive") || isTRUE(verify)) &&
      is.finite(max_threshold) &&
      (grid.size > max_threshold)) {
    stop(
      sprintf(
        "exhaustive degree search requires %s degree combinations, exceeding the configured limit of %s",
        format(grid.size, scientific = FALSE, trim = TRUE),
        format(max_threshold, scientific = FALSE, trim = TRUE)
      )
    )
  }

  if ((identical(method, "exhaustive") || isTRUE(verify)) &&
      is.finite(warn_threshold) &&
      (grid.size > warn_threshold)) {
    warning(
      sprintf(
        "exhaustive degree search will evaluate %s degree combinations; this may take considerable time",
        format(grid.size, scientific = FALSE, trim = TRUE)
      ),
      call. = FALSE
    )
  }

  grid.size
}

.np_degree_iterate_grid <- function(candidates, visit, prefix = integer(0), index = 1L) {
  if (index > length(candidates))
    return(isTRUE(visit(as.integer(prefix))))

  for (d in candidates[[index]]) {
    keep_going <- .np_degree_iterate_grid(
      candidates = candidates,
      visit = visit,
      prefix = c(prefix, as.integer(d)),
      index = index + 1L
    )
    if (!isTRUE(keep_going))
      return(FALSE)
  }

  TRUE
}

.np_degree_restart_starts <- function(candidates, restarts, exclude = NULL) {
  restarts <- npValidateNonNegativeInteger(restarts, "degree.restarts")
  if (restarts == 0L)
    return(list())

  ncon <- length(candidates)
  lower <- vapply(candidates, function(cand) as.integer(cand[1L]), integer(1))
  upper <- vapply(candidates, function(cand) as.integer(cand[length(cand)]), integer(1))
  midpoint <- vapply(candidates, function(cand) {
    cand <- as.integer(cand)
    cand[ceiling(length(cand) / 2)]
  }, integer(1))

  seen <- new.env(hash = TRUE, parent = emptyenv())
  starts <- list()

  remember <- function(degree) {
    assign(.np_degree_eval_key(degree), TRUE, envir = seen)
    invisible(NULL)
  }

  already_seen <- function(degree) {
    exists(.np_degree_eval_key(degree), envir = seen, inherits = FALSE)
  }

  if (!is.null(exclude)) {
    if (is.list(exclude)) {
      for (deg in exclude)
        remember(as.integer(deg))
    } else {
      remember(as.integer(exclude))
    }
  }

  add_start <- function(degree) {
    degree <- as.integer(degree)
    if (already_seen(degree))
      return(FALSE)
    remember(degree)
    starts[[length(starts) + 1L]] <<- degree
    TRUE
  }

  seeds <- list(lower, upper, midpoint)
  for (seed in seeds) {
    add_start(seed)
    if (length(starts) >= restarts)
      return(starts)
  }

  k <- 1L
  while (length(starts) < restarts) {
    degree <- vapply(seq_len(ncon), function(j) {
      cand <- as.integer(candidates[[j]])
      idx <- ((k * j) %% length(cand)) + 1L
      cand[idx]
    }, integer(1))
    add_start(degree)
    k <- k + 1L
  }

  starts
}

.np_degree_search_exhaustive <- function(state,
                                         candidates,
                                         phase = "exhaustive",
                                         objective_name = "objective") {
  visited <- 0L

  .np_degree_iterate_grid(candidates, function(degree) {
    if (!is.null(state$progress_state) &&
        identical(state$progress_state$renderer, "single_line")) {
      state$progress_state <- .np_degree_progress_step(
        state = state$progress_state,
        done = visited,
        detail = .np_degree_progress_detail(
          phase = phase,
          best_record = state$best_record,
          objective_name = objective_name,
          evaluating = degree
        ),
        force = TRUE
      )
    }

    state$evaluate(degree)
    visited <<- visited + 1L
    state$progress_state <- .np_degree_progress_step(
      state = state$progress_state,
      done = visited,
      detail = .np_degree_progress_detail(
        phase = phase,
        best_record = state$best_record,
        objective_name = objective_name,
        evaluating = degree
      )
    )
    !isTRUE(state$interrupted)
  })

  invisible(visited)
}

.np_degree_search_coordinate <- function(state,
                                         candidates,
                                         start_degree,
                                         restart_starts,
                                         max_cycles,
                                         direction,
                                         objective_name = "objective") {
  restart_total <- length(restart_starts) + 1L
  step_max <- .np_degree_coordinate_visit_cap(
    candidates = candidates,
    max_cycles = max_cycles,
    restart_total = restart_total
  )

  run_coordinate <- function(initial_degree, restart_index) {
    current <- .np_degree_clip_to_grid(initial_degree, candidates)
    current_record <- state$get_cached(current)
    if (!is.null(state$progress_state) &&
        identical(state$progress_state$renderer, "single_line")) {
      state$progress_state <- .np_degree_progress_step(
        state = state$progress_state,
        done = NULL,
        detail = .np_degree_progress_detail(
          phase = "coordinate",
          best_record = state$best_record,
          objective_name = objective_name,
          step = state$visit_id,
          step_max = step_max,
          restart_index = restart_index,
          restart_total = restart_total,
          evaluating = current
        ),
        force = TRUE
      )
    }
    if (is.null(current_record))
      current_record <- state$evaluate(current)
    state$progress_state <- .np_degree_progress_step(
      state = state$progress_state,
      done = NULL,
      detail = .np_degree_progress_detail(
        phase = "coordinate",
        best_record = state$best_record,
        objective_name = objective_name,
        step = max(1L, state$visit_id - 1L),
        step_max = step_max,
        restart_index = restart_index,
        restart_total = restart_total,
        evaluating = current
      )
    )
    for (cycle in seq_len(max_cycles)) {
      prior <- current
      for (j in seq_along(candidates)) {
        best_local <- current[j]
        best_local_value <- if (!is.null(current_record) &&
                                identical(current_record$status, "ok") &&
                                is.finite(current_record$objective)) {
          current_record$objective
        } else if (identical(direction, "min")) {
          Inf
        } else {
          -Inf
        }
        best_local_record <- current_record
        for (d in candidates[[j]]) {
          if (identical(as.integer(d), current[j]))
            next
          trial <- current
          trial[j] <- as.integer(d)
          if (!is.null(state$progress_state) &&
              identical(state$progress_state$renderer, "single_line")) {
            state$progress_state <- .np_degree_progress_step(
              state = state$progress_state,
              done = NULL,
              detail = .np_degree_progress_detail(
                phase = "coordinate",
                best_record = state$best_record,
                objective_name = objective_name,
                step = state$visit_id,
                step_max = step_max,
                cycle = cycle,
                max_cycles = max_cycles,
                coordinate = j,
                ncon = length(candidates),
                restart_index = restart_index,
                restart_total = restart_total,
                evaluating = trial
              ),
              force = TRUE
            )
          }
          rec <- state$evaluate(trial)
          state$progress_state <- .np_degree_progress_step(
            state = state$progress_state,
            done = NULL,
            detail = .np_degree_progress_detail(
              phase = "coordinate",
              best_record = state$best_record,
              objective_name = objective_name,
              step = max(1L, state$visit_id - 1L),
              step_max = step_max,
              cycle = cycle,
              max_cycles = max_cycles,
              coordinate = j,
              ncon = length(candidates),
              restart_index = restart_index,
              restart_total = restart_total,
              evaluating = trial
            )
          )
          if (isTRUE(state$interrupted))
            return(invisible(NULL))
          if (!identical(rec$status, "ok") || !is.finite(rec$objective))
            next
          if (.np_degree_better(rec$objective, best_local_value, direction = direction)) {
            best_local <- as.integer(d)
            best_local_value <- rec$objective
            best_local_record <- rec
          }
        }
        current[j] <- best_local
        current_record <- best_local_record
      }
      if (identical(current, prior) || isTRUE(state$interrupted))
        break
    }
    invisible(NULL)
  }

  run_coordinate(start_degree, restart_index = 1L)
  if (isTRUE(state$interrupted) || !length(restart_starts))
    return(invisible(NULL))

  for (i in seq_along(restart_starts)) {
    run_coordinate(restart_starts[[i]], restart_index = i + 1L)
    if (isTRUE(state$interrupted))
      break
  }

  invisible(NULL)
}

.np_degree_search <- function(method = c("coordinate", "exhaustive"),
                              candidates,
                              baseline_degree,
                              start_degree,
                              restarts = 0L,
                              max_cycles = 20L,
                              verify = FALSE,
                              eval_fun,
                              direction = c("min", "max"),
                              trace_level = c("full", "none"),
                              source = "explicit",
                              reason = NULL,
                              objective_name = "objective") {
  method <- match.arg(method)
  direction <- match.arg(direction)
  trace_level <- match.arg(trace_level)
  source <- as.character(source)[1L]
  if (is.na(source) || !nzchar(source))
    source <- "explicit"
  progress.label <- .np_degree_search_label(method, source = source)
  restarts <- npValidateNonNegativeInteger(restarts, "degree.restarts")
  max_cycles <- npValidatePositiveInteger(max_cycles, "degree.max.cycles")
  grid.size <- .np_degree_check_grid_budget(
    candidates = candidates,
    method = method,
    verify = verify
  )

  state <- new.env(parent = emptyenv())
  state$cache <- new.env(hash = TRUE, parent = emptyenv())
  state$trace_records <- list()
  state$trace_enabled <- identical(trace_level, "full")
  state$trace_id <- 0L
  state$eval_id <- 0L
  state$visit_id <- 0L
  state$cached_visits <- 0L
  state$nn_cache_stats <- list()
  state$best_record <- NULL
  state$best_payload <- NULL
  state$baseline_record <- NULL
  state$interrupted <- FALSE
  state$progress_state <- NULL

  state$record_trace <- function(rec) {
    if (!isTRUE(state$trace_enabled))
      return(invisible(rec))
    state$trace_id <- state$trace_id + 1L
    rec$trace_id <- state$trace_id
    state$trace_records[[length(state$trace_records) + 1L]] <- rec
    invisible(rec)
  }

  state$update_best <- function(rec, payload = NULL) {
    if (!identical(rec$status, "ok") || !is.finite(rec$objective))
      return(invisible(NULL))
    incumbent <- if (is.null(state$best_record)) {
      if (identical(direction, "min")) Inf else -Inf
    } else {
      state$best_record$objective
    }
    if (.np_degree_better(rec$objective, incumbent, direction = direction)) {
      state$best_record <- rec
      if (!is.null(payload))
        state$best_payload <- payload
    }
    invisible(NULL)
  }

  state$get_cached <- function(degree) {
    key <- .np_degree_eval_key(as.integer(degree))
    if (!exists(key, envir = state$cache, inherits = FALSE))
      return(NULL)
    get(key, envir = state$cache, inherits = FALSE)
  }

  state$evaluate <- function(degree) {
    degree <- as.integer(degree)
    state$visit_id <- state$visit_id + 1L
    key <- .np_degree_eval_key(degree)

    if (exists(key, envir = state$cache, inherits = FALSE)) {
      cached <- get(key, envir = state$cache, inherits = FALSE)
      state$cached_visits <- state$cached_visits + 1L
      cached$cached <- TRUE
      return(cached)
    }

    started <- proc.time()[3]
    payload <- NULL
    status <- "ok"
    msg <- NULL
    objective <- NA_real_
    num.feval <- NA_real_
    num.feval.fast <- NA_real_
    nn.cache <- NULL

    result <- tryCatch(
      {
        .np_progress_bandwidth_set_context(
          sprintf("deg %s", .np_degree_format_degree(degree))
        )
        on.exit(.np_progress_bandwidth_set_context(NULL), add = TRUE)
        eval_fun(degree)
      },
      interrupt = function(e) {
        state$interrupted <- TRUE
        NULL
      },
      error = function(e) {
        status <<- "error"
        msg <<- conditionMessage(e)
        NULL
      }
    )

    if (isTRUE(state$interrupted)) {
      rec <- list(
        eval_id = NA_integer_,
        degree = degree,
        objective = NA_real_,
        status = "interrupt",
        cached = FALSE,
        message = "search interrupted",
        elapsed = proc.time()[3] - started,
        num.feval = NA_real_,
        num.feval.fast = NA_real_
      )
      return(state$record_trace(rec))
    }

    if (!is.null(result)) {
      if (!is.list(result) || is.null(result$objective))
        stop("fixed-degree evaluator must return a list containing 'objective'")
      objective <- as.numeric(result$objective[1L])
      payload <- result$payload
      if (!is.null(result$num.feval))
        num.feval <- as.numeric(result$num.feval[1L])
      if (!is.null(result$num.feval.fast))
        num.feval.fast <- as.numeric(result$num.feval.fast[1L])
      if (!is.null(result$nn.cache))
        nn.cache <- result$nn.cache
    }

    state$eval_id <- state$eval_id + 1L
    rec <- list(
      eval_id = state$eval_id,
      degree = degree,
      objective = objective,
      status = status,
      cached = FALSE,
      message = msg,
      elapsed = proc.time()[3] - started,
      num.feval = num.feval,
      num.feval.fast = num.feval.fast,
      nn.cache = nn.cache
    )
    assign(key, rec, envir = state$cache)
    state$record_trace(rec)
    if (identical(status, "ok")) {
      if (!is.null(nn.cache))
        state$nn_cache_stats[[length(state$nn_cache_stats) + 1L]] <- nn.cache
      state$update_best(rec, payload = payload)
    }
    rec
  }

  baseline_degree <- as.integer(baseline_degree)
  start_degree <- as.integer(start_degree)
  restart_starts <- .np_degree_restart_starts(
    candidates = candidates,
    restarts = restarts,
    exclude = list(start_degree, baseline_degree)
  )

  search_started <- proc.time()[3]

  tryCatch({
    state$baseline_record <- state$evaluate(baseline_degree)

    if (identical(method, "exhaustive")) {
      state$progress_state <- .np_degree_progress_begin(
        total = grid.size,
        label = progress.label,
        detail = .np_degree_progress_detail(
          phase = "exhaustive",
          best_record = state$best_record,
          objective_name = objective_name
        )
      )
      .np_degree_search_exhaustive(
        state = state,
        candidates = candidates,
        phase = "exhaustive",
        objective_name = objective_name
      )
    } else {
      restart_total <- length(restart_starts) + 1L
      state$progress_state <- .np_degree_progress_begin(
        label = progress.label,
        detail = .np_degree_progress_detail(
          phase = "coordinate",
          best_record = state$best_record,
          objective_name = objective_name
        )
      )
      .np_degree_search_coordinate(
        state = state,
        candidates = candidates,
        start_degree = start_degree,
        restart_starts = restart_starts,
        max_cycles = max_cycles,
        direction = direction,
        objective_name = objective_name
      )
      if (!isTRUE(state$interrupted) && isTRUE(verify)) {
        state$progress_state <- .np_degree_progress_end(
          state = state$progress_state,
          detail = .np_degree_progress_detail(
            phase = "coordinate",
            best_record = state$best_record,
            objective_name = objective_name
          )
        )
        state$progress_state <- .np_degree_progress_begin(
          total = grid.size,
          label = progress.label,
          detail = .np_degree_progress_detail(
            phase = "verify",
            best_record = state$best_record,
            objective_name = objective_name
          )
        )
        .np_degree_search_exhaustive(
          state = state,
          candidates = candidates,
          phase = "verify",
          objective_name = objective_name
        )
      }
    }
  }, interrupt = function(e) {
    state$interrupted <- TRUE
    NULL
  })

  state$progress_state <- .np_degree_progress_end(
    state = state$progress_state,
    detail = .np_degree_progress_best_detail(
      best_record = state$best_record,
      objective_name = objective_name
    ),
    interrupted = state$interrupted
  )

  if (is.null(state$best_payload))
    stop("automatic degree search failed to obtain any admissible fitted model")

  search_elapsed <- as.numeric(proc.time()[3] - search_started)

  list(
    method = method,
    source = source,
    reason = reason,
    direction = direction,
    verify = isTRUE(verify),
    completed = !isTRUE(state$interrupted),
    certified = !isTRUE(state$interrupted) && (identical(method, "exhaustive") || isTRUE(verify)),
    interrupted = isTRUE(state$interrupted),
    baseline = state$baseline_record,
    best = state$best_record,
    best_payload = state$best_payload,
    n.unique = state$eval_id,
    n.visits = state$visit_id,
    n.cached = state$cached_visits,
    optim.time = search_elapsed,
    grid.size = grid.size,
    restart.starts = restart_starts,
    nn.cache = .np_r_nn_cache_combine_stats(state$nn_cache_stats),
    trace = .np_degree_trace_to_frame(state$trace_records, objective_name = objective_name)
  )
}

.np_nomad_require_crs <- function(version_fn = utils::packageVersion,
                                  load_namespace = loadNamespace,
                                  minimum_version = "0.15-46") {
  minimum_version_label <- as.character(minimum_version)[1L]
  minimum_version <- package_version(minimum_version_label)

  current_version <- tryCatch(
    version_fn("crs"),
    error = function(e) NULL
  )

  if (is.null(current_version)) {
    stop(
      sprintf(
        paste(
          "automatic degree search with search.engine='nomad' requires the",
          "suggested package 'crs' (>= %s) to provide the NOMAD backend;",
          "install.packages('crs')"
        ),
        minimum_version_label
      ),
      call. = FALSE
    )
  }

  current_version <- package_version(as.character(current_version))

  if (current_version < minimum_version) {
    stop(
      sprintf(
        "automatic degree search with search.engine='nomad' requires 'crs' (>= %s); installed version is %s",
        minimum_version_label,
        as.character(current_version)
      ),
      call. = FALSE
    )
  }

  load_result <- tryCatch(
    {
      suppressPackageStartupMessages(load_namespace("crs"))
      TRUE
    },
    error = function(e) e
  )

  if (inherits(load_result, "error")) {
    stop(
      sprintf(
        "automatic degree search with search.engine='nomad' requires 'crs' (>= %s); failed to load installed version %s: %s",
        minimum_version_label,
        as.character(current_version),
        conditionMessage(load_result)
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.np_degree_extract_random_seed <- function(dots, default = 42L) {
  if (is.null(dots) || length(dots) == 0L)
    return(as.integer(default)[1L])

  if (is.null(names(dots)) || !"random.seed" %in% names(dots))
    return(as.integer(default)[1L])

  npValidateNonNegativeInteger(dots[["random.seed"]], "random.seed")
}

.np_degree_reject_unknown_dots <- function(dots,
                                           where,
                                           allowed = c("random.seed")) {
  if (is.null(dots) || length(dots) == 0L)
    return(invisible(TRUE))

  dot.names <- names(dots)
  if (is.null(dot.names))
    dot.names <- rep("", length(dots))

  bad <- dot.names == "" | !(dot.names %in% allowed)
  if (any(bad))
    .np_reject_unused_dots(dots[bad], where)

  invisible(TRUE)
}

.np_nomad_coerce_start_value <- function(x, type, lb, ub) {
  if (identical(as.integer(type), 3L)) {
    x <- if (is.finite(x) && x >= 0.5) 1 else 0
  } else if (identical(as.integer(type), 1L) || identical(as.integer(type), 2L)) {
    x <- round(x)
  }

  if (is.finite(lb) && x < lb)
    x <- lb
  if (is.finite(ub) && x > ub)
    x <- ub

  as.numeric(x)
}

.np_nomad_complete_start_point <- function(point,
                                           lower,
                                           upper,
                                           ncont = 0L) {
  n <- length(lower)
  out <- rep(NA_real_, n)

  if (!is.null(point)) {
    point <- as.numeric(point)
    if (length(point) == n)
      out <- point
  }

  ncont <- max(0L, min(as.integer(ncont), n))
  if (ncont > 0L) {
    cont.default <- pmax(lower[seq_len(ncont)], pmin(upper[seq_len(ncont)], 1.5))
    bad <- !is.finite(out[seq_len(ncont)]) | out[seq_len(ncont)] <= 0
    out[seq_len(ncont)][bad] <- cont.default[bad]
  }

  if (n > ncont) {
    cat.idx <- (ncont + 1L):n
    cat.default <- pmax(lower[cat.idx], pmin(upper[cat.idx], 0.5 * upper[cat.idx]))
    bad <- !is.finite(out[cat.idx]) | out[cat.idx] < lower[cat.idx] | out[cat.idx] > upper[cat.idx]
    out[cat.idx][bad] <- cat.default[bad]
  }

  pmax(lower, pmin(upper, out))
}

.np_lp_nomad_degree_start_policy <- function(value = getOption("np.nomad.degree.start.policy", "low_first_full_random")) {
  value <- as.character(value)[1L]
  choices <- c(
    "low_first_full_random",
    "mid_first_full_random",
    "anchor_then_random",
    "spread_then_random",
    "random_full_only"
  )
  if (is.na(value) || !nzchar(value) || !(value %in% choices)) {
    stop(
      "option 'np.nomad.degree.start.policy' must be one of: ",
      paste(choices, collapse = ", "),
      call. = FALSE
    )
  }
  value
}

.np_lp_nomad_build_degree_starts <- function(initial,
                                             lower,
                                             upper,
                                             basis = c("glp", "additive", "tensor"),
                                             nobs,
                                             nmulti = 1L,
                                             random.seed = 42L,
                                             max.tries = 256L,
                                             user_supplied = FALSE) {
  basis <- match.arg(basis)
  lower <- as.integer(lower)
  upper <- as.integer(upper)
  q <- length(lower)
  nstart <- npValidateNmulti(nmulti)
  policy <- .np_lp_nomad_degree_start_policy()

  if (!q)
    return(matrix(integer(0), nrow = nstart, ncol = 0L))
  if (any(upper < lower))
    stop("invalid degree bounds in NOMAD degree-start policy", call. = FALSE)

  starts <- matrix(NA_integer_, nrow = nstart, ncol = q)
  midpoint <- as.integer(round((lower + upper) / 2))

  draw_full <- function(nrow) {
    if (nrow <= 0L)
      return(matrix(integer(0), nrow = 0L, ncol = q))
    t(vapply(seq_len(nrow), function(j) {
      vapply(seq_len(q), function(i) {
        sample.int(upper[i] - lower[i] + 1L, 1L) + lower[i] - 1L
      }, integer(1L))
    }, integer(q)))
  }

  add_unique <- function(pool, row) {
    row <- as.integer(row)
    if (!nrow(pool))
      return(matrix(row, nrow = 1L))
    key <- paste(row, collapse = ",")
    keys <- apply(pool, 1L, paste, collapse = ",")
    if (key %in% keys)
      pool
    else
      rbind(pool, row)
  }

  anchor_pool <- function() {
    if (nstart <= 1L)
      return(matrix(lower, nrow = 1L))
    if (nstart == 2L)
      return(rbind(lower, upper))
    rbind(lower, midpoint, upper)
  }

  spread_pool <- function() {
    pool <- matrix(integer(0), nrow = 0L, ncol = q)
    pool <- add_unique(pool, lower)
    pool <- add_unique(pool, midpoint)
    pool <- add_unique(pool, upper)
    if (q > 1L) {
      level.mat <- rbind(lower, midpoint, upper)
      base.pattern <- c(3L, 2L, 1L)
      for (j in seq_len(max(q, 3L))) {
        idx <- base.pattern[((seq_len(q) + j - 2L) %% 3L) + 1L]
        row <- vapply(seq_len(q), function(i) level.mat[idx[i], i], integer(1L))
        pool <- add_unique(pool, row)
      }
    }
    pool
  }

  initial <- as.integer(pmax(lower, pmin(upper, initial)))

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

  if (isTRUE(user_supplied)) {
    starts[1L, ] <- initial
    fill.from <- 2L
  } else if (identical(policy, "low_first_full_random")) {
    starts[1L, ] <- lower
    fill.from <- 2L
  } else if (identical(policy, "mid_first_full_random")) {
    starts[1L, ] <- midpoint
    fill.from <- 2L
  } else if (identical(policy, "anchor_then_random")) {
    anchors <- anchor_pool()
    take <- min(nstart, nrow(anchors))
    starts[seq_len(take), ] <- anchors[seq_len(take), , drop = FALSE]
    fill.from <- take + 1L
  } else if (identical(policy, "spread_then_random")) {
    pool <- if (nstart <= 2L) anchor_pool() else spread_pool()
    take <- min(nstart, nrow(pool))
    starts[seq_len(take), ] <- pool[seq_len(take), , drop = FALSE]
    fill.from <- take + 1L
  } else if (identical(policy, "random_full_only")) {
    fill.from <- 1L
  } else {
    stop(sprintf("unknown np.nomad.degree.start.policy: %s", policy), call. = FALSE)
  }

  if (fill.from <= nstart)
    starts[fill.from:nstart, ] <- draw_full(nstart - fill.from + 1L)

  if (any(is.na(starts)) ||
      any(t(t(starts) < lower)) ||
      any(t(t(starts) > upper)))
    stop("generated NOMAD degree start outside bounds", call. = FALSE)

  storage.mode(starts) <- "integer"
  starts
}

.np_nomad_build_starts <- function(x0,
                                   bbin,
                                   lb,
                                   ub,
                                   nmulti = 1L,
                                   random.seed = 42L,
                                   degree_spec = NULL,
                                   start.lower = NULL,
                                   start.upper = NULL) {
  n <- length(lb)
  if (length(ub) != n || length(bbin) != n)
    stop("NOMAD start construction requires matching lengths for x0/bbin/lb/ub")

  nstart <- npValidateNmulti(nmulti)
  starts <- matrix(0, nrow = nstart, ncol = n)
  normalize.start.bound <- function(value, fallback, name) {
    if (is.null(value))
      return(as.numeric(fallback))
    value <- as.numeric(value)
    if (length(value) == 1L)
      value <- rep.int(value, n)
    if (length(value) != n)
      stop(sprintf("NOMAD start construction requires '%s' to have length 1 or %d", name, n),
           call. = FALSE)
    value
  }
  start.lower <- normalize.start.bound(start.lower, lb, "start.lower")
  start.upper <- normalize.start.bound(start.upper, ub, "start.upper")

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

  for (j in seq_len(nstart)) {
    for (i in seq_len(n)) {
      lo <- if (is.finite(start.lower[i])) start.lower[i] else if (is.finite(lb[i])) lb[i] else -1
      hi <- if (is.finite(start.upper[i])) start.upper[i] else if (is.finite(ub[i])) ub[i] else 1
      if (hi < lo) {
        tmp <- lo
        lo <- hi
        hi <- tmp
      }
      starts[j, i] <- .np_nomad_coerce_start_value(
        runif(1L, min = lo, max = hi),
        type = bbin[i],
        lb = lb[i],
        ub = ub[i]
      )
    }
  }

  if (!is.null(degree_spec)) {
    q <- length(degree_spec$lower)
    if (q > 0L) {
      degree.idx <- (n - q + 1L):n
      degree.starts <- .np_lp_nomad_build_degree_starts(
        initial = if (is.null(degree_spec$initial)) rep.int(1L, q) else degree_spec$initial,
        lower = degree_spec$lower,
        upper = degree_spec$upper,
        basis = if (is.null(degree_spec$basis)) "glp" else degree_spec$basis,
        nobs = degree_spec$nobs,
        nmulti = nstart,
        random.seed = random.seed,
        user_supplied = isTRUE(degree_spec$user_supplied)
      )
      starts[, degree.idx] <- degree.starts
    }
  }

  if (is.null(x0))
    return(starts)

  x0 <- as.numeric(x0)
  if (length(x0) == n) {
    starts[1L, ] <- vapply(
      seq_len(n),
      function(i) .np_nomad_coerce_start_value(x0[i], type = bbin[i], lb = lb[i], ub = ub[i]),
      numeric(1L)
    )
    if (!is.null(degree_spec)) {
      q <- length(degree_spec$lower)
      if (q > 0L) {
        degree.idx <- (n - q + 1L):n
        starts[1L, degree.idx] <- .np_lp_nomad_build_degree_starts(
          initial = if (is.null(degree_spec$initial)) rep.int(1L, q) else degree_spec$initial,
          lower = degree_spec$lower,
          upper = degree_spec$upper,
          basis = if (is.null(degree_spec$basis)) "glp" else degree_spec$basis,
          nobs = degree_spec$nobs,
          nmulti = 1L,
          random.seed = random.seed,
          user_supplied = isTRUE(degree_spec$user_supplied)
        )[1L, ]
      }
    }
    return(starts)
  }

  if (length(x0) >= n * nstart) {
    for (j in seq_len(nstart)) {
      starts[j, ] <- vapply(
        seq_len(n),
        function(i) .np_nomad_coerce_start_value(
          x0[j + (i - 1L) * nstart],
          type = bbin[i],
          lb = lb[i],
          ub = ub[i]
        ),
        numeric(1L)
      )
    }
  }

  starts
}

.np_nomad_normalize_user_opts <- function(nomad.opts, where = "NOMAD") {
  if (is.null(nomad.opts))
    return(list())

  if (!is.list(nomad.opts))
    stop(sprintf("%s: 'nomad.opts' must be a list", where), call. = FALSE)

  nomad.opts
}

.np_nomad_default_opts <- function(random.seed, nomad.opts = list()) {
  base <- list(
    SEED = as.integer(random.seed),
    RNG_ALT_SEEDING = TRUE,
    DIRECTION_TYPE = "ORTHO 2N",
    QUAD_MODEL_SEARCH = "no",
    NM_SEARCH = "no",
    SPECULATIVE_SEARCH = "no",
    EVAL_OPPORTUNISTIC = "no"
  )
  utils::modifyList(base, if (is.null(nomad.opts)) list() else nomad.opts)
}

.np_nomad_prepare_solver_opts <- function(random.seed,
                                          nomad.opts = list(),
                                          coordinate.roles = NULL,
                                          expected.length = NULL,
                                          geometry.policy = c("user-only",
                                                              "generate-central"),
                                          where = "NOMAD source geometry") {
  geometry.policy <- match.arg(geometry.policy)
  opts <- .np_nomad_default_opts(random.seed, nomad.opts)
  if (!identical(geometry.policy, "generate-central"))
    return(opts)

  .np_nomad_apply_source_geometry(
    opts,
    user.opts = nomad.opts,
    roles = coordinate.roles,
    expected.length = expected.length,
    where = where
  )
}

.np_nomad_native_false_option <- function(value) {
  if (is.logical(value))
    return(length(value) >= 1L && !isTRUE(value[1L]))
  token <- tolower(trimws(as.character(value[1L])))
  token %in% c("false", "f", "no", "n", "0")
}

.np_nomad_native_reject_unsupported_options <- function(opts, route) {
  if (is.null(opts) || !length(opts))
    return(invisible(TRUE))
  option.names <- names(opts)
  if (is.null(option.names))
    return(invisible(TRUE))
  idx <- which(toupper(trimws(option.names)) == "EVAL_USE_CACHE")
  if (length(idx) && any(vapply(opts[idx], .np_nomad_native_false_option, logical(1)))) {
    stop(sprintf(
      "%s requires EVAL_USE_CACHE = TRUE; cache-off native solves are not supported",
      route
    ), call. = FALSE)
  }
  invisible(TRUE)
}

.np_nomad_native_route_uses_solver <- function(nomad = FALSE,
                                               degree.select = NULL,
                                               search.engine = NULL,
                                               bwsolver = NULL) {
  match_one <- function(value, choices, default) {
    if (is.null(value))
      return(default)
    out <- tryCatch(match.arg(value, choices), error = function(e) NULL)
    if (is.null(out) || !length(out))
      return(default)
    as.character(out[1L])
  }

  degree.value <- match_one(degree.select, c("manual", "coordinate", "exhaustive"), "manual")
  engine.value <- match_one(search.engine, c("nomad+powell", "cell", "nomad"), "nomad+powell")
  solver.value <- match_one(bwsolver, c("powell", "mads", "mads+powell"), "powell")

  nomad.value <- tryCatch(
    npValidateNomadControl(nomad, "nomad"),
    error = function(e) "false"
  )

  nomad.value %in% c("true", "auto") ||
    (!identical(degree.value, "manual") && engine.value %in% c("nomad", "nomad+powell")) ||
    solver.value %in% c("mads", "mads+powell")
}

.np_nomad_native_reject_unsupported_options_for_route <- function(opts,
                                                                  route,
                                                                  nomad = FALSE,
                                                                  degree.select = NULL,
                                                                  search.engine = NULL,
                                                                  bwsolver = NULL) {
  if (.np_nomad_native_route_uses_solver(
    nomad = nomad,
    degree.select = degree.select,
    search.engine = search.engine,
    bwsolver = bwsolver
  )) {
    .np_nomad_native_reject_unsupported_options(opts, route)
  }
  invisible(TRUE)
}

.np_nomad_native_reject_unsupported_options_from_dots <- function(dots, route) {
  if (is.null(dots) || !length(dots) || !("nomad.opts" %in% names(dots)))
    return(invisible(TRUE))
  .np_nomad_native_reject_unsupported_options_for_route(
    opts = dots$nomad.opts,
    route = route,
    nomad = if ("nomad" %in% names(dots)) dots$nomad else FALSE,
    degree.select = if ("degree.select" %in% names(dots)) dots$degree.select else NULL,
    search.engine = if ("search.engine" %in% names(dots)) dots$search.engine else NULL,
    bwsolver = if ("bwsolver" %in% names(dots)) dots$bwsolver else NULL
  )
}

.np_nomad_native_option_vectors <- function(opts) {
  if (is.null(opts) || !length(opts))
    return(list(names = character(), values = character()))

  .np_nomad_native_reject_unsupported_options(opts, "native NOMAD route")

  option.names <- names(opts)
  if (is.null(option.names) || any(!nzchar(option.names)))
    stop("native NOMAD route received unnamed NOMAD options", call. = FALSE)

  option.values <- vapply(opts, function(value) {
    if (is.logical(value)) {
      if (isTRUE(value[1L])) "true" else "false"
    } else if (length(value) > 1L) {
      paste0("( ", paste(as.character(value), collapse = " "), " )")
    } else {
      as.character(value[1L])
    }
  }, character(1L))

  list(names = as.character(option.names), values = option.values)
}

.np_native_nomad_field <- function(x, name, default = NA) {
  if (!is.list(x) || is.null(names(x)) || !(name %in% names(x)))
    return(default)
  x[[name, exact = TRUE]]
}

.np_native_nomad_cache_diagnostics <- function(native) {
  list(
    cache.hits = as.integer(.np_native_nomad_field(native, "cache_hits", NA_integer_)[1L]),
    cache.size = as.integer(.np_native_nomad_field(native, "cache_size", NA_integer_)[1L]),
    total.evaluations = as.integer(.np_native_nomad_field(native, "total_evaluations", NA_integer_)[1L])
  )
}

.np_nomad_native_r_callback_search <- function(eval.f,
                                               x0,
                                               bbin,
                                               lb,
                                               ub,
                                               random.seed = 42L,
                                               inner.start.count = 0L,
                                               option.names = character(),
                                               option.values = character(),
                                               display.nomad.progress = FALSE) {
  if (!is.function(eval.f))
    stop("native NOMAD R callback route requires a function", call. = FALSE)

  native.call <- .np_nomad_capture_solver_output(.Call(
    "C_np_nomad_r_callback_native_search",
    eval.f,
    environment(eval.f),
    as.double(x0),
    as.integer(bbin),
    as.double(lb),
    as.double(ub),
    as.integer(random.seed),
    as.integer(inner.start.count),
    as.character(option.names),
    as.character(option.values),
    as.integer(if (isTRUE(display.nomad.progress)) 0L else 1L),
    PACKAGE = "np"
  ), capture.output = !isTRUE(display.nomad.progress))

  native.call
}

.np_nomad_progress_detail <- function(current_degree,
                                      best_record,
                                      iteration = NULL,
                                      cumulative_iteration = NULL,
                                      restart_index = NULL,
                                      nmulti = 1L,
                                      restart_durations = numeric(),
                                      elapsed = NULL) {
  fields <- .np_degree_progress_context_fields()
  nmulti <- suppressWarnings(as.integer(nmulti)[1L])
  restart_index <- suppressWarnings(as.integer(restart_index)[1L])
  if (!is.na(nmulti) && nmulti > 1L) {
    if (!is.na(restart_index) && restart_index >= 1L) {
      fields <- c(fields, sprintf("multistart %s/%s", format(restart_index), format(nmulti)))
    } else {
      fields <- c(fields, sprintf("multistarts %s", format(nmulti)))
    }
  }

  cumulative_iteration <- suppressWarnings(as.integer(cumulative_iteration)[1L])
  if (!is.null(iteration) && is.finite(iteration) && iteration >= 1L) {
    if (!is.na(cumulative_iteration) &&
        cumulative_iteration >= 1L &&
        cumulative_iteration != iteration) {
      fields <- c(
        fields,
        sprintf("iteration %s (%s)", format(iteration), format(cumulative_iteration))
      )
    } else {
      fields <- c(fields, sprintf("iteration %s", format(iteration)))
    }
  }

  if (is.finite(elapsed) && !is.na(elapsed) && elapsed >= 0)
    fields <- c(fields, sprintf("elapsed %ss", .np_progress_fmt_num(elapsed)))

  if (!is.null(current_degree))
    fields <- c(fields, sprintf("deg %s", .np_degree_format_degree(current_degree)))

  fields <- c(fields, .np_degree_progress_best_degree_detail(best_record = best_record))

  paste(fields, collapse = ", ")
}

.np_nomad_progress_start_detail <- function(baseline_degree,
                                            restart_index = 1L,
                                            nmulti,
                                            best_record,
                                            restart_durations = numeric(),
                                            elapsed = NULL) {
  .np_nomad_progress_detail(
    current_degree = baseline_degree,
    best_record = best_record,
    restart_index = restart_index,
    nmulti = nmulti,
    restart_durations = restart_durations,
    elapsed = elapsed
  )
}

.np_nomad_progress_fields <- function(state,
                                      done = NULL,
                                      detail = NULL,
                                      now = .np_progress_now()) {
  .np_nomad_progress_detail(
    current_degree = state$nomad_current_degree,
    best_record = state$nomad_best_record,
    iteration = done,
    cumulative_iteration = state$nomad_eval_id,
    restart_index = state$nomad_restart_index,
    nmulti = state$nomad_nmulti,
    restart_durations = state$nomad_restart_durations,
    elapsed = max(0, now - state$started)
  )
}

.np_nomad_progress_begin <- function(nmulti,
                                     baseline_degree,
                                     best_record,
                                     label = NULL) {
  state <- .np_progress_begin(
    label = .np_degree_progress_label(label),
    domain = "general",
    surface = "bandwidth"
  )

  state$unknown_total_fields <- .np_nomad_progress_fields
  state$nomad_nmulti <- npValidateNmulti(nmulti)
  state$nomad_restart_index <- 1L
  state$nomad_restart_durations <- numeric()
  state$nomad_current_degree <- as.integer(baseline_degree)
  state$nomad_best_record <- best_record
  state$nomad_eval_id <- 0L
  state$start_note <- if (state$nomad_nmulti > 1L) {
    sprintf(
      "%s %s (multistart 1/%s)",
      state$pkg_prefix,
      state$label,
      format(state$nomad_nmulti)
    )
  } else {
    sprintf("%s %s...", state$pkg_prefix, state$label)
  }

  .np_progress_show_now(state)
}

.np_nomad_progress_configure <- function(state,
                                         nmulti,
                                         baseline_degree,
                                         best_record,
                                         label = NULL) {
  if (is.null(state))
    return(state)

  state$label <- .np_degree_progress_label(label)
  state$unknown_total_fields <- .np_nomad_progress_fields
  state$nomad_nmulti <- npValidateNmulti(nmulti)
  state$nomad_restart_index <- 1L
  state$nomad_restart_durations <- numeric()
  state$nomad_current_degree <- as.integer(baseline_degree)
  state$nomad_best_record <- best_record
  state$nomad_eval_id <- 0L
  state$start_note <- if (state$nomad_nmulti > 1L) {
    sprintf(
      "%s %s (multistart 1/%s)",
      state$pkg_prefix,
      state$label,
      format(state$nomad_nmulti)
    )
  } else {
    sprintf("%s %s...", state$pkg_prefix, state$label)
  }

  .np_progress_show_now(state)
}

.np_nomad_baseline_note <- function(degree) {
  invisible(NULL)
}

.np_nomad_powell_note <- function(degree) {
  invisible(NULL)
}

.np_nomad_powell_progress_label <- function() {
  "Refining bandwidth"
}

.np_nomad_powell_hotstart_nmulti <- function(strategy = c("disable_multistart",
                                                          "single_iteration")) {
  strategy <- match.arg(strategy)
  switch(strategy,
         disable_multistart = 1L,
         single_iteration = 1L)
}

.np_nomad_powell_hotstart_opt_args <- function(opt.args,
                                               strategy = c("disable_multistart",
                                                            "single_iteration"),
                                               remin = FALSE) {
  out <- opt.args
  out$nmulti <- .np_nomad_powell_hotstart_nmulti(strategy)
  out$powell.remin <- isTRUE(remin)
  out$bwsolver <- "powell"
  out
}

.np_nomad_powell_context_label <- function(degree) {
  sprintf(
    "Refining NOMAD solution with one Powell hot start at degree %s",
    .np_degree_format_degree(degree)
  )
}

.np_nomad_powell_progress_detail <- function(current_degree,
                                             best_record,
                                             iteration = NULL,
                                             elapsed = NULL) {
  fields <- .np_degree_progress_context_fields()

  if (is.finite(elapsed) && !is.na(elapsed) && elapsed >= 0)
    fields <- c(fields, sprintf("elapsed %ss", .np_progress_fmt_num(elapsed)))

  if (!is.null(iteration)) {
    iteration <- suppressWarnings(as.integer(iteration)[1L])
    if (!is.na(iteration) && iteration >= 1L)
      fields <- c(fields, sprintf("iter %s", format(iteration)))
  }

  if (!is.null(current_degree))
    fields <- c(fields, sprintf("deg %s", .np_degree_format_degree(current_degree)))

  fields <- c(fields, .np_degree_progress_best_degree_detail(best_record = best_record))

  paste(fields, collapse = ", ")
}

.np_nomad_powell_progress_fields <- function(state,
                                             done = NULL,
                                             detail = NULL,
                                             now = .np_progress_now()) {
  .np_nomad_powell_progress_detail(
    current_degree = state$nomad_current_degree,
    best_record = state$nomad_best_record,
    iteration = done,
    elapsed = max(0, now - state$started)
  )
}

.np_nomad_progress_enter_powell <- function(state,
                                            degree,
                                            best_record) {
  if (is.null(state))
    return(state)

  now <- .np_progress_now()
  state$label <- .np_nomad_powell_progress_label()
  state$unknown_total_fields <- .np_nomad_powell_progress_fields
  state$nomad_current_degree <- as.integer(degree)
  state$nomad_best_record <- best_record
  state$started <- now
  state$last_done <- NULL
  state$last_emitted_done <- NULL
  state$last_emitted_detail <- NULL
  state$last_emit <- now - state$throttle_sec
  .np_progress_step_at(
    state = state,
    now = now,
    done = NULL,
    force = TRUE
  )
}

.np_nomad_with_powell_progress <- function(degree, expr, best_record = NULL) {
  old.context <- .np_progress_runtime$bandwidth_context_label
  old.state <- .np_progress_runtime$bandwidth_state
  local.state <- NULL
  progress.enabled <- isTRUE(.np_progress_enabled(domain = "general"))

  powell.context <- if (!is.null(old.context) && nzchar(old.context)) old.context else NULL
  .np_progress_bandwidth_set_context(powell.context)
  on.exit({
    if (!is.null(local.state) && !is.null(.np_progress_runtime$bandwidth_state)) {
      .np_progress_end(.np_progress_runtime$bandwidth_state)
    }
    if (!is.null(local.state)) {
      .np_progress_runtime$bandwidth_state <- old.state
    }
    .np_progress_bandwidth_set_context(old.context)
  }, add = TRUE)

  if (is.null(old.state) && isTRUE(progress.enabled)) {
    local.state <- .np_progress_begin(
      label = .np_nomad_powell_progress_label(),
      domain = "general",
      surface = "bandwidth"
    )
    local.state$unknown_total_fields <- .np_nomad_powell_progress_fields
    local.state$nomad_current_degree <- as.integer(degree)
    local.state$nomad_best_record <- best_record
    local.state$nomad_nmulti <- 1L
    local.state <- .np_progress_show_now(local.state)
    .np_progress_runtime$bandwidth_state <- local.state
  }

  value <- force(expr)

  if (!is.null(local.state) && is.list(value) && !is.null(value$num.feval)) {
    done <- suppressWarnings(as.integer(value$num.feval[1L]))
    if (!is.na(done) && done >= 1L && !is.null(.np_progress_runtime$bandwidth_state)) {
      .np_progress_runtime$bandwidth_state <- .np_progress_step_at(
        state = .np_progress_runtime$bandwidth_state,
        now = .np_progress_now(),
        done = done,
        force = TRUE
      )
    }
  }

  value
}

.np_nomad_capture_solver_output <- function(expr, capture.output = TRUE) {
  if (!isTRUE(capture.output)) {
    return(list(value = force(expr), output = character()))
  }

  value <- NULL
  output <- utils::capture.output(
    value <- force(expr),
    type = "output"
  )

  list(value = value, output = output)
}

.np_nomad_external_output_lines <- function(lines) {
  lines <- as.character(lines)
  lines <- trimws(lines)
  lines <- lines[!is.na(lines) & nzchar(lines)]
  unique(lines)
}

.np_nomad_render_external_output <- function(line, progress_state) {
  if (is.null(progress_state) ||
      !isTRUE(progress_state$enabled) ||
      !isTRUE(progress_state$visible) ||
      !identical(progress_state$renderer, "single_line")) {
    return(FALSE)
  }

  line <- .np_io_prefix_text(paste0("NOMAD: ", line))
  render_line <- .np_progress_ellipsize_middle(
    line,
    max_width = max(20L, getOption("width", 80L) - 1L)
  )
  snapshot <- .np_progress_make_snapshot(
    state = progress_state,
    line = line,
    render_line = render_line,
    event = "render",
    now = .np_progress_now(),
    done = progress_state$last_done,
    detail = progress_state$last_emitted_detail
  )
  .np_progress_render_single_line(snapshot = snapshot, event = "render")

  transient_state <- progress_state
  transient_state$rendered <- TRUE
  transient_state$last_render_width <- nchar(render_line, type = "width")
  transient_state$last_line <- line
  .np_progress_prepare_for_external_output(transient_state)
  TRUE
}

.np_nomad_emit_external_output <- function(lines, progress_state, seen) {
  lines <- .np_nomad_external_output_lines(lines)
  lines <- lines[!(lines %in% seen)]
  if (!length(lines) || !isTRUE(.np_progress_enabled())) {
    return(list(progress_state = progress_state, seen = seen))
  }

  progress_state <- .np_progress_prepare_for_external_output(progress_state)
  for (line in lines) {
    if (!isTRUE(.np_nomad_render_external_output(line, progress_state))) {
      .np_message("NOMAD: ", line)
    }
  }

  list(progress_state = progress_state, seen = c(seen, lines))
}

.np_nomad_native_call_value <- function(native.call, seen = character()) {
  progress_state <- .np_progress_runtime$bandwidth_state
  if (!length(seen) && !is.null(progress_state$external_output_seen)) {
    seen <- progress_state$external_output_seen
  }
  emitted <- .np_nomad_emit_external_output(
    lines = native.call$output,
    progress_state = progress_state,
    seen = seen
  )
  if (!is.null(emitted$progress_state)) {
    emitted$progress_state$external_output_seen <- emitted$seen
  }
  .np_progress_runtime$bandwidth_state <- emitted$progress_state
  native.call$value
}

.np_nomad_native_progress_begin <- function(nmulti,
                                            baseline_degree,
                                            best_record,
                                            label = NULL) {
  handle <- new.env(parent = emptyenv())
  handle$old_state <- .np_progress_runtime$bandwidth_state
  handle$closed <- FALSE
  handle$state <- .np_nomad_progress_begin(
    nmulti = nmulti,
    baseline_degree = baseline_degree,
    best_record = best_record,
    label = label
  )
  handle$state$nomad_native_progress <- TRUE
  handle$state$nomad_eval_offset <- 0L
  .np_progress_runtime$bandwidth_state <- handle$state
  handle
}

.np_nomad_native_progress_state <- function(handle) {
  if (is.null(handle) || isTRUE(handle$closed)) {
    return(NULL)
  }

  state <- .np_progress_runtime$bandwidth_state
  if (is.null(state)) {
    state <- handle$state
  }
  state
}

.np_nomad_progress_best_record <- function(incumbent,
                                           candidate,
                                           direction = c("min", "max")) {
  direction <- match.arg(direction)
  if (is.null(candidate))
    return(incumbent)
  if (is.null(incumbent))
    return(candidate)

  candidate.objective <- suppressWarnings(as.numeric(candidate$objective)[1L])
  incumbent.objective <- suppressWarnings(as.numeric(incumbent$objective)[1L])
  if (.np_degree_better(candidate.objective, incumbent.objective, direction = direction))
    return(candidate)

  incumbent
}

.np_progress_nomad_native_step_from_c <- function(iteration,
                                                  current.degree,
                                                  best.degree = integer(),
                                                  best.objective = NA_real_,
                                                  force = FALSE) {
  state <- .np_progress_runtime$bandwidth_state
  if (is.null(state) || !isTRUE(state$nomad_native_progress)) {
    return(invisible(FALSE))
  }

  iteration <- suppressWarnings(as.integer(iteration)[1L])
  if (is.na(iteration) || iteration < 1L) {
    return(invisible(FALSE))
  }

  current.degree <- as.integer(current.degree)
  if (length(current.degree) && !anyNA(current.degree)) {
    state$nomad_current_degree <- current.degree
  }

  eval.offset <- suppressWarnings(as.integer(state$nomad_eval_offset)[1L])
  if (is.na(eval.offset) || eval.offset < 0L) {
    eval.offset <- 0L
  }
  state$nomad_eval_id <- eval.offset + iteration

  best.degree <- as.integer(best.degree)
  if (length(best.degree) && !anyNA(best.degree)) {
    best.objective <- suppressWarnings(as.numeric(best.objective)[1L])
    c_best_record <- list(
      eval_id = state$nomad_eval_id,
      degree = best.degree,
      objective = best.objective,
      status = "ok",
      cached = FALSE
    )
    state$nomad_best_record <- .np_nomad_progress_best_record(
      incumbent = state$nomad_best_record,
      candidate = c_best_record,
      direction = "min"
    )
  }

  state <- .np_degree_progress_step(
    state = state,
    done = iteration,
    detail = NULL,
    force = isTRUE(force)
  )
  .np_progress_runtime$bandwidth_state <- state
  invisible(TRUE)
}

.np_progress_nomad_native_observer_config <- function() {
  state <- .np_progress_runtime$bandwidth_state
  interval <- if (!is.null(state$throttle_sec)) {
    suppressWarnings(as.numeric(state$throttle_sec)[1L])
  } else {
    .np_progress_interval_sec(known_total = FALSE, domain = "general")
  }
  if (!is.finite(interval) || is.na(interval) || interval < 0) {
    interval <- .np_progress_interval_sec(known_total = FALSE, domain = "general")
  }

  c(
    enabled = as.numeric(!is.null(state) && isTRUE(state$enabled)),
    interval = interval
  )
}

.np_progress_nomad_native_observer_dispatch <- function(iteration,
                                                         current.degree,
                                                         best.degree,
                                                         best.objective) {
  tryCatch(
    {
      state <- .np_progress_runtime$bandwidth_state
      if (isTRUE(state$nomad_native_progress)) {
        .np_progress_nomad_native_step_from_c(
          iteration = iteration,
          current.degree = current.degree,
          best.degree = best.degree,
          best.objective = best.objective,
          force = TRUE
        )
      } else {
        .np_progress_bandwidth_activity_step(done = iteration, force = TRUE)
      }
      list(0L, "")
    },
    interrupt = function(e) {
      list(2L, conditionMessage(e))
    },
    error = function(e) {
      list(1L, conditionMessage(e))
    }
  )
}

.np_progress_nomad_native_observer_report <- function(message) {
  .np_message("NOMAD progress observer disabled: ", message)
  invisible(NULL)
}

.np_nomad_native_progress_restart <- function(handle,
                                              restart_index,
                                              degree,
                                              best_record,
                                              restart_durations = numeric(),
                                              eval_offset = 0L) {
  state <- .np_nomad_native_progress_state(handle)
  if (is.null(state)) {
    return(invisible(NULL))
  }

  state$nomad_current_degree <- as.integer(degree)
  state$nomad_best_record <- .np_nomad_progress_best_record(
    incumbent = state$nomad_best_record,
    candidate = best_record,
    direction = "min"
  )
  state$nomad_restart_index <- as.integer(restart_index)
  state$nomad_restart_durations <- restart_durations
  eval_offset <- suppressWarnings(as.integer(eval_offset)[1L])
  if (is.na(eval_offset) || eval_offset < 0L) {
    eval_offset <- 0L
  }
  state$nomad_eval_offset <- eval_offset
  state$nomad_eval_id <- state$nomad_eval_offset
  state$last_done <- NULL
  state <- .np_degree_progress_step(
    state = state,
    done = NULL,
    detail = NULL,
    force = TRUE
  )
  handle$state <- state
  .np_progress_runtime$bandwidth_state <- state
  invisible(NULL)
}

.np_nomad_native_progress_end <- function(handle,
                                          degree,
                                          best_record,
                                          interrupted = FALSE) {
  if (is.null(handle) || isTRUE(handle$closed)) {
    return(invisible(NULL))
  }

  state <- .np_nomad_native_progress_state(handle)
  if (!is.null(state)) {
    state$nomad_current_degree <- as.integer(degree)
    state$nomad_best_record <- .np_nomad_progress_best_record(
      incumbent = state$nomad_best_record,
      candidate = best_record,
      direction = "min"
    )
    state <- .np_degree_progress_end(
      state = state,
      detail = NULL,
      interrupted = isTRUE(interrupted)
    )
  }
  .np_progress_runtime$bandwidth_state <- handle$old_state
  handle$closed <- TRUE
  invisible(NULL)
}

.np_nomad_native_progress_abort <- function(handle, detail = NULL) {
  if (is.null(handle) || isTRUE(handle$closed)) {
    return(invisible(NULL))
  }

  state <- .np_nomad_native_progress_state(handle)
  if (!is.null(state)) {
    .np_progress_abort(state = state, detail = detail)
  }
  .np_progress_runtime$bandwidth_state <- handle$old_state
  handle$closed <- TRUE
  invisible(NULL)
}

.np_nomad_search <- function(engine = c("nomad", "nomad+powell"),
                             baseline_record,
                             start_degree = NULL,
                             x0,
                             bbin,
                             lb,
                             ub,
                             eval_fun,
                             build_payload,
                             direction = c("min", "max"),
                             objective_name = "objective",
                             nmulti = 1L,
                             nomad.inner.nmulti = 0L,
                             random.seed = 42L,
                             display.nomad.progress = FALSE,
                             progress_state = NULL,
                             manage_progress_lifecycle = is.null(progress_state),
                             bind_bandwidth_runtime = FALSE,
                             handoff_before_build = FALSE,
                             remin = FALSE,
                             degree_spec = NULL,
                             start.lower = NULL,
                             start.upper = NULL,
                             coordinate.roles = NULL,
                             nomad.opts = list(),
                             native.r.bridge = FALSE,
                             source = "explicit",
                             reason = NULL,
                             progress_label = NULL) {
  engine <- match.arg(engine)
  direction <- match.arg(direction)
  source <- as.character(source)[1L]
  if (is.na(source) || !nzchar(source))
    source <- "explicit"
  if (is.null(progress_label))
    progress_label <- .np_degree_search_label(engine, source = source)
  .np_nomad_require_crs()

  state <- new.env(parent = emptyenv())
  state$trace_records <- list()
  state$trace_id <- 0L
  state$eval_id <- 0L
  state$visit_id <- 0L
  state$best_record <- baseline_record
  state$baseline_record <- baseline_record
  state$best_point <- NULL
  state$best_payload <- NULL
  state$interrupted <- FALSE
  state$error <- NULL
  state$progress_state <- NULL
  state$nomad.time <- NA_real_
  state$powell.time <- NA_real_
  state$restart_starts <- NULL
  state$restart_results <- NULL
  state$nomad.remin <- isTRUE(remin)
  state$nomad.remin.index <- NA_integer_
  state$nomad.remin.roundtrip <- NULL
  state$current_restart <- NA_integer_
  state$best_restart_index <- NA_integer_
  state$restart_durations <- numeric()
  state$restart_eval_id <- 0L
  state$external_output_seen <- character()
  state$native.r.bridge <- isTRUE(native.r.bridge)

  set_progress_state <- function(value) {
    state$progress_state <- value
    if (isTRUE(bind_bandwidth_runtime)) {
      .np_progress_runtime$bandwidth_state <- value
    }
    invisible(value)
  }

  nomad.nmulti <- npValidateNmulti(nmulti)
  nomad.inner.nmulti <- npValidateNonNegativeInteger(nomad.inner.nmulti, "nomad.inner.nmulti")
  nomad.opts <- if (is.null(nomad.opts)) list() else nomad.opts
  start_matrix <- .np_nomad_build_starts(
    x0 = x0,
    bbin = bbin,
    lb = lb,
    ub = ub,
    nmulti = nomad.nmulti,
    random.seed = random.seed,
    degree_spec = degree_spec,
    start.lower = start.lower,
    start.upper = start.upper
  )
  coordinate.roles <- .np_nomad_validate_coordinate_roles(
    coordinate.roles,
    ncol(start_matrix),
    where = ".np_nomad_search"
  )
  state$restart_starts <- lapply(
    seq_len(nrow(start_matrix)),
    function(i) as.numeric(start_matrix[i, ])
  )
  if (!is.null(degree_spec) && length(degree_spec$lower) > 0L) {
    q <- length(degree_spec$lower)
    degree.idx <- (ncol(start_matrix) - q + 1L):ncol(start_matrix)
    state$restart_degree_starts <- lapply(
      seq_len(nrow(start_matrix)),
      function(i) as.integer(start_matrix[i, degree.idx])
    )
    if (degree.idx[1L] > 1L) {
      state$restart_bandwidth_starts <- lapply(
        seq_len(nrow(start_matrix)),
        function(i) as.numeric(start_matrix[i, seq_len(degree.idx[1L] - 1L)])
      )
    } else {
      state$restart_bandwidth_starts <- replicate(nrow(start_matrix), numeric(0), simplify = FALSE)
    }
    state$restart_start_info <- list(
      basis = if (is.null(degree_spec$basis)) "glp" else degree_spec$basis,
      degree.start.policy = .np_lp_nomad_degree_start_policy(),
      lower = as.integer(degree_spec$lower),
      upper = as.integer(degree_spec$upper),
      user_supplied_start = isTRUE(degree_spec$user_supplied)
    )
  } else {
    state$restart_degree_starts <- NULL
    state$restart_bandwidth_starts <- NULL
    state$restart_start_info <- NULL
  }

  state$record_trace <- function(rec) {
    state$trace_id <- state$trace_id + 1L
    rec$trace_id <- state$trace_id
    state$trace_records[[length(state$trace_records) + 1L]] <- rec
    invisible(rec)
  }

  state$update_best <- function(rec, point) {
    if (!identical(rec$status, "ok") || !is.finite(rec$objective))
      return(invisible(NULL))

    incumbent <- if (is.null(state$best_record)) {
      if (identical(direction, "min")) Inf else -Inf
    } else {
      state$best_record$objective
    }

    if (.np_degree_better(rec$objective, incumbent, direction = direction)) {
      state$best_record <- rec
      state$best_point <- point
      state$best_restart_index <- state$current_restart
    }

    invisible(NULL)
  }

  wrapped_eval <- function(point) {
    point <- as.numeric(point)
    state$visit_id <- state$visit_id + 1L
    started <- proc.time()[3L]
    status <- "ok"
    msg <- NULL
    objective <- NA_real_
    degree <- integer(0)
    num.feval <- NA_real_
    num.feval.fast <- NA_real_

    result <- tryCatch(
      eval_fun(point),
      interrupt = function(e) {
        state$interrupted <- TRUE
        NULL
      },
      error = function(e) {
        status <<- "error"
        msg <<- conditionMessage(e)
        NULL
      }
    )

    if (isTRUE(state$interrupted))
      stop(structure(list(message = "interrupt"), class = c("interrupt", "condition")))

    if (!is.null(result)) {
      if (!is.list(result) || is.null(result$objective) || is.null(result$degree))
        stop("NOMAD fixed-point evaluator must return a list containing 'objective' and 'degree'")

      objective <- as.numeric(result$objective[1L])
      degree <- as.integer(result$degree)
      if (!is.null(result$num.feval))
        num.feval <- as.numeric(result$num.feval[1L])
      if (!is.null(result$num.feval.fast))
        num.feval.fast <- as.numeric(result$num.feval.fast[1L])
    }

    state$eval_id <- state$eval_id + 1L
    rec <- list(
      eval_id = state$eval_id,
      degree = degree,
      objective = objective,
      status = status,
      cached = FALSE,
      message = msg,
      elapsed = proc.time()[3L] - started,
      num.feval = num.feval,
      num.feval.fast = num.feval.fast
    )
    if (is.null(state$baseline_record))
      state$baseline_record <- rec
    state$record_trace(rec)
    state$update_best(rec, point = point)
    state$restart_eval_id <- state$restart_eval_id + 1L
    state$progress_state$nomad_current_degree <- degree
    state$progress_state$nomad_best_record <- state$best_record
    state$progress_state$nomad_eval_id <- state$eval_id
    state$progress_state$nomad_restart_index <- state$current_restart
    state$progress_state$nomad_restart_durations <- state$restart_durations
    state$progress_state <- .np_degree_progress_step(
      state = state$progress_state,
      done = state$restart_eval_id,
      detail = NULL
    )

    if (identical(status, "ok")) {
      if (identical(direction, "min")) {
        return(objective)
      } else {
        return(-objective)
      }
    }

    if (identical(direction, "min")) Inf else .Machine$double.xmax
  }

  baseline.degree <- if (is.null(start_degree)) {
    if (is.null(baseline_record)) integer(0) else baseline_record$degree
  } else {
    as.integer(start_degree)
  }
  if (is.null(progress_state)) {
    set_progress_state(.np_nomad_progress_begin(
      nmulti = nomad.nmulti,
      baseline_degree = baseline.degree,
      best_record = state$best_record,
      label = progress_label
    ))
  } else {
    set_progress_state(.np_nomad_progress_configure(
      state = progress_state,
      nmulti = nomad.nmulti,
      baseline_degree = baseline.degree,
      best_record = state$best_record,
      label = progress_label
    ))
  }

  run_nomad_solver <- function(start) {
    solver.opts <- .np_nomad_prepare_solver_opts(
      random.seed = random.seed,
      nomad.opts = nomad.opts,
      coordinate.roles = coordinate.roles,
      expected.length = length(start),
      geometry.policy = if (is.null(coordinate.roles)) "user-only" else "generate-central",
      where = ".np_nomad_search source geometry"
    )
    start <- as.numeric(start)

    if (isTRUE(state$native.r.bridge)) {
      native.eval <- function(point) {
        value <- as.numeric(wrapped_eval(point)[1L])
        if (!is.finite(value)) .Machine$double.xmax else value
      }
      native.option.vectors <- .np_nomad_native_option_vectors(solver.opts)
      native.call <- .np_nomad_native_r_callback_search(
        eval.f = native.eval,
        x0 = start,
        bbin = bbin,
        lb = lb,
        ub = ub,
        random.seed = random.seed,
        inner.start.count = nomad.inner.nmulti,
        option.names = native.option.vectors$names,
        option.values = native.option.vectors$values,
        display.nomad.progress = display.nomad.progress
      )
      native.value <- native.call$value
      native.status.raw <- if (!is.null(native.value$status)) native.value$status[1L] else 0L
      native.status.integer <- suppressWarnings(as.integer(native.status.raw))
      native.status.failed <- if (!is.na(native.status.integer)) {
        !identical(native.status.integer, 0L)
      } else {
        !(tolower(as.character(native.status.raw)) %in% c("", "ok", "success"))
      }
      native.result.status <- if (!is.null(native.value$result_status)) {
        suppressWarnings(as.integer(native.value$result_status[1L]))
      } else {
        0L
      }
      native.result.failed <- !is.na(native.result.status) && !identical(native.result.status, 0L)
      if (native.status.failed || native.result.failed) {
        native.message <- if (!is.null(native.value$message) &&
                              nzchar(as.character(native.value$message[1L]))) {
          as.character(native.value$message[1L])
        } else {
          "NOMAD rejected the supplied parameters"
        }
        stop(sprintf(
          "native NOMAD R-callback route failed (status=%s, result_status=%s): %s",
          as.character(native.status.raw),
          native.result.status,
          native.message
        ), call. = FALSE)
      }
      return(native.call)
    }

    stop("legacy R-level NOMAD fallback is retired for np NOMAD routes; this route must use the native crs C API or fail earlier as unsupported",
         call. = FALSE)
  }

  restart_results <- vector("list", nomad.nmulti)
  best_solution <- NULL
  nomad.elapsed <- 0
  for (i in seq_len(nomad.nmulti)) {
    state$current_restart <- as.integer(i)
    state$restart_eval_id <- 0L
    restart_degree <- if (!is.null(state$restart_degree_starts) &&
                          length(state$restart_degree_starts) >= i) {
      state$restart_degree_starts[[i]]
    } else if (is.null(start_degree)) {
      if (is.null(baseline_record)) integer(0) else baseline_record$degree
    } else {
      as.integer(start_degree)
    }

    if (i > 1L) {
      state$progress_state$nomad_current_degree <- restart_degree
      state$progress_state$nomad_best_record <- state$best_record
      state$progress_state$nomad_restart_index <- i
      state$progress_state$nomad_restart_durations <- state$restart_durations
      state$progress_state <- .np_degree_progress_step(
        state = state$progress_state,
        done = NULL,
        detail = NULL,
        force = TRUE
      )
    }

    pre_restart_best <- if (is.null(state$best_record) || is.null(state$best_record$objective)) {
      if (identical(direction, "min")) Inf else -Inf
    } else {
      state$best_record$objective
    }

    nomad.start <- proc.time()[3L]
    solution_i <- tryCatch(
      {
        nomad.call <- run_nomad_solver(as.numeric(start_matrix[i, ]))
        emitted <- .np_nomad_emit_external_output(
          lines = nomad.call$output,
          progress_state = state$progress_state,
          seen = state$external_output_seen
        )
        state$progress_state <- emitted$progress_state
        state$external_output_seen <- emitted$seen
        nomad.call$value
      },
      interrupt = function(e) {
        state$interrupted <- TRUE
        NULL
      },
      error = function(e) {
        state$error <- conditionMessage(e)
        NULL
      }
    )
    restart.elapsed <- proc.time()[3L] - nomad.start
    nomad.elapsed <- nomad.elapsed + restart.elapsed
    state$restart_durations <- c(state$restart_durations, restart.elapsed)
    state$progress_state$nomad_restart_durations <- state$restart_durations

    restart_results[[i]] <- list(
      restart = i,
      start = as.numeric(start_matrix[i, ]),
      degree.start = restart_degree,
      elapsed = restart.elapsed,
      status = if (is.null(solution_i)) {
        if (isTRUE(state$interrupted)) "interrupt" else "error"
      } else if (!is.null(solution_i$status)) {
        as.character(solution_i$status)
      } else {
        "ok"
      },
      message = if (is.null(solution_i)) state$error else solution_i$message,
      objective = if (!is.null(solution_i) && !is.null(solution_i$objective)) {
        raw.objective <- as.numeric(solution_i$objective[1L])
        if (identical(direction, "min")) raw.objective else -raw.objective
      } else {
        NA_real_
      },
      bbe = if (!is.null(solution_i) && !is.null(solution_i$bbe)) {
        as.numeric(solution_i$bbe[1L])
      } else {
        NA_real_
      },
      iterations = if (!is.null(solution_i) && !is.null(solution_i$iterations)) {
        as.numeric(solution_i$iterations[1L])
      } else {
        NA_real_
      },
      solution = if (!is.null(solution_i) && !is.null(solution_i$solution)) {
        as.numeric(solution_i$solution)
      } else {
        NULL
      }
    )

    if (isTRUE(state$interrupted))
      break

    if (is.null(solution_i) && !is.null(state$error)) {
      if (!is.null(state$progress_state)) {
        state$progress_state <- .np_progress_abort(
          state = state$progress_state,
          detail = state$error
        )
      }
      stop(state$error, call. = FALSE)
    }

    post_restart_best <- if (is.null(state$best_record) || is.null(state$best_record$objective)) {
      if (identical(direction, "min")) Inf else -Inf
    } else {
      state$best_record$objective
    }
    if (!is.null(solution_i) &&
        .np_degree_better(post_restart_best, pre_restart_best, direction = direction)) {
      best_solution <- solution_i
    }
  }

  if (isTRUE(state$nomad.remin) && !isTRUE(state$interrupted) &&
      !is.null(state$best_point) && length(state$best_point) == length(x0)) {
    remin.index <- length(restart_results) + 1L
    remin.start <- as.numeric(state$best_point)
    remin.degree <- if (!is.null(state$best_record$degree)) {
      as.integer(state$best_record$degree)
    } else {
      integer(0)
    }
    state$nomad.remin.roundtrip <- list(
      objective = as.numeric(state$best_record$objective[1L]),
      degree = remin.degree,
      note = "accepted best point reused; generic eval_fun is not side-effect-free"
    )
    state$current_restart <- as.integer(remin.index)
    state$restart_eval_id <- 0L
    state$nomad.remin.index <- as.integer(remin.index)
    state$restart_starts[[remin.index]] <- remin.start
    if (!is.null(state$restart_degree_starts))
      state$restart_degree_starts[[remin.index]] <- remin.degree
    if (!is.null(state$restart_bandwidth_starts)) {
      q <- length(remin.degree)
      state$restart_bandwidth_starts[[remin.index]] <- if (q > 0L) {
        as.numeric(remin.start[seq_len(length(remin.start) - q)])
      } else {
        as.numeric(remin.start)
      }
    }
    if (!is.null(state$progress_state)) {
      state$progress_state$nomad_current_degree <- remin.degree
      state$progress_state$nomad_best_record <- state$best_record
      state$progress_state$nomad_restart_index <- remin.index
      state$progress_state$nomad_restart_durations <- state$restart_durations
      state$progress_state <- .np_degree_progress_step(
        state = state$progress_state,
        done = NULL,
        detail = NULL,
        force = TRUE
      )
    }
    pre_restart_best <- state$best_record$objective
    nomad.start <- proc.time()[3L]
    solution_i <- tryCatch(
      {
        nomad.call <- run_nomad_solver(remin.start)
        emitted <- .np_nomad_emit_external_output(
          lines = nomad.call$output,
          progress_state = state$progress_state,
          seen = state$external_output_seen
        )
        state$progress_state <- emitted$progress_state
        state$external_output_seen <- emitted$seen
        nomad.call$value
      },
      interrupt = function(e) {
        state$interrupted <- TRUE
        NULL
      },
      error = function(e) {
        state$error <- conditionMessage(e)
        NULL
      }
    )
    restart.elapsed <- proc.time()[3L] - nomad.start
    nomad.elapsed <- nomad.elapsed + restart.elapsed
    state$restart_durations <- c(state$restart_durations, restart.elapsed)
    state$progress_state$nomad_restart_durations <- state$restart_durations
    restart_results[[remin.index]] <- list(
      restart = remin.index,
      remin = TRUE,
      start = remin.start,
      degree.start = remin.degree,
      elapsed = restart.elapsed,
      status = if (is.null(solution_i)) {
        if (isTRUE(state$interrupted)) "interrupt" else "error"
      } else if (!is.null(solution_i$status)) {
        as.character(solution_i$status)
      } else {
        "ok"
      },
      message = if (is.null(solution_i)) state$error else solution_i$message,
      objective = if (!is.null(solution_i) && !is.null(solution_i$objective)) {
        raw.objective <- as.numeric(solution_i$objective[1L])
        if (identical(direction, "min")) raw.objective else -raw.objective
      } else {
        NA_real_
      },
      bbe = if (!is.null(solution_i) && !is.null(solution_i$bbe)) {
        as.numeric(solution_i$bbe[1L])
      } else {
        NA_real_
      },
      iterations = if (!is.null(solution_i) && !is.null(solution_i$iterations)) {
        as.numeric(solution_i$iterations[1L])
      } else {
        NA_real_
      },
      solution = if (!is.null(solution_i) && !is.null(solution_i$solution)) {
        as.numeric(solution_i$solution)
      } else {
        NULL
      }
    )
    if (is.null(solution_i) && !is.null(state$error)) {
      if (!is.null(state$progress_state)) {
        state$progress_state <- .np_progress_abort(
          state = state$progress_state,
          detail = state$error
        )
      }
      stop(state$error, call. = FALSE)
    }
    post_restart_best <- if (is.null(state$best_record) || is.null(state$best_record$objective)) {
      if (identical(direction, "min")) Inf else -Inf
    } else {
      state$best_record$objective
    }
    if (!is.null(solution_i) &&
        .np_degree_better(post_restart_best, pre_restart_best, direction = direction)) {
      best_solution <- solution_i
    }
  }
  state$nomad.time <- nomad.elapsed
  state$restart_results <- restart_results

  if (!is.null(best_solution) &&
      is.null(state$best_point) &&
      length(best_solution$solution) == length(x0)) {
    state$best_point <- as.numeric(best_solution$solution)
  }

  if (is.null(state$best_point))
    stop(if (is.null(state$error)) {
      "automatic degree search failed to obtain any admissible fitted model"
    } else {
      state$error
    })

  if (identical(engine, "nomad+powell") && !is.null(state$progress_state)) {
    state$progress_state <- .np_nomad_progress_enter_powell(
      state = state$progress_state,
      degree = state$best_record$degree,
      best_record = state$best_record
    )
  }

  if (isTRUE(handoff_before_build) && !is.null(state$progress_state)) {
    state$progress_state$nomad_current_degree <- state$best_record$degree
    state$progress_state$nomad_best_record <- state$best_record
    state$progress_state <- .np_degree_progress_end(
      state = state$progress_state,
      detail = NULL,
      interrupted = state$interrupted
    )
    state$progress_state <- NULL
  }

  payload_result <- tryCatch(
    build_payload(
      point = state$best_point,
      best_record = state$best_record,
      solution = best_solution,
      interrupted = state$interrupted
    ),
    error = function(e) {
      if (!is.null(state$progress_state)) {
        state$progress_state <- .np_progress_abort(
          state = state$progress_state,
          detail = conditionMessage(e)
        )
      }
      stop(e)
    }
  )
  if (is.list(payload_result) && !is.null(payload_result$payload)) {
    state$best_payload <- payload_result$payload
    if (!is.null(payload_result$powell.time))
      state$powell.time <- as.numeric(payload_result$powell.time[1L])
    if (!is.null(payload_result$objective) &&
        .np_degree_better(payload_result$objective, state$best_record$objective, direction = direction)) {
      state$best_record$objective <- as.numeric(payload_result$objective[1L])
    }
  } else {
    state$best_payload <- payload_result
  }

  if (!is.null(state$progress_state)) {
    state$progress_state$nomad_current_degree <- state$best_record$degree
    state$progress_state$nomad_best_record <- state$best_record
    if (isTRUE(manage_progress_lifecycle)) {
      set_progress_state(.np_degree_progress_end(
        state = state$progress_state,
        detail = NULL,
        interrupted = state$interrupted
      ))
    }
  }

  list(
    method = engine,
    source = source,
    reason = reason,
    direction = direction,
    verify = FALSE,
    completed = !isTRUE(state$interrupted),
    certified = FALSE,
    interrupted = isTRUE(state$interrupted),
    baseline = state$baseline_record,
    best = state$best_record,
    best_payload = state$best_payload,
    best_point = state$best_point,
    n.unique = state$eval_id,
    n.visits = state$visit_id,
    n.cached = 0L,
    nomad.time = state$nomad.time,
    powell.time = state$powell.time,
    optim.time = sum(c(state$nomad.time, state$powell.time), na.rm = TRUE),
    grid.size = NA_integer_,
    best.restart = state$best_restart_index,
    nomad.remin = state$nomad.remin,
    nomad.remin.index = state$nomad.remin.index,
    nomad.remin.roundtrip = state$nomad.remin.roundtrip,
    restart.starts = state$restart_starts,
    restart.degree.starts = state$restart_degree_starts,
    restart.bandwidth.starts = state$restart_bandwidth_starts,
    restart.start.info = state$restart_start_info,
    restart.results = state$restart_results,
    trace = .np_degree_trace_to_frame(state$trace_records, objective_name = objective_name),
    native = isTRUE(state$native.r.bridge),
    native.diagnostics = if (isTRUE(state$native.r.bridge)) {
      c(list(
        route_native = TRUE,
        callback_mode = "R",
        native_symbol = "C_np_nomad_r_callback_native_search",
        crs.callback.evaluations = as.integer(.np_native_nomad_field(best_solution, "callback_evaluations", NA_integer_)[1L]),
        blackbox.evaluations = as.integer(.np_native_nomad_field(best_solution, "bbe", NA_integer_)[1L])
      ), .np_native_nomad_cache_diagnostics(best_solution))
    } else {
      NULL
    }
  )
}
