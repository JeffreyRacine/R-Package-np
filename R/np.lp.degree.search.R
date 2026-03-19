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

  candidates <- lapply(seq_len(ncon), function(j) seq.int(lower[j], upper[j]))

  list(lower = lower, upper = upper, candidates = candidates)
}

.np_degree_search_engine_controls <- function(search.engine) {
  match.arg(search.engine, c("nomad+powell", "cell", "nomad"))
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

.np_degree_progress_label <- function() {
  "Selecting polynomial degree and bandwidth"
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
  fields <- c(phase)

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

.np_degree_progress_begin <- function(total = NULL, detail = NULL) {
  state <- .np_progress_begin(
    label = .np_degree_progress_label(),
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
                              objective_name = "objective") {
  method <- match.arg(method)
  direction <- match.arg(direction)
  trace_level <- match.arg(trace_level)
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
        num.feval = NA_real_
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
      num.feval = num.feval
    )
    assign(key, rec, envir = state$cache)
    state$record_trace(rec)
    if (identical(status, "ok"))
      state$update_best(rec, payload = payload)
    rec
  }

  baseline_degree <- as.integer(baseline_degree)
  start_degree <- as.integer(start_degree)
  restart_starts <- .np_degree_restart_starts(
    candidates = candidates,
    restarts = restarts,
    exclude = list(start_degree, baseline_degree)
  )

  tryCatch({
    .np_progress_note(sprintf(
      "Automatic polynomial degree search baseline %s",
      .np_degree_format_degree(baseline_degree)
    ))
    state$baseline_record <- state$evaluate(baseline_degree)

    if (identical(method, "exhaustive")) {
      .np_progress_note(sprintf(
        "Exhaustive automatic polynomial degree search over %s degree combinations on %s",
        format(grid.size, scientific = FALSE, trim = TRUE),
        .np_degree_format_search_space(candidates)
      ))
      state$progress_state <- .np_degree_progress_begin(
        total = grid.size,
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
      .np_progress_note(sprintf(
        "Coordinate automatic polynomial degree search over %s (max %s search evaluations)",
        .np_degree_format_search_space(candidates),
        format(
          .np_degree_coordinate_visit_cap(
            candidates = candidates,
            max_cycles = max_cycles,
            restart_total = restart_total
          ),
          scientific = FALSE,
          trim = TRUE
        )
      ))
      state$progress_state <- .np_degree_progress_begin(
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
        .np_progress_note(sprintf(
          "Exhaustively certifying automatic polynomial degree search over %s degree combinations (re-optimizing bandwidths)",
          format(grid.size, scientific = FALSE, trim = TRUE)
        ))
        state$progress_state <- .np_degree_progress_begin(
          total = grid.size,
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

  list(
    method = method,
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
    grid.size = grid.size,
    restart.starts = restart_starts,
    trace = .np_degree_trace_to_frame(state$trace_records, objective_name = objective_name)
  )
}

.np_nomad_require_crs <- function() {
  if (!requireNamespace("crs", quietly = TRUE)) {
    stop(
      "automatic degree search with search.engine='nomad' requires the 'crs' package; install.packages('crs')",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.np_nomad_progress_detail <- function(engine,
                                      eval_id,
                                      current_degree,
                                      best_record,
                                      objective_name = "objective") {
  fields <- c(engine)

  if (!is.null(eval_id) && is.finite(eval_id) && eval_id >= 0L)
    fields <- c(fields, sprintf("eval %s", format(eval_id)))

  if (!is.null(current_degree))
    fields <- c(fields, sprintf("deg %s", .np_degree_format_degree(current_degree)))

  fields <- c(fields, .np_degree_progress_best_detail(
    best_record = best_record,
    objective_name = objective_name
  ))

  paste(fields, collapse = ", ")
}

.np_nomad_search <- function(engine = c("nomad", "nomad+powell"),
                             baseline_record,
                             x0,
                             bbin,
                             lb,
                             ub,
                             eval_fun,
                             build_payload,
                             direction = c("min", "max"),
                             objective_name = "objective",
                             random.seed = 42L,
                             display.nomad.progress = FALSE) {
  engine <- match.arg(engine)
  direction <- match.arg(direction)
  .np_nomad_require_crs()

  state <- new.env(parent = emptyenv())
  state$trace_records <- list()
  state$trace_id <- 0L
  state$eval_id <- 0L
  state$visit_id <- 0L
  state$best_record <- baseline_record
  state$best_point <- NULL
  state$best_payload <- NULL
  state$interrupted <- FALSE
  state$error <- NULL
  state$progress_state <- NULL

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
      num.feval = num.feval
    )
    state$record_trace(rec)
    state$update_best(rec, point = point)
    state$progress_state <- .np_degree_progress_step(
      state = state$progress_state,
      done = state$eval_id,
      detail = .np_nomad_progress_detail(
        engine = engine,
        eval_id = state$eval_id,
        current_degree = degree,
        best_record = state$best_record,
        objective_name = objective_name
      )
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

  .np_progress_note(sprintf(
    "NOMAD automatic polynomial degree search starting from %s",
    .np_degree_format_degree(baseline_record$degree)
  ))
  state$progress_state <- .np_degree_progress_begin(
    detail = .np_nomad_progress_detail(
      engine = engine,
      eval_id = 0L,
      current_degree = baseline_record$degree,
      best_record = state$best_record,
      objective_name = objective_name
    )
  )

  solution <- tryCatch(
    {
      crs::snomadr(
        eval.f = wrapped_eval,
        n = length(x0),
        bbin = as.integer(bbin),
        bbout = 0L,
        x0 = as.numeric(x0),
        lb = as.double(lb),
        ub = as.double(ub),
        nmulti = 0L,
        random.seed = as.integer(random.seed),
        opts = list(),
        display.nomad.progress = display.nomad.progress
        ,
        snomadr.environment = environment(wrapped_eval)
      )
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

  if (!is.null(solution) && is.null(state$best_point) && length(solution$solution) == length(x0)) {
    state$best_point <- as.numeric(solution$solution)
  }

  if (is.null(state$best_point))
    stop(if (is.null(state$error)) {
      "automatic degree search failed to obtain any admissible fitted model"
    } else {
      state$error
    })

  if (identical(engine, "nomad+powell")) {
    .np_progress_note("Refining NOMAD solution with one Powell hot start")
  }

  payload_result <- build_payload(
    point = state$best_point,
    best_record = state$best_record,
    solution = solution,
    interrupted = state$interrupted
  )
  if (is.list(payload_result) && !is.null(payload_result$payload)) {
    state$best_payload <- payload_result$payload
    if (!is.null(payload_result$objective) &&
        .np_degree_better(payload_result$objective, state$best_record$objective, direction = direction)) {
      state$best_record$objective <- as.numeric(payload_result$objective[1L])
    }
  } else {
    state$best_payload <- payload_result
  }

  state$progress_state <- .np_degree_progress_end(
    state = state$progress_state,
    detail = .np_degree_progress_best_detail(
      best_record = state$best_record,
      objective_name = objective_name
    ),
    interrupted = state$interrupted
  )

  list(
    method = engine,
    verify = FALSE,
    completed = !isTRUE(state$interrupted),
    certified = FALSE,
    interrupted = isTRUE(state$interrupted),
    baseline = baseline_record,
    best = state$best_record,
    best_payload = state$best_payload,
    best_point = state$best_point,
    n.unique = state$eval_id,
    n.visits = state$visit_id,
    n.cached = 0L,
    grid.size = NA_integer_,
    restart.starts = list(),
    trace = .np_degree_trace_to_frame(state$trace_records, objective_name = objective_name)
  )
}
