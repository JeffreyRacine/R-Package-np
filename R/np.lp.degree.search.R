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
  "Selecting degree and bandwidth"
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

.np_nomad_require_crs <- function(version_fn = utils::packageVersion,
                                  load_namespace = loadNamespace,
                                  minimum_version = "0.15-41") {
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

.np_lp_nomad_dim_budget <- function(nobs) {
  nobs <- as.integer(nobs[1L])
  if (!is.finite(nobs) || is.na(nobs) || nobs <= 1L)
    return(1L)
  max(1L, as.integer(floor(0.25 * (nobs - 1L))))
}

.np_lp_nomad_proposal_upper <- function(lower, upper) {
  q <- length(lower)
  if (!q)
    return(integer(0))

  lower <- as.integer(lower)
  upper <- as.integer(upper)
  as.integer(pmin(
    upper,
    pmax(lower + 1L, ceiling(upper / max(1L, q)))
  ))
}

.np_lp_nomad_reduce_degree_start <- function(degree,
                                             lower,
                                             basis = c("glp", "additive", "tensor"),
                                             dim_budget = Inf) {
  basis <- match.arg(basis)
  out <- as.integer(degree)
  lower <- as.integer(lower)

  if (!is.finite(dim_budget))
    return(out)

  repeat {
    current.dim <- tryCatch(
      dim_basis(basis = basis, degree = out),
      error = function(e) Inf
    )
    if (is.finite(current.dim) && current.dim <= dim_budget)
      break

    idx <- which(out > lower)
    if (!length(idx))
      break

    drop.idx <- idx[which.max(out[idx])][1L]
    out[drop.idx] <- out[drop.idx] - 1L
  }

  out
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

  if (!q)
    return(matrix(integer(0), nrow = nstart, ncol = 0L))

  starts <- matrix(0L, nrow = nstart, ncol = q)
  dim_budget <- .np_lp_nomad_dim_budget(nobs)
  proposal.upper <- .np_lp_nomad_proposal_upper(lower, upper)

  initial <- as.integer(pmax(lower, pmin(upper, initial)))
  if (!isTRUE(user_supplied))
    initial <- .np_lp_nomad_reduce_degree_start(
      degree = initial,
      lower = lower,
      basis = basis,
      dim_budget = dim_budget
    )
  starts[1L, ] <- initial

  if (nstart <= 1L)
    return(starts)

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

  max.tries <- max(1L, as.integer(max.tries[1L]))
  for (j in 2:nstart) {
    accepted <- FALSE
    fallback <- starts[1L, ]
    for (k in seq_len(max.tries)) {
      candidate <- vapply(
        seq_len(q),
        function(i) sample.int(proposal.upper[i] - lower[i] + 1L, 1L) + lower[i] - 1L,
        integer(1L)
      )
      fallback <- candidate
      candidate.dim <- tryCatch(
        dim_basis(basis = basis, degree = candidate),
        error = function(e) Inf
      )
      if (is.finite(candidate.dim) && candidate.dim <= dim_budget) {
        starts[j, ] <- as.integer(candidate)
        accepted <- TRUE
        break
      }
    }

    if (!accepted) {
      starts[j, ] <- .np_lp_nomad_reduce_degree_start(
        degree = fallback,
        lower = lower,
        basis = basis,
        dim_budget = dim_budget
      )
    }
  }

  starts
}

.np_nomad_build_starts <- function(x0,
                                   bbin,
                                   lb,
                                   ub,
                                   nmulti = 1L,
                                   random.seed = 42L,
                                   degree_spec = NULL) {
  n <- length(lb)
  if (length(ub) != n || length(bbin) != n)
    stop("NOMAD start construction requires matching lengths for x0/bbin/lb/ub")

  nstart <- npValidateNmulti(nmulti)
  starts <- matrix(0, nrow = nstart, ncol = n)

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

  for (j in seq_len(nstart)) {
    for (i in seq_len(n)) {
      lo <- if (is.finite(lb[i])) lb[i] else -1
      hi <- if (is.finite(ub[i])) ub[i] else 1
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

.np_nomad_progress_detail <- function(current_degree,
                                      best_record,
                                      iteration = NULL,
                                      cumulative_iteration = NULL,
                                      restart_index = NULL,
                                      nmulti = 1L,
                                      restart_durations = numeric(),
                                      elapsed = NULL) {
  fields <- character()
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
                                     best_record) {
  state <- .np_progress_begin(
    label = .np_degree_progress_label(),
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

  state
}

.np_nomad_progress_configure <- function(state,
                                         nmulti,
                                         baseline_degree,
                                         best_record) {
  if (is.null(state))
    return(state)

  state$label <- .np_degree_progress_label()
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

  state
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

.np_nomad_powell_context_label <- function(degree) {
  sprintf(
    "Refining NOMAD solution with one Powell hot start at degree %s",
    .np_degree_format_degree(degree)
  )
}

.np_nomad_powell_progress_detail <- function(current_degree,
                                             best_record,
                                             elapsed = NULL) {
  fields <- character()

  if (is.finite(elapsed) && !is.na(elapsed) && elapsed >= 0)
    fields <- c(fields, sprintf("elapsed %ss", .np_progress_fmt_num(elapsed)))

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
    elapsed = max(0, now - state$started)
  )
}

.np_nomad_progress_enter_powell <- function(state,
                                            degree,
                                            best_record) {
  if (is.null(state))
    return(state)

  state$label <- .np_nomad_powell_progress_label()
  state$unknown_total_fields <- .np_nomad_powell_progress_fields
  state$nomad_current_degree <- as.integer(degree)
  state$nomad_best_record <- best_record
  .np_progress_step_at(
    state = state,
    now = .np_progress_now(),
    done = state$last_done,
    force = TRUE
  )
}

.np_nomad_with_powell_progress <- function(degree, expr) {
  .np_progress_bandwidth_set_context(.np_nomad_powell_context_label(degree))
  on.exit(.np_progress_bandwidth_set_context(NULL), add = TRUE)
  force(expr)
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
                             degree_spec = NULL,
                             nomad.opts = list()) {
  engine <- match.arg(engine)
  direction <- match.arg(direction)
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
  state$current_restart <- NA_integer_
  state$best_restart_index <- NA_integer_
  state$restart_durations <- numeric()
  state$restart_eval_id <- 0L

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
    degree_spec = degree_spec
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
      dim_budget = .np_lp_nomad_dim_budget(degree_spec$nobs),
      proposal.upper = .np_lp_nomad_proposal_upper(
        lower = degree_spec$lower,
        upper = degree_spec$upper
      ),
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
      extra.fields <- result[setdiff(names(result), c("objective", "degree", "num.feval"))]
    } else {
      extra.fields <- list()
    }

    state$eval_id <- state$eval_id + 1L
    rec <- c(list(
      eval_id = state$eval_id,
      degree = degree,
      objective = objective,
      status = status,
      cached = FALSE,
      message = msg,
      elapsed = proc.time()[3L] - started,
      num.feval = num.feval
    ), extra.fields)
    if (is.null(state$baseline_record))
      state$baseline_record <- rec
    state$record_trace(rec)
    state$update_best(rec, point = point)
    state$restart_eval_id <- state$restart_eval_id + 1L
    if (!is.null(state$progress_state)) {
      state$progress_state$nomad_current_degree <- degree
      state$progress_state$nomad_best_record <- state$best_record
      state$progress_state$nomad_eval_id <- state$eval_id
      state$progress_state$nomad_restart_index <- state$current_restart
      state$progress_state$nomad_restart_durations <- state$restart_durations
      set_progress_state(.np_degree_progress_step(
        state = state$progress_state,
        done = state$restart_eval_id,
        detail = NULL
      ))
    }

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
      best_record = state$best_record
    ))
  } else {
    set_progress_state(.np_nomad_progress_configure(
      state = progress_state,
      nmulti = nomad.nmulti,
      baseline_degree = baseline.degree,
      best_record = state$best_record
    ))
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
      if (!is.null(state$progress_state)) {
        state$progress_state$nomad_current_degree <- restart_degree
        state$progress_state$nomad_best_record <- state$best_record
        state$progress_state$nomad_restart_index <- i
        state$progress_state$nomad_restart_durations <- state$restart_durations
        set_progress_state(.np_degree_progress_step(
          state = state$progress_state,
          done = NULL,
          detail = NULL,
          force = TRUE
        ))
      }
    }

    pre_restart_best <- if (is.null(state$best_record) || is.null(state$best_record$objective)) {
      if (identical(direction, "min")) Inf else -Inf
    } else {
      state$best_record$objective
    }

    nomad.start <- proc.time()[3L]
    solution_i <- tryCatch(
      {
        solver.opts <- utils::modifyList(
          list(
            SEED = as.integer(random.seed),
            RNG_ALT_SEEDING = TRUE
          ),
          nomad.opts
        )
        crs::snomadr(
          eval.f = wrapped_eval,
          n = length(x0),
          bbin = as.integer(bbin),
          bbout = 0L,
          x0 = as.numeric(start_matrix[i, ]),
          lb = as.double(lb),
          ub = as.double(ub),
          nmulti = nomad.inner.nmulti,
          random.seed = as.integer(random.seed),
          opts = solver.opts,
          display.nomad.progress = display.nomad.progress,
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
    restart.elapsed <- proc.time()[3L] - nomad.start
    nomad.elapsed <- nomad.elapsed + restart.elapsed
    state$restart_durations <- c(state$restart_durations, restart.elapsed)
    if (!is.null(state$progress_state))
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

  if (!is.null(best_solution) && is.null(state$best_point) && length(best_solution$solution) == length(x0)) {
    state$best_point <- as.numeric(best_solution$solution)
  }

  if (is.null(state$best_point))
    stop(if (is.null(state$error)) {
      "automatic degree search failed to obtain any admissible fitted model"
    } else {
      state$error
    })

  if (identical(engine, "nomad+powell") && !is.null(state$progress_state)) {
    set_progress_state(.np_nomad_progress_enter_powell(
      state = state$progress_state,
      degree = state$best_record$degree,
      best_record = state$best_record
    ))
  }

  if (isTRUE(handoff_before_build) && !is.null(state$progress_state)) {
    state$progress_state$nomad_current_degree <- state$best_record$degree
    state$progress_state$nomad_best_record <- state$best_record
    set_progress_state(.np_degree_progress_end(
      state = state$progress_state,
      detail = NULL,
      interrupted = state$interrupted
    ))
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
        set_progress_state(.np_progress_abort(
          state = state$progress_state,
          detail = conditionMessage(e)
        ))
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
    restart.starts = state$restart_starts,
    restart.degree.starts = state$restart_degree_starts,
    restart.bandwidth.starts = state$restart_bandwidth_starts,
    restart.start.info = state$restart_start_info,
    restart.results = state$restart_results,
    trace = .np_degree_trace_to_frame(state$trace_records, objective_name = objective_name)
  )
}
