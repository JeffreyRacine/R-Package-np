test_that("npindex LP via-npreg evaluators forward heartbeat but not child counters", {
  runtime <- getFromNamespace(".np_progress_runtime", "npRmpi")
  old.state <- runtime$bandwidth_state
  old.forward.active <- runtime$bandwidth_forward_active
  old.forward.state <- runtime$bandwidth_forward_state
  on.exit({
    runtime$bandwidth_state <- old.state
    runtime$bandwidth_forward_active <- old.forward.active
    runtime$bandwidth_forward_state <- old.forward.state
  }, add = TRUE)

  eval.ichimura <- getFromNamespace(".npindexbw_eval_ichimura_lp_via_npreg", "npRmpi")
  eval.kleinspady <- getFromNamespace(".npindexbw_eval_kleinspady_lp_via_npreg", "npRmpi")
  build.leaf <- function(descriptor, index, ydat, h, bws, spec) {
    list(
      xdat = data.frame(index = index),
      bws = list()
    )
  }

  outer.state <- list(
    id = "outer-npindex",
    enabled = FALSE,
    known_total = FALSE,
    nomad_native_progress = FALSE,
    label = "Bandwidth selection",
    progress_provider = "npindex",
    degree = c(1L, 2L),
    bandwidth_multistart_current = 1L,
    current = 2L,
    last_done = 2L,
    total = 5L
  )

  seen.done <- list()
  heartbeat.step <- function(state, now, done = NULL, detail = NULL, force = FALSE) {
    seen.done <<- c(seen.done, list(done))
    state$heartbeat_calls <- if (is.null(state$heartbeat_calls)) 1L else state$heartbeat_calls + 1L
    state
  }

  with_npRmpi_degree_bindings(
    list(
      .npindexbw_lp_regression_leaf = build.leaf,
      .np_progress_step_at = heartbeat.step,
      .npregbw_eval_only = function(...) {
        getFromNamespace(".np_progress_bandwidth_activity_step", "npRmpi")(done = 3L)
        getFromNamespace(".np_progress_bandwidth_activity_step", "npRmpi")(done = 99L)
        runtime$bandwidth_state <- list(id = "inner-npreg", current = 1L)
        list(objective = 1.25, num.feval.fast = 3L)
      }
    ),
    {
      runtime$bandwidth_state <- outer.state
      out <- eval.ichimura(
        index = c(0.1, 0.2, 0.3),
        ydat = c(1, 2, 3),
        h = 1,
        bws = list(),
        spec = list(),
        invalid.penalty = 99
      )

      expect_identical(runtime$bandwidth_state$id, "outer-npindex")
      expect_identical(runtime$bandwidth_state$last_done, 2L)
      expect_identical(runtime$bandwidth_state[names(outer.state)], outer.state)
      expect_identical(runtime$bandwidth_state$heartbeat_calls, 2L)
      expect_identical(seen.done, list(NULL, NULL))
      expect_false(runtime$bandwidth_forward_active)
      expect_equal(out$objective, 1.25)
      expect_equal(out$num.feval.fast, 3)
    }
  )

  with_npRmpi_degree_bindings(
    list(
      .npindexbw_lp_regression_leaf = build.leaf,
      .np_progress_step_at = heartbeat.step,
      .npregbw_eval_only = function(...) {
        getFromNamespace(".np_progress_bandwidth_activity_step", "npRmpi")(done = 4L)
        runtime$bandwidth_state <- list(id = "inner-npreg", current = 1L)
        stop("inner bandwidth failure", call. = FALSE)
      }
    ),
    {
      runtime$bandwidth_state <- outer.state
      out <- eval.kleinspady(
        index = c(0.1, 0.2, 0.3),
        ydat = c(0, 1, 1),
        h = 1,
        bws = list(),
        spec = list(),
        invalid.penalty = 77
      )

      expect_identical(runtime$bandwidth_state$id, "outer-npindex")
      expect_identical(runtime$bandwidth_state$last_done, 2L)
      expect_identical(runtime$bandwidth_state[names(outer.state)], outer.state)
      expect_false(runtime$bandwidth_forward_active)
      expect_equal(out$objective, 77)
      expect_equal(out$num.feval.fast, 0)
    }
  )
})
