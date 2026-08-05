test_that("npindex LP via-npreg evaluators preserve outer bandwidth progress state", {
  runtime <- getFromNamespace(".np_progress_runtime", "np")
  old.state <- runtime$bandwidth_state
  old.forward.active <- runtime$bandwidth_forward_active
  old.forward.state <- runtime$bandwidth_forward_state
  on.exit({
    runtime$bandwidth_state <- old.state
    runtime$bandwidth_forward_active <- old.forward.active
    runtime$bandwidth_forward_state <- old.forward.state
  }, add = TRUE)

  eval.ichimura <- getFromNamespace(".npindexbw_eval_ichimura_lp_via_npreg", "np")
  eval.kleinspady <- getFromNamespace(".npindexbw_eval_kleinspady_lp_via_npreg", "np")

  build.leaf <- function(index, ydat, h, bws, spec) {
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
    current = 2L,
    total = 5L
  )

  with_np_progress_bindings(
    list(
      .npindexbw_build_lp_regression_leaf = build.leaf,
      .npregbw_eval_only = function(...) {
        getFromNamespace(".np_progress_bandwidth_activity_step", "np")(done = 3L)
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
      expect_identical(runtime$bandwidth_state$last_done, 3L)
      expect_false(runtime$bandwidth_forward_active)
      expect_equal(out$objective, 1.25)
      expect_equal(out$num.feval.fast, 3)
    }
  )

  with_np_progress_bindings(
    list(
      .npindexbw_build_lp_regression_leaf = build.leaf,
      .npregbw_eval_only = function(...) {
        getFromNamespace(".np_progress_bandwidth_activity_step", "np")(done = 4L)
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
      expect_identical(runtime$bandwidth_state$last_done, 4L)
      expect_false(runtime$bandwidth_forward_active)
      expect_equal(out$objective, 77)
      expect_equal(out$num.feval.fast, 0)
    }
  )
})
