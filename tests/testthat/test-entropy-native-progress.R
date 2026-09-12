# Native-to-R entropy progress must be observational only.
with_entropy_progress_fixture <- function(code) {
  ns <- asNamespace("np")
  names <- c(".np_progress_is_interactive", ".np_progress_render_single_line",
             ".np_progress_signal_from_c")
  old <- lapply(names, get, envir = ns, inherits = FALSE)
  set <- function(name, value) {
    unlockBinding(name, ns); assign(name, value, ns); lockBinding(name, ns)
  }
  on.exit(for (i in seq_along(names)) set(names[i], old[[i]]), add = TRUE)
  opt <- options(np.messages = TRUE, np.progress.interval.sec = 0,
                 np.progress.start.grace.known.sec = 0,
                 np.progress.start.grace.unknown.sec = 0)
  on.exit(options(opt), add = TRUE)
  events <- new.env(parent = emptyenv())
  events$native <- list(); events$render <- list()
  set(names[1L], function() TRUE)
  set(names[2L], function(snapshot, event) {
    events$render[[length(events$render) + 1L]] <- list(snapshot, event)
  })
  signal <- old[[3L]]
  set(names[3L], function(...) {
    events$native[[length(events$native) + 1L]] <- list(...)
    gc(FALSE)
    signal(...)
  })
  code(events, ns)
}

test_that("entropy native heartbeat preserves work and does not advance replications", {
  with_entropy_progress_fixture(function(events, ns) {
    run <- if ("np" == "npRmpi") get(".npRmpi_with_local_regression", ns) else
      function(expr) force(expr)
    run({
      x <- seq(-2, 2, length.out = 300L); y <- sin(x)
      counts <- matrix(1, nrow = length(x), ncol = 2L)
      cases <- list(
        list("C_np_entropy_bivariate_summation", x, y, c(.7,.8,.9,1)),
        list("C_np_entropy_bivariate_summation_xindex", x, y,
             matrix(rep(seq_along(x), each = 2L), 2L), c(.7,.8,.9,1)),
        list("C_np_entropy_symmetric_summation_counts", x, counts, .7),
        list("C_np_entropy_univariate_summation_counts", x, counts, counts, c(.7,.8)),
        list("C_np_entropy_gaussian_integrand", rbind(x,y), x, y, c(.7,.8,.9,1)))
      context <- .np_bootstrap_progress_begin(9L, "Bootstrap replications")
      on.exit(.np_progress_activity_end(context), add = TRUE)
      owner <- .np_progress_registry$active_id
      for (case in cases) {
        before <- length(events$native)
        expected <- do.call(.Call, c(case, list(PACKAGE = "np")))
        expect_identical(length(events$native), before)
        actual <- .np_entropy_compute(expr =
          do.call(.Call, c(case, list(PACKAGE = "np"))))
        expect_identical(actual, expected)
        expect_gt(length(events$native), before)
        expect_identical(context$done, 0L)
        expect_identical(.np_progress_registry$active_id, owner)
      }
      .np_progress_activity_end(context, completed = TRUE)
      expect_null(.np_progress_registry$active_id)
      expect_null(.np_progress_runtime$fit_forward)
      expect_null(.np_progress_runtime$fit_state)
    })
  })
})

test_that("native entropy errors clear activation and the outer owner", {
  with_entropy_progress_fixture(function(events, ns) {
    run <- if ("np" == "npRmpi") get(".npRmpi_with_local_regression", ns) else
      function(expr) force(expr)
    run({
      expect_error(.np_progress_activity_run("Computing entropy statistic",
        .np_entropy_compute(expr = .Call("C_np_entropy_bivariate_summation",
          1:3, 1:3, c(.7,.8,.9,1), PACKAGE = "np"))),
        "numeric")
      expect_null(.np_progress_registry$active_id)
      expect_null(.np_progress_runtime$fit_forward)
      expect_null(.np_progress_runtime$fit_state)
      x <- seq(-2, 2, length.out = 300L)
      before <- length(events$native)
      invisible(.Call("C_np_entropy_bivariate_summation",
        x, sin(x), c(.7,.8,.9,1), PACKAGE = "np"))
      expect_identical(length(events$native), before)
      expect_true(is.finite(.np_progress_activity_run("Computing entropy statistic",
        .np_entropy_bivariate_gaussian_summation(x, sin(x), .7, .8, c(.9,1)))))
      expect_gt(length(events$native), before)
      expect_null(.np_progress_registry$active_id)
    })
  })
})
