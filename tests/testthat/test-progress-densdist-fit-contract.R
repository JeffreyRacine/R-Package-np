densdist_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

densdist_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

capture_densdist_fit_progress_trace <- function(expr,
                                                force_renderer = "single_line",
                                                now = function() 0,
                                                interactive = TRUE,
                                                master = TRUE) {
  expr_env <- parent.frame()
  expr <- substitute(expr)
  trace <- list()

  recorder <- function(snapshot, event = c("render", "finish", "abort")) {
    event <- match.arg(event)
    if (isTRUE(getFromNamespace(".np_progress_is_message_muffled", "npRmpi")())) {
      return(invisible(snapshot))
    }
    trace[[length(trace) + 1L]] <<- list(
      event = event,
      renderer = force_renderer,
      id = snapshot$id,
      kind = snapshot$kind,
      current = snapshot$current,
      total = snapshot$total,
      detail = snapshot$detail,
      line = snapshot$line,
      started_at = snapshot$started_at,
      now = snapshot$now,
      last_width = snapshot$last_width
    )
    invisible(snapshot)
  }

  value <- with_nprmpi_progress_bindings(
    c(
      if (!is.null(force_renderer)) {
        list(.np_progress_renderer_for_surface = function(surface, capability) force_renderer)
      } else {
        list()
      },
      list(
        .np_progress_render_legacy = recorder,
        .np_progress_render_single_line = recorder,
        .np_progress_is_interactive = function() interactive,
        .np_progress_is_master = function() master,
        .np_progress_now = now,
        .np_progress_output_width = function() 500L
      )
    ),
    {
      reset <- getFromNamespace(".np_progress_reset_registry", "npRmpi")
      reset()
      on.exit(reset(), add = TRUE)
      eval(expr, envir = expr_env)
    }
  )

  list(
    value = value,
    trace = trace,
    final_line = if (length(trace)) trace[[length(trace)]]$line else NULL
  )
}

make_densdist_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 24L
  dat <- data.frame(x = seq(-0.9, 0.9, length.out = n))

  list(
    dat = dat,
    tdat = dat["x"],
    n = n,
    dens.bw = npudensbw(
      dat = dat["x"],
      bws = 0.24,
      bandwidth.compute = FALSE
    ),
    dist.bw = npudistbw(
      dat = dat["x"],
      bws = 0.24,
      bandwidth.compute = FALSE
    )
  )
}

test_that("npudens direct bws fit emits single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_densdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_densdist_fit_progress_trace(
    npudens(
      bws = fixture$dens.bw,
      tdat = fixture$tdat
    ),
    now = densdist_fit_progress_time_counter()
  )

  lines <- densdist_fit_progress_lines(actual)

  expect_s3_class(actual$value, "npdensity")
  expect_true(any(grepl(
    sprintf("^\\[npRmpi\\] Fitting density 1/%d \\([0-9]+\\.[0-9]%%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", fixture$n),
    lines
  )))
  expect_true(any(grepl(
    sprintf("^\\[npRmpi\\] Fitting density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )))
})

test_that("npudens bw to fit route hands off immediately into single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_densdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_densdist_fit_progress_trace(
    npudens(
      tdat = fixture$tdat,
      nmulti = 2L
    ),
    now = densdist_fit_progress_time_counter()
  )

  lines <- densdist_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting density 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", fixture$n),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )

  expect_s3_class(actual$value, "npdensity")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npudist direct bws fit emits single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_densdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.unknown.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_densdist_fit_progress_trace(
    npudist(
      bws = fixture$dist.bw,
      tdat = fixture$tdat
    ),
    now = densdist_fit_progress_time_counter()
  )

  lines <- densdist_fit_progress_lines(actual)

  expect_s3_class(actual$value, "npdistribution")
  expect_true(any(grepl(
    sprintf("^\\[npRmpi\\] Fitting distribution 1/%d \\([0-9]+\\.[0-9]%%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", fixture$n),
    lines
  )))
  expect_true(any(grepl(
    sprintf("^\\[npRmpi\\] Fitting distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )))
})

test_that("npudist bw to fit route hands off immediately into single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_densdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_densdist_fit_progress_trace(
    npudist(
      tdat = fixture$tdat,
      nmulti = 2L
    ),
    now = densdist_fit_progress_time_counter()
  )

  lines <- densdist_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting distribution 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", fixture$n),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )

  expect_s3_class(actual$value, "npdistribution")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})
