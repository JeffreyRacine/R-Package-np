condensdist_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

condensdist_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

capture_condensdist_fit_progress_trace <- function(expr,
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

make_condensdist_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 18L
  dat <- data.frame(
    x = seq(-0.9, 0.9, length.out = n),
    y = 0.4 * seq(-0.9, 0.9, length.out = n) + sin(seq(-0.9, 0.9, length.out = n))
  )

  list(
    dat = dat,
    tx = dat["x"],
    ty = dat["y"],
    n = n,
    cdens.bw = npcdensbw(
      xdat = dat["x"],
      ydat = dat["y"],
      bws = c(0.24, 0.24),
      bandwidth.compute = FALSE
    ),
    cdist.bw = npcdistbw(
      xdat = dat["x"],
      ydat = dat["y"],
      bws = c(0.24, 0.24),
      bandwidth.compute = FALSE
    )
  )
}

test_that("npcdens direct bws fit emits single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_condensdist_fit_progress_fixture()
  total <- 2L * fixture$n

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_condensdist_fit_progress_trace(
    npcdens(
      bws = fixture$cdens.bw,
      txdat = fixture$tx,
      tydat = fixture$ty
    ),
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)

  expect_s3_class(actual$value, "condensity")
  expect_true(any(grepl(
    sprintf("^\\[npRmpi\\] Fitting conditional density 1/%d \\([0-9]+\\.[0-9]%%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", total),
    lines
  )))
  expect_true(any(grepl(
    sprintf("^\\[npRmpi\\] Fitting conditional density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", total, total),
    lines
  )))
})

test_that("npcdens bw to fit route hands off immediately into single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_condensdist_fit_progress_fixture()
  total <- 2L * fixture$n

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_condensdist_fit_progress_trace(
    npcdens(
      txdat = fixture$tx,
      tydat = fixture$ty,
      nmulti = 1L
    ),
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting conditional density 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", total),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting conditional density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", total, total),
    lines
  )

  expect_s3_class(actual$value, "condensity")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npcdens nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_condensdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_condensdist_fit_progress_trace(
    npcdens(
      y ~ x,
      data = fixture$dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    ),
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)
  fit.start.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting conditional density 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", fixture$n),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting conditional density %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )
  powell.pos <- grep("^\\[npRmpi\\] Refining bandwidth \\(", lines)
  bandwidth.pos <- grep("^\\[npRmpi\\] Selecting degree and bandwidth", lines)

  expect_s3_class(actual$value, "condensity")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npcdist direct bws fit emits single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_condensdist_fit_progress_fixture()
  total <- 2L * fixture$n

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_condensdist_fit_progress_trace(
    npcdist(
      bws = fixture$cdist.bw,
      txdat = fixture$tx,
      tydat = fixture$ty
    ),
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)

  expect_s3_class(actual$value, "condistribution")
  expect_true(any(grepl(
    sprintf("^\\[npRmpi\\] Fitting conditional distribution 1/%d \\([0-9]+\\.[0-9]%%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$", total),
    lines
  )))
  expect_true(any(grepl(
    sprintf("^\\[npRmpi\\] Fitting conditional distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", total, total),
    lines
  )))
})

test_that("npcdist bw to fit route hands off immediately into single-line fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_condensdist_fit_progress_fixture()
  total <- 2L * fixture$n

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_condensdist_fit_progress_trace(
    npcdist(
      txdat = fixture$tx,
      tydat = fixture$ty,
      nmulti = 1L
    ),
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting conditional distribution 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", total),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting conditional distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", total, total),
    lines
  )

  expect_s3_class(actual$value, "condistribution")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npcdist nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_condensdist_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_condensdist_fit_progress_trace(
    npcdist(
      y ~ x,
      data = fixture$dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    ),
    now = condensdist_fit_progress_time_counter()
  )

  lines <- condensdist_fit_progress_lines(actual)
  fit.start.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting conditional distribution 0/%d \\(0\\.0%%, elapsed 0\\.0s, eta 0\\.0s\\): starting$", fixture$n),
    lines
  )
  fit.finish.pos <- grep(
    sprintf("^\\[npRmpi\\] Fitting conditional distribution %d/%d \\(100\\.0%%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$", fixture$n, fixture$n),
    lines
  )
  powell.pos <- grep("^\\[npRmpi\\] Refining bandwidth \\(", lines)
  bandwidth.pos <- grep("^\\[npRmpi\\] Selecting degree and bandwidth", lines)

  expect_s3_class(actual$value, "condistribution")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})
