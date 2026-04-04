npplreg_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

npplreg_fit_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

npplreg_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

capture_npplreg_fit_progress_trace <- function(expr,
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

make_npplreg_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 32L
  dat <- data.frame(
    x = seq(-0.9, 0.9, length.out = n),
    z = seq(0.2, 1.8, length.out = n)
  )
  dat$y <- 1.5 * dat$x + sin(2 * pi * dat$z) + 0.1 * cos(pi * dat$x)

  list(
    dat = dat,
    tx = dat["x"],
    tz = dat["z"],
    y = dat$y,
    bw = npplregbw(
      xdat = dat["x"],
      zdat = dat["z"],
      ydat = dat$y,
      bws = matrix(c(0.35, 0.35), nrow = 2L),
      bandwidth.compute = FALSE
    )
  )
}

test_that("npplreg fit progress uses the bandwidth single-line surface", {
  begin <- getFromNamespace(".np_plreg_fit_progress_begin", "npRmpi")

  old_opts <- options(np.messages = TRUE)
  on.exit(options(old_opts), add = TRUE)

  state <- with_nprmpi_progress_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_now = function() 0
    ),
    begin(xnames = "x")
  )

  expect_identical(state$surface, "bandwidth")
  expect_identical(state$renderer, "single_line")
  expect_equal(state$total, 2L)
})

test_that("npplreg direct bws fit emits known-total fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npplreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npplreg_fit_progress_trace(
    npplreg(
      bws = fixture$bw,
      txdat = fixture$tx,
      tzdat = fixture$tz,
      tydat = fixture$y
    ),
    now = npplreg_fit_progress_time_counter()
  )

  lines <- npplreg_fit_progress_lines(actual)

  expect_s3_class(actual$value, "plregression")
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting partially linear regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): y~z$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting partially linear regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\): x~z$",
    lines
  )))
})

test_that("npplreg direct bws fit stays silent below start grace without handoff", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npplreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0.75,
    np.progress.interval.known.sec = 0.5
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npplreg_fit_progress_trace(
    npplreg(
      bws = fixture$bw,
      txdat = fixture$tx,
      tzdat = fixture$tz,
      tydat = fixture$y
    ),
    now = npplreg_fit_progress_time_values(c(0, 0.2, 0.4, 0.6))
  )

  expect_length(actual$trace, 0L)
})

test_that("npplreg formula bw to fit route hands off immediately into fit progress", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npplreg_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npplreg_fit_progress_trace(
    npplreg(
      y ~ x | z,
      data = fixture$dat,
      nmulti = 1L
    ),
    now = npplreg_fit_progress_time_counter()
  )

  lines <- npplreg_fit_progress_lines(actual)
  fit.zero.pos <- grep(
    "^\\[npRmpi\\] Fitting partially linear regression 0/2 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting y~z$",
    lines
  )
  fit.one.pos <- grep(
    "^\\[npRmpi\\] Fitting partially linear regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): y~z$",
    lines
  )
  fit.two.pos <- grep(
    "^\\[npRmpi\\] Fitting partially linear regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\): x~z$",
    lines
  )
  bandwidth.pos <- grep("^\\[npRmpi\\] Bandwidth selection \\(", lines)

  expect_s3_class(actual$value, "plregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.zero.pos) == 1L)
  expect_true(length(fit.one.pos) >= 1L)
  expect_true(length(fit.two.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.zero.pos[[1L]])
  expect_lt(fit.zero.pos[[1L]], fit.one.pos[[1L]])
  expect_lt(fit.one.pos[[1L]], fit.two.pos[[1L]])
})

test_that("npplreg nomad to powell to fit route preserves single-line fit handoff", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260404)
  n <- 18L
  dat <- data.frame(x = rnorm(n), z = sort(runif(n)))
  dat$y <- 1 + 0.5 * dat$x + sin(2 * pi * dat$z) + rnorm(n, sd = 0.05)

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_npplreg_fit_progress_trace(
    npplreg(
      y ~ x | z,
      data = dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    ),
    now = npplreg_fit_progress_time_counter()
  )

  lines <- npplreg_fit_progress_lines(actual)
  fit.zero.pos <- grep(
    "^\\[npRmpi\\] Fitting partially linear regression 0/2 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting y~z$",
    lines
  )
  fit.one.pos <- grep(
    "^\\[npRmpi\\] Fitting partially linear regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\): y~z$",
    lines
  )
  fit.two.pos <- grep(
    "^\\[npRmpi\\] Fitting partially linear regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\): x~z$",
    lines
  )
  powell.pos <- grep("^\\[npRmpi\\] Refining bandwidth \\(", lines)
  bandwidth.pos <- grep("^\\[npRmpi\\] Selecting degree and bandwidth", lines)

  expect_s3_class(actual$value, "plregression")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.zero.pos) == 1L)
  expect_true(length(fit.one.pos) >= 1L)
  expect_true(length(fit.two.pos) >= 1L)
  expect_lt(max(powell.pos), fit.zero.pos[[1L]])
  expect_lt(max(bandwidth.pos), fit.zero.pos[[1L]])
  expect_lt(fit.zero.pos[[1L]], fit.one.pos[[1L]])
  expect_lt(fit.one.pos[[1L]], fit.two.pos[[1L]])
})
