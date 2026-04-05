npindex_fit_progress_time_counter <- function(start = 0, by = 1.7) {
  current <- start
  function() {
    current <<- current + by
    current
  }
}

npindex_fit_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

npindex_fit_progress_lines <- function(shadow) {
  vapply(shadow$trace, `[[`, character(1L), "line")
}

make_npindex_fit_progress_fixture <- function() {
  set.seed(20260404)
  n <- 24L
  dat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  index <- dat$x1 + 0.6 * dat$x2
  dat$y <- sin(index) + 0.2 * index^2 + rnorm(n, sd = 0.05)

  list(
    dat = dat,
    tx = dat[c("x1", "x2")],
    y = dat$y,
    bw = npindexbw(
      xdat = dat[c("x1", "x2")],
      ydat = dat$y,
      bws = c(1, 0.6, 0.35),
      method = "ichimura",
      regtype = "lp",
      degree = 1L,
      bernstein.basis = TRUE,
      bwtype = "fixed",
      bandwidth.compute = FALSE
    ),
    bw_lc_fixed = npindexbw(
      xdat = dat[c("x1", "x2")],
      ydat = dat$y,
      bws = c(1, 0.6, 0.35),
      method = "ichimura",
      regtype = "lc",
      bwtype = "fixed",
      bandwidth.compute = FALSE
    )
  )
}

test_that("npindex direct lp fit emits fit progress without handoff", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npindex(
      bws = fixture$bw,
      txdat = fixture$tx,
      tydat = fixture$y
    ),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_counter()
  )

  lines <- npindex_fit_progress_lines(actual)

  expect_s3_class(actual$value, "singleindex")
  expect_false(any(grepl(": starting$", lines)))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 1/24 \\([0-9]+\\.[0-9]%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})

test_that("npindex direct lp fit stays silent below start grace without handoff", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0.75,
    np.progress.interval.known.sec = 0.5
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npindex(
      bws = fixture$bw,
      txdat = fixture$tx,
      tydat = fixture$y
    ),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_values(c(0, 0.2, 0.4, 0.6))
  )

  expect_length(actual$trace, 0L)
})

test_that("npindex lp bw to fit route hands off into the regression fit surface", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npindex(
      y ~ x1 + x2,
      data = fixture$dat,
      method = "ichimura",
      regtype = "lp",
      degree.select = "coordinate",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_counter()
  )

  lines <- npindex_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    lines
  )
  fit.finish.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )

  expect_s3_class(actual$value, "singleindex")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("npindex lp nomad to powell to fit route preserves single-line handoff", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npindex(
      y ~ x1 + x2,
      data = fixture$dat,
      method = "ichimura",
      nomad = TRUE,
      degree.min = 0L,
      degree.max = 1L,
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_counter()
  )

  lines <- npindex_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Selecting degree and bandwidth", lines)
  powell.pos <- grep("^\\[npRmpi\\] Refining bandwidth \\(", lines)
  fit.start.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    lines
  )
  fit.finish.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )

  expect_s3_class(actual$value, "singleindex")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(powell.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(max(powell.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("predict.singleindex lp re-entry emits evaluation and training fit progress without handoff", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  fit <- npindex(
    bws = fixture$bw,
    txdat = fixture$tx,
    tydat = fixture$y
  )

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    predict(fit, newdata = fixture$dat[c(2L, 7L), c("x1", "x2")]),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_counter()
  )

  lines <- npindex_fit_progress_lines(actual)

  expect_false(any(grepl(": starting$", lines)))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 1/24 \\([0-9]+\\.[0-9]%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})

test_that("npindex direct fixed lc fit emits fit progress without handoff", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npindex(
      bws = fixture$bw_lc_fixed,
      txdat = fixture$tx,
      tydat = fixture$y
    ),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_counter()
  )

  lines <- npindex_fit_progress_lines(actual)

  expect_s3_class(actual$value, "singleindex")
  expect_false(any(grepl(": starting$", lines)))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 1/24 \\([0-9]+\\.[0-9]%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})

test_that("npindex direct fixed lc fit stays silent below start grace without handoff", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0.75,
    np.progress.interval.known.sec = 0.5
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npindex(
      bws = fixture$bw_lc_fixed,
      txdat = fixture$tx,
      tydat = fixture$y
    ),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_values(c(0, 0.2, 0.4, 0.6))
  )

  expect_length(actual$trace, 0L)
})

test_that("npindex fixed lc bw to fit route hands off into the regression fit surface", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    npindex(
      y ~ x1 + x2,
      data = fixture$dat,
      method = "ichimura",
      regtype = "lc",
      bwtype = "fixed",
      nmulti = 1L
    ),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_counter()
  )

  lines <- npindex_fit_progress_lines(actual)
  bandwidth.pos <- grep("^\\[npRmpi\\] Bandwidth selection \\(", lines)
  fit.start.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 0/24 \\(0\\.0%, elapsed 0\\.0s, eta 0\\.0s\\): starting$",
    lines
  )
  fit.finish.pos <- grep(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )

  expect_s3_class(actual$value, "singleindex")
  expect_true(length(bandwidth.pos) > 0L)
  expect_true(length(fit.start.pos) == 1L)
  expect_true(length(fit.finish.pos) >= 1L)
  expect_lt(max(bandwidth.pos), fit.start.pos[[1L]])
  expect_lt(fit.start.pos[[1L]], fit.finish.pos[[1L]])
})

test_that("predict.singleindex fixed lc re-entry emits evaluation and training fit progress without handoff", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_npindex_fit_progress_fixture()

  fit <- npindex(
    bws = fixture$bw_lc_fixed,
    txdat = fixture$tx,
    tydat = fixture$y
  )

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    np.progress.start.grace.known.sec = 0,
    np.progress.interval.known.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  actual <- capture_progress_shadow_trace(
    predict(fit, newdata = fixture$dat[c(2L, 7L), c("x1", "x2")]),
    force_renderer = "single_line",
    now = npindex_fit_progress_time_counter()
  )

  lines <- npindex_fit_progress_lines(actual)

  expect_false(any(grepl(": starting$", lines)))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 1/2 \\(50\\.0%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 2/2 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 1/24 \\([0-9]+\\.[0-9]%, elapsed [0-9]+\\.[0-9]s, eta [0-9]+\\.[0-9]s\\)$",
    lines
  )))
  expect_true(any(grepl(
    "^\\[npRmpi\\] Fitting regression 24/24 \\(100\\.0%, elapsed [0-9]+\\.[0-9]s, eta 0\\.0s\\)$",
    lines
  )))
})
