expect_npregbw_powell_progress_surface <- function(lines) {
  powell.lines <- lines[grepl("^\\[(np|npRmpi)\\] Refining bandwidth \\(", lines)]
  info <- paste(lines, collapse = "\n")

  expect_true(length(powell.lines) > 0L, info = info)
  expect_false(any(grepl(
    "^\\[(np|npRmpi)\\] Refining NOMAD solution with one Powell hot start at degree ",
    lines
  )), info = info)
  expect_false(any(grepl("best (", powell.lines, fixed = TRUE)), info = info)
  powell.iter.lines <- powell.lines[grepl("iter [0-9]+", powell.lines)]
  expect_true(length(powell.iter.lines) > 0L, info = info)
  expect_true(any(grepl("^\\[(np|npRmpi)\\] Refining bandwidth \\(elapsed ", powell.iter.lines)), info = info)
  expect_true(any(grepl(", degree \\(", powell.iter.lines)), info = info)
}

test_that("npregbw NOMAD plus Powell progress keeps lines compact and restart-oriented", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(
    np.messages = TRUE,
    np.tree = FALSE,
    npRmpi.autodispatch = TRUE,
    np.progress.start.grace.known.sec = 0,
    np.progress.start.grace.unknown.sec = 0,
    np.progress.interval.known.sec = 0,
    np.progress.interval.unknown.sec = 0
  )
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260321)
  dat <- data.frame(x = sort(runif(16)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  msgs <- with_npRmpi_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = degree_progress_time_values(seq(0, 40, by = 0.25))
    ),
    capture_degree_messages_only(
      get("npregbw", envir = asNamespace("npRmpi"), inherits = FALSE)(
        y ~ x,
        data = dat,
        nomad = TRUE,
        degree.max = 1L,
        nmulti = 2L,
        max.bb.eval = 8
      )
    )
  )

  expect_false(any(grepl("^\\[(np|npRmpi)\\] Starting NOMAD automatic polynomial degree search at degree ", msgs)))
  expect_false(any(grepl("^\\[(np|npRmpi)\\] Refining NOMAD solution with one Powell hot start at degree ", msgs)))
  expect_false(any(grepl("starting at degree \\(", msgs)))
  expect_false(any(grepl("nomad\\+powell", msgs, ignore.case = TRUE)))
  expect_false(any(grepl("eval [0-9]+", msgs)))
  expect_false(any(grepl("fval=", msgs, fixed = TRUE)))
  expect_false(any(grepl("%|eta ", msgs)))
  expect_true(any(grepl("^\\[(np|npRmpi)\\] Selecting degree and bandwidth \\(", msgs)))
  expect_true(any(grepl("^\\[(np|npRmpi)\\] Refining bandwidth \\(", msgs)))
  expect_npregbw_powell_progress_surface(msgs)
  expect_true(any(grepl("multistart [12]/2", msgs)))
  expect_true(any(grepl("iteration [0-9]+", msgs)))
  expect_true(any(grepl("deg \\(", msgs)))
  expect_true(any(grepl("best \\(", msgs)))
})
