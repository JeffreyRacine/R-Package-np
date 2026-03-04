test_that("bootstrap phase mark writes optional trace breadcrumbs", {
  mark <- getFromNamespace(".npRmpi_bootstrap_phase_mark", "npRmpi")
  tf <- tempfile("npRmpi-phase-", fileext = ".tsv")
  old <- getOption("npRmpi.bootstrap.phase.file")
  options(npRmpi.bootstrap.phase.file = tf)
  on.exit(options(npRmpi.bootstrap.phase.file = old), add = TRUE)

  ok <- mark(what = "wild", phase = "preflight", where = "mpi.applyLB:wild")
  expect_true(isTRUE(ok))
  expect_true(file.exists(tf))

  lines <- readLines(tf, warn = FALSE)
  expect_true(length(lines) >= 1L)
  expect_true(any(grepl("what=wild", lines, fixed = TRUE)))
  expect_true(any(grepl("phase=preflight", lines, fixed = TRUE)))
})

test_that("bootstrap phase mark is no-op for empty phase", {
  mark <- getFromNamespace(".npRmpi_bootstrap_phase_mark", "npRmpi")
  tf <- tempfile("npRmpi-phase-empty-", fileext = ".tsv")
  old <- getOption("npRmpi.bootstrap.phase.file")
  options(npRmpi.bootstrap.phase.file = tf)
  on.exit(options(npRmpi.bootstrap.phase.file = old), add = TRUE)

  ok <- mark(what = "wild", phase = "", where = "none")
  expect_false(isTRUE(ok))
  expect_false(file.exists(tf))
})

test_that("bootstrap transport trace writes optional breadcrumbs", {
  trace_fun <- getFromNamespace(".npRmpi_bootstrap_transport_trace", "npRmpi")
  tf <- tempfile("npRmpi-boot-transport-", fileext = ".tsv")
  old <- getOption("npRmpi.bootstrap.transport.trace.file")
  options(npRmpi.bootstrap.transport.trace.file = tf)
  on.exit(options(npRmpi.bootstrap.transport.trace.file = old), add = TRUE)

  ok <- trace_fun(
    what = "wild",
    event = "fanout.start",
    fields = list(workers = 1L, tasks = 4L, master_local = TRUE)
  )
  expect_true(isTRUE(ok))
  expect_true(file.exists(tf))
  lines <- readLines(tf, warn = FALSE)
  expect_true(any(grepl("event=fanout.start", lines, fixed = TRUE)))
})
