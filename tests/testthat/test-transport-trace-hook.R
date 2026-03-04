test_that("transport trace hook writes optional breadcrumb lines", {
  trace_fun <- getFromNamespace(".npRmpi_transport_trace", "npRmpi")
  tf <- tempfile("npRmpi-transport-", fileext = ".tsv")
  old <- getOption("npRmpi.transport.trace.file")
  options(npRmpi.transport.trace.file = tf)
  on.exit(options(npRmpi.transport.trace.file = old), add = TRUE)

  ok <- trace_fun(
    role = "master",
    event = "remote.exec.start",
    fields = list(tag = 42L, ret = TRUE, simplify = TRUE, comm = 1L)
  )
  expect_true(isTRUE(ok))
  expect_true(file.exists(tf))

  lines <- readLines(tf, warn = FALSE)
  expect_true(length(lines) >= 1L)
  expect_true(any(grepl("event=remote.exec.start", lines, fixed = TRUE)))
  expect_true(any(grepl("tag=42", lines, fixed = TRUE)))
})

test_that("transport trace hook is no-op when unset", {
  trace_fun <- getFromNamespace(".npRmpi_transport_trace", "npRmpi")
  old <- getOption("npRmpi.transport.trace.file")
  options(npRmpi.transport.trace.file = "")
  on.exit(options(npRmpi.transport.trace.file = old), add = TRUE)
  ok <- trace_fun(role = "master", event = "noop", fields = list())
  expect_false(isTRUE(ok))
})
