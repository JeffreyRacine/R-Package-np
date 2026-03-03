run_spmd_subprocess <- function(lines, timeout = 60L, env = character()) {
  script <- tempfile("npRmpi-spmd-", fileext = ".R")
  writeLines(lines, script, useBytes = TRUE)
  on.exit(unlink(script), add = TRUE)

  cmd <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(
    cmd,
    c("--no-save", script),
    stdout = TRUE,
    stderr = TRUE,
    timeout = timeout,
    env = env
  ))
  status <- attr(out, "status")
  if (is.null(status))
    status <- 0L
  list(status = as.integer(status), output = out)
}

test_that("SPMD envelope sequence contract fails fast on mismatch", {
  make.env <- getFromNamespace(".npRmpi_spmd_make_envelope", "npRmpi")
  exec.env <- getFromNamespace(".npRmpi_spmd_execute_local", "npRmpi")
  try.exec <- getFromNamespace(".npRmpi_spmd_try_execute_local", "npRmpi")
  set.seq <- getFromNamespace(".npRmpi_spmd_seq_set", "npRmpi")

  old.seq <- getOption("npRmpi.spmd.seq_id", 0L)
  on.exit(options(npRmpi.spmd.seq_id = old.seq), add = TRUE)

  set.seq(0L)
  env1 <- make.env(opcode = "spmd.ping", timeout_class = "smoke")
  ok1 <- exec.env(env1, payload = list(label = "first"))
  expect_true(is.list(ok1) && isTRUE(ok1$ok))
  expect_identical(ok1$ack$status, "ACK")
  expect_identical(ok1$ack$seq_id, 1L)

  env.bad <- make.env(opcode = "spmd.ping", seq_id = 3L, timeout_class = "smoke")
  expect_error(
    exec.env(env.bad, payload = list(label = "bad")),
    "sequence mismatch"
  )

  err <- try.exec(env.bad, payload = list(label = "bad"))
  expect_true(is.list(err) && !isTRUE(err$ok))
  expect_identical(err$ack$status, "ERR")
  expect_match(err$error, "sequence mismatch")
})

test_that("SPMD tiny smoke opcode runs in session mode subprocess", {
  skip_on_cran()
  res <- run_spmd_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves=1, quiet=TRUE)",
      "on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE)",
      "f <- getFromNamespace('.npRmpi_spmd_tiny_smoke', 'npRmpi')",
      "ans <- f(label='session-smoke')",
      "stopifnot(is.list(ans), isTRUE(ans$ok), identical(ans$ack$status, 'ACK'))",
      "cat('SPMD_TINY_SESSION_OK\\n')"
    ),
    timeout = 90L
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("SPMD_TINY_SESSION_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
