test_that("within-task completion is scoped and credited exactly once", {
  ns <- asNamespace("npRmpi")
  evaluate <- get(".npRmpi_bootstrap_task_evaluate", ns)
  checkpoint <- get(".npRmpi_bootstrap_task_checkpoint", ns)
  runtime <- get(".npRmpi_bootstrap_task_progress", ns)
  old <- runtime$checkpoint
  on.exit(runtime$checkpoint <- old)
  sentinel <- function(done) stop("outer callback used")
  runtime$checkpoint <- sentinel
  seen <- new.env(parent = emptyenv())
  seen$deltas <- integer()
  report <- function(delta) seen$deltas <- c(seen$deltas, delta)
  ans <- evaluate({
    checkpoint(8L)
    checkpoint(8L)
    checkpoint(16L)
    42
  }, 19L, report)
  expect_identical(ans, 42)
  expect_identical(seen$deltas, c(8L, 8L, 3L))
  expect_identical(runtime$checkpoint, sentinel)
  for (bad in list(-1L, 20L, NA_integer_, Inf, 1.5, c(1L, 2L))) {
    seen$deltas <- integer()
    expect_error(evaluate(checkpoint(bad), 19L, report),
                 "invalid internal bootstrap completion")
    expect_identical(seen$deltas, integer())
    expect_identical(runtime$checkpoint, sentinel)
  }
  seen$deltas <- integer()
  expect_error(evaluate({checkpoint(8L); checkpoint(7L)}, 19L, report),
               "invalid internal bootstrap completion")
  expect_identical(seen$deltas, 8L)
  seen$deltas <- integer()
  expect_error(evaluate({checkpoint(8L); stop("numerical failure")}, 19L, report),
               "numerical failure")
  expect_identical(seen$deltas, 8L)
  expect_identical(runtime$checkpoint, sentinel)
  seen$deltas <- integer()
  failed <- structure("failed", class = "try-error")
  expect_identical(evaluate(failed, 19L, report), failed)
  expect_identical(seen$deltas, integer())
  expect_identical(evaluate(41, 19L), 41)
  expect_identical(runtime$checkpoint, sentinel)
})

test_that("fanout rethrows the original error or message-less interrupt", {
  rethrow <- getFromNamespace(".npRmpi_fanout_rethrow", "npRmpi")
  original <- simpleError("original transport condition")
  caught <- tryCatch(rethrow(original), error = identity)
  expect_identical(caught, original)
  interruption <- structure(list(), class = c("interrupt", "condition"))
  caught <- tryCatch(rethrow(interruption), error = identity, interrupt = identity)
  expect_identical(caught, interruption)
  expect_false(inherits(caught, "error"))
})
