test_that(".npRmpi_eval_scmd handles symbol, call, and argument paths", {
  env <- new.env(parent = baseenv())
  env$x <- 5L
  env$plus1 <- function(v) v + 1L

  expect_identical(.npRmpi_eval_scmd(quote(x), envir = env), 5L)
  expect_identical(.npRmpi_eval_scmd(quote(plus1), arg = list(v = 2L), envir = env), 3L)
  expect_identical(.npRmpi_eval_scmd(quote({ y <- 9L; y }), envir = env), 9L)
})

test_that("session attach loop executes messages through helper", {
  fn.body <- paste(deparse(body(.npRmpi_session_attach_worker_loop), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.npRmpi_eval_scmd\\(msg, envir = \\.GlobalEnv\\)")
})
