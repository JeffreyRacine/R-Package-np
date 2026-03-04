test_that(".npRmpi_eval_scmd handles symbol, call, and argument paths", {
  env <- new.env(parent = baseenv())
  env$x <- 5L
  env$plus1 <- function(v) v + 1L

  expect_identical(.npRmpi_eval_scmd(quote(x), envir = env), 5L)
  expect_identical(.npRmpi_eval_scmd(quote(plus1), arg = list(v = 2L), envir = env), 3L)
  expect_identical(.npRmpi_eval_scmd(quote({ y <- 9L; y }), envir = env), 9L)
})

test_that(".npRmpi_eval_scmd preserves NULL symbol values", {
  env <- new.env(parent = baseenv())
  env$x <- NULL
  expect_null(.npRmpi_eval_scmd(quote(x), envir = env))
})

test_that(".npRmpi_eval_scmd evaluates through shared strict helper", {
  fn.body <- paste(deparse(body(.npRmpi_eval_scmd), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.np_try_eval_in_frames\\(scmd, eval_env = envir, search_frames = FALSE\\)")
})

test_that("session attach loop delegates to shared worker loop helper", {
  fn.body <- paste(deparse(body(.npRmpi_session_attach_worker_loop), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.npRmpi_worker_loop\\(")
  expect_match(fn.body, "loop\\.label = \"attach slave\"")
})

test_that("mpi.bcast.cmd caller path executes no-arg commands through helper", {
  fn.body <- paste(deparse(body(mpi.bcast.cmd), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.npRmpi_eval_scmd\\(tcmd, envir = parent\\.frame\\(\\)\\)")
})

test_that(".mpi.worker.exec no longer uses .mpi.err side-channel", {
  fn.body <- paste(deparse(body(.mpi.worker.exec), width.cutoff = 500L), collapse = " ")
  expect_no_match(fn.body, "\\.mpi\\.err")
  expect_no_match(fn.body, "myerrcode")
  expect_match(fn.body, "type <- \\.typeindex\\(out\\)")
})

test_that("applyLB short-circuit paths use shared cache-clear helper", {
  body.lb <- paste(deparse(body(mpi.applyLB), width.cutoff = 500L), collapse = " ")
  body.ilb <- paste(deparse(body(mpi.iapplyLB), width.cutoff = 500L), collapse = " ")

  expect_match(body.lb, "\\.npRmpi_clear_applylb_cache\\(\\)")
  expect_match(body.ilb, "\\.npRmpi_clear_applylb_cache\\(\\)")
})
