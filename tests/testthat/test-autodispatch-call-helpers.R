test_that(".npRmpi_bcast_cmd_expr forwards command expression structurally", {
  env <- new.env(parent = baseenv())
  env$seen <- NULL
  env$mpi.bcast.cmd <- function(cmd, comm = 1L, caller.execute = TRUE) {
    env$seen <- list(expr = substitute(cmd),
                     value = cmd,
                     comm = comm,
                     caller.execute = caller.execute)
    "OK"
  }

  out <- evalq(.npRmpi_bcast_cmd_expr(quote(x <- 1L), comm = 3L, caller.execute = FALSE), envir = env)

  expect_identical(out, "OK")
  expect_true(is.list(env$seen))
  expect_true(is.call(env$seen$expr))
  expect_identical(env$seen$expr, quote(x <- 1L))
  expect_identical(env$seen$value, 1L)
  expect_identical(env$seen$comm, 3L)
  expect_identical(env$seen$caller.execute, FALSE)
})

test_that(".npRmpi_autodispatch_call delegates to shared distributed-call helper", {
  fn.body <- paste(deparse(body(.npRmpi_autodispatch_call), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.npRmpi_distributed_call_impl\\(mc = mc, caller_env = caller_env, comm = comm, warn_nested = TRUE\\)")
})

test_that(".npRmpi_manual_distributed_call delegates to shared distributed-call helper", {
  fn.body <- paste(deparse(body(.npRmpi_manual_distributed_call), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "\\.npRmpi_distributed_call_impl\\(mc = mc, caller_env = caller_env, comm = comm, warn_nested = FALSE\\)")
})

test_that("autodispatch eval helpers route through shared command executor", {
  bcast.body <- paste(deparse(body(.npRmpi_bcast_robj_by_name), width.cutoff = 500L), collapse = " ")
  nodisp.body <- paste(deparse(body(.npRmpi_eval_without_dispatch), width.cutoff = 500L), collapse = " ")
  arg.body <- paste(deparse(body(.npRmpi_autodispatch_eval_arg), width.cutoff = 500L), collapse = " ")

  expect_match(bcast.body, "\\.npRmpi_eval_scmd\\(call, envir = caller_env\\)")
  expect_match(nodisp.body, "\\.npRmpi_eval_scmd\\(mc\\.eval, envir = caller_env\\)")
  expect_match(arg.body, "\\.npRmpi_eval_scmd\\(expr, envir = caller_env\\)")
})
