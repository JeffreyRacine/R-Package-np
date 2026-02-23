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
