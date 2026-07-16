test_that("native NOMAD status decoding preserves R condition semantics", {
  decode <- getFromNamespace(".np_nomad_native_status", "npRmpi")

  ok <- list(status = 0L, result_status = 0L, message = "")
  expect_invisible(decode(ok, "native test route"))

  expect_error(
    decode(
      list(status = 1L, result_status = 0L, message = "ordinary failure"),
      "native test route"
    ),
    "native test route failed.*ordinary failure"
  )
  expect_error(
    decode(list(status = 0L, result_status = NA_integer_), "native test route"),
    "native test route failed"
  )

  interrupted <- tryCatch(
    decode(
      list(status = 4L, result_status = 4L, message = "user interrupt"),
      "native test route"
    ),
    interrupt = identity
  )
  expect_s3_class(interrupted, "interrupt")
  expect_identical(conditionMessage(interrupted), "user interrupt")
})

test_that("MPI interrupt state is inactive without a distributed pool", {
  state <- getFromNamespace(".npRmpi_mpi_interrupt_state", "npRmpi")
  scope <- getFromNamespace(".npRmpi_mpi_interrupt_scope", "npRmpi")
  with.scope <- getFromNamespace(
    ".npRmpi_with_command_interrupt_scope",
    "npRmpi"
  )
  expect_false(state(3L))
  expect_false(state(2L))
  expect_false(state(1L))
  expect_false(state(4L))
  expect_false(scope(1L))
  expect_false(scope(2L))
  expect_identical(with.scope(42L), 42L)
  expect_error(with.scope(stop("scope error")), "scope error")
})
