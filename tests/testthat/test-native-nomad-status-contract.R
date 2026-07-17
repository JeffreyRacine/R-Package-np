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
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "subprocess library setup unavailable")

  res <- npRmpi_run_rscript_subprocess(
    c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "state <- getFromNamespace('.npRmpi_mpi_interrupt_state', 'npRmpi')",
      "scope <- getFromNamespace('.npRmpi_mpi_interrupt_scope', 'npRmpi')",
      "with.scope <- getFromNamespace('.npRmpi_with_command_interrupt_scope', 'npRmpi')",
      "stopifnot(!state(3L), !state(2L), !state(1L), !state(4L))",
      "stopifnot(!scope(1L), !scope(2L))",
      "stopifnot(identical(with.scope(42L), 42L))",
      "err <- tryCatch(with.scope(stop('scope error')), error = conditionMessage)",
      "stopifnot(identical(err, 'scope error'))",
      "cat('NO_DISTRIBUTED_POOL_INTERRUPT_STATE_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  expect_identical(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(
    "NO_DISTRIBUTED_POOL_INTERRUPT_STATE_OK",
    res$output,
    fixed = TRUE
  )))
})
