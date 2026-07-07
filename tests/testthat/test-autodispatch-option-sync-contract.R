test_that("autodispatch option sync treats a null snapshot as unsynchronized", {
  keys <- getFromNamespace(".npRmpi_autodispatch_option_keys", "npRmpi")()
  snapshot <- getFromNamespace(".npRmpi_autodispatch_option_snapshot", "npRmpi")(keys)
  should_sync <- getFromNamespace(".npRmpi_autodispatch_should_sync_options", "npRmpi")

  old <- options(
    npRmpi.autodispatch.option.sync = "onchange",
    npRmpi.autodispatch.option.snapshot = NULL
  )
  on.exit(options(old), add = TRUE)

  expect_true(should_sync(snapshot))
  expect_identical(getOption("npRmpi.autodispatch.option.snapshot"), snapshot)
  expect_false(should_sync(snapshot))
})

test_that("npRmpi.init keeps autodispatch option snapshot unprimed", {
  npRmpi.init <- getFromNamespace("npRmpi.init", "npRmpi")
  init_body <- paste(deparse(body(npRmpi.init)), collapse = "\n")
  expect_match(init_body, "npRmpi\\.autodispatch\\.option\\.snapshot = NULL")
  expect_false(grepl("\\.npRmpi_autodispatch_prime_options\\(\\)", init_body))
})
