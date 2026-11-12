test_that("autodispatch option sync treats a null snapshot as unsynchronized", {
  keys <- getFromNamespace(".npRmpi_autodispatch_option_keys", "npRmpi")()
  snapshot <- getFromNamespace(".npRmpi_autodispatch_option_snapshot", "npRmpi")(keys)
  should_sync <- getFromNamespace(".npRmpi_autodispatch_should_sync_options", "npRmpi")

  state <- getFromNamespace(".npRmpi_option_agreement", "npRmpi")
  saved <- state$snapshot; on.exit(state$snapshot <- saved, add = TRUE)
  state$snapshot <- NULL
  old <- options(
    npRmpi.autodispatch.option.sync = "onchange"
  )
  on.exit(options(old), add = TRUE)

  expect_true(should_sync(snapshot))
  expect_null(state$snapshot)
  expect_true(should_sync(snapshot))
  # Only a completed writer may publish agreement; inspection is read-only.
  state$snapshot <- snapshot
  expect_identical(state$snapshot, snapshot)
  expect_false(should_sync(snapshot))
})

test_that("npRmpi.init keeps autodispatch option snapshot unprimed", {
  npRmpi.init <- getFromNamespace("npRmpi.init", "npRmpi")
  init_body <- paste(deparse(body(npRmpi.init)), collapse = "\n")
  expect_match(init_body, ".npRmpi_option_agreement$snapshot <- NULL", fixed = TRUE)
  expect_false(grepl("\\.npRmpi_autodispatch_prime_options\\(\\)", init_body))
})


test_that("ordinary options restoration cannot roll back worker agreement", {
  state <- getFromNamespace(".npRmpi_option_agreement", "npRmpi")
  unrelated <- as.list(state); unrelated$snapshot <- NULL
  saved <- state$snapshot; on.exit(state$snapshot <- saved, add = TRUE)
  sync <- getFromNamespace(".npRmpi_autodispatch_should_sync_options", "npRmpi")
  old <- options(); on.exit(options(old), add = TRUE)
  options(npRmpi.autodispatch.option.sync="onchange")
  state$snapshot <- NULL
  expect_true(sync(list(np.tree=FALSE)))
  state$snapshot <- list(np.tree=FALSE)
  op <- options()
  expect_true(sync(list(np.tree=TRUE)))
  expect_identical(state$snapshot,list(np.tree=FALSE))
  state$snapshot <- list(np.tree=TRUE)
  options(op)
  expect_identical(state$snapshot,list(np.tree=TRUE))
  expect_true(sync(list(np.tree=FALSE)))
  state$snapshot <- list(np.tree=FALSE)
  expect_false(sync(list(np.tree=FALSE)))
  expect_true("snapshot" %in% ls(state))
  remaining <- as.list(state); remaining$snapshot <- NULL
  expect_identical(remaining, unrelated)
})
