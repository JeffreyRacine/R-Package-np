test_that("SPMD sequence follows the session despite options restoration", {
  get.seq <- getFromNamespace(".npRmpi_spmd_seq_get", "npRmpi")
  set.seq <- getFromNamespace(".npRmpi_spmd_seq_set", "npRmpi")
  state <- getFromNamespace(".npRmpi_spmd_state", "npRmpi")
  lifecycle <- state$lifecycle
  saved <- get.seq(); on.exit(set.seq(saved), add=TRUE)
  old <- options(); on.exit(options(old), add=TRUE)
  set.seq(12L); op <- options(); set.seq(13L); options(op)
  expect_identical(get.seq(),13L)
  expect_error(set.seq(-1L),"non-negative integer",fixed=TRUE)
  expect_identical(get.seq(),13L)
  expect_identical(ls(state), c("lifecycle", "seq_id"))
  expect_identical(state$lifecycle, lifecycle)
})
