test_that("mpi.setup.rngstream uses explicit .GlobalEnv seed assignment", {
  fn.body <- paste(deparse(body(mpi.setup.rngstream), width.cutoff = 500L), collapse = " ")
  expect_match(fn.body, "assign\\(\"\\.Random.seed\", seed, envir = \\.GlobalEnv\\)")
  expect_false(grepl("<<-", fn.body, fixed = TRUE))
})
