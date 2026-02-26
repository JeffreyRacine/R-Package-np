test_that("bcast command resolver resolves internal worker functions", {
  ns <- asNamespace("npRmpi")
  funref <- get(".npRmpi_bcast_cmd_funref", envir = ns, inherits = FALSE)

  w.apply <- get(".mpi.worker.apply", envir = ns, inherits = FALSE)
  w.applylb <- get(".mpi.worker.applyLB", envir = ns, inherits = FALSE)

  got.apply <- funref(as.name(".mpi.worker.apply"))
  got.applylb <- funref(as.name(".mpi.worker.applyLB"))

  expect_true(is.function(got.apply))
  expect_true(is.function(got.applylb))
  expect_identical(got.apply, w.apply)
  expect_identical(got.applylb, w.applylb)
})

