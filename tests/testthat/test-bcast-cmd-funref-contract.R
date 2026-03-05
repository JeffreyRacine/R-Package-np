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

test_that("bcast plot-call detector identifies plot calls", {
  detect.plot <- getFromNamespace(".npRmpi_bcast_cmd_is_plot_call", "npRmpi")

  expect_true(isTRUE(detect.plot(quote(plot(x, y)))))
  expect_true(isTRUE(detect.plot(quote(graphics::plot(x, y)))))
  expect_false(isTRUE(detect.plot(quote(npregbw(y ~ x)))))
})

test_that("mpi.bcast.cmd rejects nested plot calls in canonical SPMD mode", {
  bcast.cmd <- getFromNamespace("mpi.bcast.cmd", "npRmpi")

  local_mocked_bindings(
    mpi.comm.rank = function(comm = 1) 0L,
    .package = "npRmpi"
  )

  expect_error(
    bcast.cmd(plot(1, 1), rank = 0, comm = 1, caller.execute = FALSE),
    "plot\\(\\.\\.\\.\\) inside mpi\\.bcast\\.cmd\\(\\.\\.\\.\\) is unsupported",
    fixed = FALSE
  )
})
