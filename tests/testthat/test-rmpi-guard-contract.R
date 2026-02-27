test_that("runtime guard blocks execution when package:Rmpi is attached", {
  added.fake <- FALSE
  if (!("package:Rmpi" %in% search())) {
    attach(list(), name = "package:Rmpi")
    added.fake <- TRUE
  }
  on.exit({
    if (added.fake && ("package:Rmpi" %in% search()))
      detach("package:Rmpi")
  }, add = TRUE)

  expect_error(
    npRmpi:::.npRmpi_abort_if_rmpi_attached(where = "unit-test"),
    "package 'Rmpi' is attached"
  )
  expect_error(
    npRmpi:::.npRmpi_require_active_slave_pool(where = "npreg\\(\\)"),
    "package 'Rmpi' is attached"
  )
  expect_error(
    npRmpi:::.npRmpi_autodispatch_preflight(),
    "package 'Rmpi' is attached"
  )
})
