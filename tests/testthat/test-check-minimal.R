test_that("check-minimal verifies native registration and package metadata", {
  dll <- getLoadedDLLs()[["npRmpi"]]
  regs <- getDLLRegisteredRoutines(dll)
  dll_info <- unclass(dll)

  expect_false(isTRUE(dll_info[["dynamicLookup"]]))
  expect_true(all(c("mpi_finalize", "mpi_gather", "mpi_sendrecv", "mpi_wait") %in% names(regs$.Call)))

  desc <- utils::packageDescription("npRmpi")
  expect_identical(desc$SystemRequirements, "MPI")
})

test_that("check-minimal verifies the load hook helper contract", {
  ns <- asNamespace("npRmpi")
  load_hook <- get(".onLoad", envir = ns, inherits = FALSE)
  try_dynload <- get(".npRmpi_try_dynload", envir = ns, inherits = FALSE)
  load_body <- paste(deparse(body(load_hook), width.cutoff = 500L), collapse = " ")
  expect_match(load_body, "\\.npRmpi_try_dynload\\(lib = lib, pkg = pkg\\)")
  expect_false(try_dynload(tempdir(), "npRmpi_definitely_missing_pkg"))
})

test_that("check-minimal verifies deterministic bandwidth name metadata updates", {
  bws <- list(varnames = list(x = "oldx", y = "oldy"))
  names_in <- list(xnames = "x_new", ynames = "y_new")

  out <- updateBwNameMetadata(nameList = names_in, bws = bws)

  expect_equal(out$xnames, "x_new")
  expect_equal(out$ynames, "y_new")
  expect_equal(out$varnames$x, "x_new")
  expect_equal(out$varnames$y, "y_new")
  expect_equal(bws$varnames$x, "oldx")
  expect_equal(bws$varnames$y, "oldy")
})
