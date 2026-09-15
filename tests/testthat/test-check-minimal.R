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
test_that("copula grid coordinates cannot fall back to equal-length training data", {
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  coordinates <- get(".npcopula_eval_xgrid", ns)
  x <- structure(list(evaluation = "grid", xnames = c("x", "log(y)"),
    copula = c(.2, .4, .6, .8), grid.dim = c(2L, 2L),
    data = data.frame(x = 1:4, check.names = FALSE, "log(y)" = 5:8),
    eval = data.frame(copula = c(.2, .4, .6, .8), u1 = c(.2, .8, .2, .8),
      u2 = c(.2, .2, .8, .8), x = 9:12, log.y. = 13:16)),
    class = "npcopula")
  expect_identical(coordinates(x)[[1L]], 9:12)
  expect_identical(names(coordinates(x)), x$xnames)
  x$evaluation <- "sample"
  expect_identical(coordinates(x)[[1L]], 1:4)
})
