test_that("namespace hygiene scanner distinguishes defects from allowed patterns", {
  mkroot <- function() {
    root <- tempfile("npRmpi-namespace-hygiene-")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    writeLines(c("Package: npRmpi", "Version: 0.0.0"), file.path(root, "DESCRIPTION"))
    root
  }

  good <- mkroot()
  dir.create(file.path(good, "R"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "f <- function(x) stats::na.omit(x)",
      "msg <- 'parser literal ::: should not fail the gate'"
    ),
    file.path(good, "R", "good.R")
  )
  good.scan <- npRmpi_namespace_hygiene_scan(good)
  expect_identical(nrow(good.scan$runtime_same_package), 0L)
  expect_identical(nrow(good.scan$runtime_external_triple), 0L)
  expect_identical(nrow(good.scan$runtime_external_double), 1L)
  expect_identical(nrow(good.scan$runtime_parser_literals), 1L)

  bad.same <- mkroot()
  dir.create(file.path(bad.same, "inst"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    "npRmpi:::.npRmpi_worker_loop(comm = 1)",
    file.path(bad.same, "inst", "Rprofile")
  )
  bad.same.scan <- npRmpi_namespace_hygiene_scan(bad.same)
  expect_identical(nrow(bad.same.scan$runtime_same_package), 1L)
  expect_match(npRmpi_namespace_hygiene_format(bad.same.scan$runtime_same_package), "inst/Rprofile")

  bad.external <- mkroot()
  dir.create(file.path(bad.external, "R"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    "f <- function(x) stats:::na.omit(x)",
    file.path(bad.external, "R", "bad.R")
  )
  bad.external.scan <- npRmpi_namespace_hygiene_scan(bad.external)
  expect_identical(nrow(bad.external.scan$runtime_external_triple), 1L)
})

test_that("runtime package surfaces avoid forbidden namespace qualification", {
  scan <- npRmpi_namespace_hygiene_scan()

  expect_identical(
    nrow(scan$runtime_same_package),
    0L,
    info = npRmpi_namespace_hygiene_format(scan$runtime_same_package)
  )
  expect_identical(
    nrow(scan$runtime_external_triple),
    0L,
    info = npRmpi_namespace_hygiene_format(scan$runtime_external_triple)
  )
})

test_that("profile startup template resolves the worker loop through a namespace binding", {
  profile <- file.path(npRmpi_namespace_hygiene_root(), "inst", "Rprofile")
  tracker <- new.env(parent = emptyenv())
  fake.ns <- new.env(parent = baseenv())

  assign(".npRmpi_worker_loop", function(...) {
    tracker$worker_called <- TRUE
    tracker$args <- list(...)
    invisible(NULL)
  }, envir = fake.ns)

  env <- new.env(parent = baseenv())
  env$library <- function(package, ..., logical.return = FALSE) TRUE
  env$asNamespace <- function(ns) {
    tracker$namespace_name <- ns
    fake.ns
  }
  env$mpi.comm.size <- function(comm = 0) if (identical(comm, 0)) 2L else 1L
  env$mpi.comm.rank <- function(comm = 0) 1L
  env$mpi.comm.dup <- function(...) invisible(NULL)
  env$mpi.barrier <- function(...) invisible(NULL)
  env$np.mpi.initialize <- function(...) {
    tracker$init_called <- TRUE
    invisible(TRUE)
  }
  env$mpi.comm.free <- function(...) {
    tracker$free_called <- TRUE
    invisible(NULL)
  }
  env$mpi.quit <- function(...) {
    tracker$quit_called <- TRUE
    invisible(NULL)
  }
  env$setHook <- function(...) invisible(NULL)
  env$packageEvent <- function(...) NULL
  env$search <- function() c(".GlobalEnv", "package:stats")
  env$suppressPackageStartupMessages <- function(expr) expr
  env$slave.hostinfo <- function(...) invisible(NULL)
  env$quit <- function(...) stop("quit should not be called in the profile source probe", call. = FALSE)

  sys.source(profile, envir = env)

  expect_identical(tracker$namespace_name, "npRmpi")
  expect_true(isTRUE(tracker$init_called))
  expect_true(isTRUE(tracker$worker_called))
  expect_true(isTRUE(tracker$free_called))
  expect_true(isTRUE(tracker$quit_called))
  expect_identical(tracker$args$comm, 1)
  expect_identical(tracker$args$loop.label, "profile slave")
})
