test_that("namespace hygiene scanner distinguishes defects from allowed patterns", {
  mkroot <- function() {
    root <- tempfile("np-namespace-hygiene-")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    writeLines(c("Package: np", "Version: 0.0.0"), file.path(root, "DESCRIPTION"))
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
  good.scan <- np_namespace_hygiene_scan(good)
  expect_identical(nrow(good.scan$runtime_same_package), 0L)
  expect_identical(nrow(good.scan$runtime_external_triple), 0L)
  expect_identical(nrow(good.scan$runtime_external_double), 1L)
  expect_identical(nrow(good.scan$runtime_parser_literals), 1L)

  bad.same <- mkroot()
  dir.create(file.path(bad.same, "R"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    "np:::.hidden_helper()",
    file.path(bad.same, "R", "bad_same.R")
  )
  bad.same.scan <- np_namespace_hygiene_scan(bad.same)
  expect_identical(nrow(bad.same.scan$runtime_same_package), 1L)
  expect_match(np_namespace_hygiene_format(bad.same.scan$runtime_same_package), "R/bad_same.R")

  bad.external <- mkroot()
  dir.create(file.path(bad.external, "demo"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    "f <- function(x) stats:::na.omit(x)",
    file.path(bad.external, "demo", "bad_external.R")
  )
  bad.external.scan <- np_namespace_hygiene_scan(bad.external)
  expect_identical(nrow(bad.external.scan$runtime_external_triple), 1L)
  expect_match(np_namespace_hygiene_format(bad.external.scan$runtime_external_triple), "demo/bad_external.R")
})

test_that("runtime package surfaces avoid forbidden namespace qualification", {
  scan <- np_namespace_hygiene_scan()

  expect_identical(
    nrow(scan$runtime_same_package),
    0L,
    info = np_namespace_hygiene_format(scan$runtime_same_package)
  )
  expect_identical(
    nrow(scan$runtime_external_triple),
    0L,
    info = np_namespace_hygiene_format(scan$runtime_external_triple)
  )
})
