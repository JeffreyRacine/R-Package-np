nprmpi_internal <- function(name) {
  get(name, envir = asNamespace("npRmpi"), inherits = FALSE)
}

test_that(".np_nomad_require_crs fails fast when crs is missing", {
  require_crs <- nprmpi_internal(".np_nomad_require_crs")
  loaded <- FALSE

  expect_error(
    require_crs(
      version_fn = function(pkg) stop(sprintf("missing %s", pkg), call. = FALSE),
      load_namespace = function(pkg) {
        loaded <<- TRUE
        NULL
      }
    ),
    paste(
      "requires the suggested package 'crs'",
      "\\(>= 0\\.15-41\\) to provide the NOMAD backend;",
      "install\\.packages\\('crs'\\)"
    )
  )

  expect_false(loaded)
})

test_that(".np_nomad_require_crs fails fast when crs is too old", {
  require_crs <- nprmpi_internal(".np_nomad_require_crs")
  loaded <- FALSE

  expect_error(
    require_crs(
      version_fn = function(pkg) package_version("0.15-40"),
      load_namespace = function(pkg) {
        loaded <<- TRUE
        NULL
      }
    ),
    "requires 'crs' \\(>= 0\\.15-41\\); installed version is 0\\.15\\.40"
  )

  expect_false(loaded)
})

test_that(".np_nomad_require_crs succeeds when crs version is sufficient", {
  require_crs <- nprmpi_internal(".np_nomad_require_crs")
  loaded <- FALSE

  observed <- withVisible(
    require_crs(
      version_fn = function(pkg) package_version("0.15-41"),
      load_namespace = function(pkg) {
        loaded <<- identical(pkg, "crs")
        baseenv()
      }
    )
  )

  expect_true(isTRUE(observed$value))
  expect_false(observed$visible)
  expect_true(loaded)
})

test_that(".np_nomad_require_crs reports namespace load failures after version check", {
  require_crs <- nprmpi_internal(".np_nomad_require_crs")

  expect_error(
    require_crs(
      version_fn = function(pkg) package_version("0.15-41"),
      load_namespace = function(pkg) stop("broken namespace", call. = FALSE)
    ),
    "requires 'crs' \\(>= 0\\.15-41\\); failed to load installed version 0\\.15\\.41: broken namespace"
  )
})
