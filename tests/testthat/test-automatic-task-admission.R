admission_fixture <- function() {
  e <- new.env(parent = asNamespace("npRmpi"))
  for (name in c(".npRmpi_automatic_task_owner", ".npRmpi_admit_automatic_task",
                 ".npRmpi_with_local_bootstrap")) {
    f <- get(name, asNamespace("npRmpi"))
    environment(f) <- e
    assign(name, f, e)
  }
  e$.npRmpi_task_origin <- new.env(parent = emptyenv())
  e$.npRmpi_task_origin$automatic <- FALSE
  e$manual <- FALSE; e$pool <- TRUE; e$calls <- 0L; e$fail <- FALSE
  e$.npRmpi_autodispatch_called_from_bcast <- function() e$manual
  e$.npRmpi_has_active_slave_pool <- function(comm) e$pool
  e$.npRmpi_bcast_cmd_expr <- function(expr, comm, caller.execute) {
    stopifnot(!caller.execute)
    e$calls <- e$calls + 1L
    if (e$fail) stop("transport sentinel")
    e$command <- expr
    invisible(NULL)
  }
  e
}

test_that("automatic task admission survives local wrappers and unwinds", {
  withr::local_options(npRmpi.autodispatch = TRUE,
    npRmpi.autodispatch.disable = FALSE, npRmpi.autodispatch.context = FALSE,
    npRmpi.autodispatch.option.sync = "never")
  e <- admission_fixture()
  run <- e$.npRmpi_admit_automatic_task; local <- e$.npRmpi_with_local_bootstrap
  expect_error(run(), "automatic numerical calls require verified option agreement", fixed = TRUE)
  expect_error(local(local(run())), "automatic numerical calls require verified option agreement", fixed = TRUE)
  expect_identical(e$calls, 0L)
  expect_false(e$.npRmpi_task_origin$automatic)
  expect_false(getOption("npRmpi.autodispatch.context"))
  expect_false(getOption("npRmpi.autodispatch.disable"))
  expect_identical(local(local(17L)), 17L) # no task, no refusal
  expect_false(e$.npRmpi_task_origin$automatic)
  for (mode in c("manual", "no-pool", "disabled", "context", "auto-off")) {
    e$manual <- mode == "manual"; e$pool <- mode != "no-pool"
    options(npRmpi.autodispatch = mode != "auto-off",
      npRmpi.autodispatch.disable = mode == "disabled",
      npRmpi.autodispatch.context = mode == "context")
    expect_false(local(local(run())), info = mode)
    expect_false(e$.npRmpi_task_origin$automatic, info = mode)
  }
  options(npRmpi.autodispatch = TRUE, npRmpi.autodispatch.disable = FALSE,
    npRmpi.autodispatch.context = FALSE)
  e$pool <- TRUE; e$manual <- FALSE
  expect_false(local({ e$manual <- TRUE; run() }))
  expect_false(e$.npRmpi_task_origin$automatic)
})

test_that("admission does not draw or alter supported option policies", {
  withr::local_options(npRmpi.autodispatch = TRUE,
    npRmpi.autodispatch.disable = FALSE, npRmpi.autodispatch.context = FALSE,
    npRmpi.autodispatch.option.snapshot = NULL)
  e <- admission_fixture(); set.seed(6001); seed <- .Random.seed
  for (mode in c("onchange", "always")) {
    options(npRmpi.autodispatch.option.sync = mode)
    expect_true(e$.npRmpi_with_local_bootstrap(e$.npRmpi_admit_automatic_task()))
    expect_identical(getOption("npRmpi.autodispatch.option.sync"), mode)
    expect_null(getOption("npRmpi.autodispatch.option.snapshot"))
    expect_identical(.Random.seed, seed)
    expect_identical(e$calls, 0L)
  }
})
