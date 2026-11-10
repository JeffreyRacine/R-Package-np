admission_fixture <- function() {
  e <- new.env(parent = asNamespace("npRmpi"))
  for (name in c(".npRmpi_automatic_task_owner", ".npRmpi_admit_automatic_task",
                 ".npRmpi_with_local_bootstrap", ".npRmpi_autodispatch_should_sync_options")) {
    f <- get(name, asNamespace("npRmpi"))
    environment(f) <- e
    assign(name, f, e)
  }
  e$.npRmpi_option_agreement <- new.env(parent=emptyenv())
  e$.npRmpi_option_agreement$snapshot <- NULL
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

test_that("automatic task options reuse the incumbent cache without draws", {
  withr::local_options(npRmpi.autodispatch = TRUE,
    npRmpi.autodispatch.disable = FALSE, npRmpi.autodispatch.context = FALSE,
    npRmpi.autodispatch.option.sync = "onchange",
    np.tree = FALSE)
  e <- admission_fixture(); run <- e$.npRmpi_admit_automatic_task
  set.seed(6001); seed <- .Random.seed
  expect_true(run()); expect_identical(e$calls, 1L)
  expect_true(run()); expect_identical(e$calls, 1L)
  options(np.tree = TRUE)
  expect_true(e$.npRmpi_with_local_bootstrap(run()))
  expect_identical(e$calls, 2L)
  expect_true(e$.npRmpi_option_agreement$snapshot[["np.tree"]])
  options(npRmpi.autodispatch.option.sync = "always")
  expect_true(run()); expect_true(run()); expect_identical(e$calls, 4L)
  expect_identical(.Random.seed, seed)
  e$fail <- TRUE
  expect_error(e$.npRmpi_with_local_bootstrap(run()), "transport sentinel", fixed = TRUE)
  expect_null(e$.npRmpi_option_agreement$snapshot)
  expect_false(e$.npRmpi_task_origin$automatic)
  expect_false(getOption("npRmpi.autodispatch.context"))
  e$fail <- FALSE; options(npRmpi.autodispatch.option.sync = "onchange")
  expect_true(run()); expect_identical(e$calls, 6L)
})


test_that("transported option command cannot resolve worker-global masks", {
  withr::local_options(npRmpi.autodispatch = TRUE,
    npRmpi.autodispatch.disable = FALSE, npRmpi.autodispatch.context = FALSE,
    npRmpi.autodispatch.option.sync = "always",
    npRmpi.autodispatch.verify.options = TRUE,
    np.tree = FALSE)
  e <- admission_fixture()
  expect_true(e$.npRmpi_admit_automatic_task())
  shadow <- new.env(parent = baseenv())
  shadow$options <- shadow$getOption <- function(...) stop("ambient mask")
  options(np.tree = TRUE)
  expect_null(eval(e$command, envir = shadow))
  expect_false(getOption("np.tree"))
  expect_identical(ls(shadow), c("getOption", "options"))
})
