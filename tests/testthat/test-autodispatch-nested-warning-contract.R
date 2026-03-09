warn_nested <- getFromNamespace(".npRmpi_autodispatch_warn_nested", "npRmpi")

test_that("nested autodispatch warning still fires outside profile/manual context", {
  old.warned <- getOption("npRmpi.autodispatch.warned.nested")
  old.profile <- getOption("npRmpi.profile.active")
  old.manual <- getOption("npRmpi.manual.bcast.context")
  on.exit(options(npRmpi.autodispatch.warned.nested = old.warned), add = TRUE)
  on.exit(options(npRmpi.profile.active = old.profile), add = TRUE)
  on.exit(options(npRmpi.manual.bcast.context = old.manual), add = TRUE)

  options(
    npRmpi.autodispatch.warned.nested = FALSE,
    npRmpi.profile.active = FALSE,
    npRmpi.manual.bcast.context = FALSE
  )

  expect_warning(
    warn_nested(),
    "detected active mpi.bcast.cmd context; skipping nested auto-dispatch for this call"
  )
  expect_true(isTRUE(getOption("npRmpi.autodispatch.warned.nested")))
})

test_that("nested autodispatch warning is suppressed in canonical profile/manual context", {
  old.warned <- getOption("npRmpi.autodispatch.warned.nested")
  old.profile <- getOption("npRmpi.profile.active")
  old.manual <- getOption("npRmpi.manual.bcast.context")
  on.exit(options(npRmpi.autodispatch.warned.nested = old.warned), add = TRUE)
  on.exit(options(npRmpi.profile.active = old.profile), add = TRUE)
  on.exit(options(npRmpi.manual.bcast.context = old.manual), add = TRUE)

  options(
    npRmpi.autodispatch.warned.nested = FALSE,
    npRmpi.profile.active = TRUE,
    npRmpi.manual.bcast.context = TRUE
  )

  expect_warning(warn_nested(), NA)
  expect_false(isTRUE(getOption("npRmpi.autodispatch.warned.nested")))
})
