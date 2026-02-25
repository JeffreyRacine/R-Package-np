test_that("npRmpi profiling helper defaults to basic mode and records timing", {
  begin.fun <- getFromNamespace(".npRmpi_profile_bootstrap_begin", "npRmpi")
  end.fun <- getFromNamespace(".npRmpi_profile_bootstrap_end", "npRmpi")
  clear.fun <- getFromNamespace(".npRmpi_profile_clear", "npRmpi")
  last.fun <- getFromNamespace(".npRmpi_profile_last", "npRmpi")

  old.enable <- getOption("npRmpi.profile.enable")
  old.level <- getOption("npRmpi.profile.level")
  on.exit(options(npRmpi.profile.enable = old.enable), add = TRUE)
  on.exit(options(npRmpi.profile.level = old.level), add = TRUE)
  options(npRmpi.profile.enable = TRUE, npRmpi.profile.level = "basic")
  clear.fun()

  ctx <- begin.fun(
    where = "compute.bootstrap.errors.rbandwidth",
    method = "inid",
    B = 9L,
    ntrain = 100L,
    neval = 25L
  )
  expect_true(is.list(ctx))

  out <- end.fun(list(foo = 1L), ctx)
  rec <- last.fun()

  expect_true(is.list(out))
  expect_true(is.list(rec))
  expect_identical(rec$level, "basic")
  expect_identical(rec$where, "compute.bootstrap.errors.rbandwidth")
  expect_true(is.finite(rec$wall_elapsed_sec))
  expect_true(rec$wall_elapsed_sec >= 0)
  expect_true("timing.profile" %in% names(out))
})

test_that("npRmpi profiling helper supports detailed mode", {
  begin.fun <- getFromNamespace(".npRmpi_profile_bootstrap_begin", "npRmpi")
  add.fun <- getFromNamespace(".npRmpi_profile_add_comm_elapsed", "npRmpi")
  end.fun <- getFromNamespace(".npRmpi_profile_bootstrap_end", "npRmpi")
  clear.fun <- getFromNamespace(".npRmpi_profile_clear", "npRmpi")
  last.fun <- getFromNamespace(".npRmpi_profile_last", "npRmpi")

  old.enable <- getOption("npRmpi.profile.enable")
  old.level <- getOption("npRmpi.profile.level")
  on.exit(options(npRmpi.profile.enable = old.enable), add = TRUE)
  on.exit(options(npRmpi.profile.level = old.level), add = TRUE)
  options(npRmpi.profile.enable = TRUE, npRmpi.profile.level = "detailed")
  clear.fun()

  ctx <- begin.fun("compute.bootstrap.errors.sibandwidth", "wild", 5L, 80L, 80L)
  add.fun(0.002, "mpi.applyLB:wild")
  out <- end.fun(list(bar = 2L), ctx)
  rec <- last.fun()

  expect_true(is.list(out))
  expect_true(is.list(rec))
  expect_identical(rec$level, "detailed")
  expect_true(is.finite(rec$comm_elapsed_sec))
  expect_true(rec$comm_elapsed_sec >= 0.002)
  expect_true(is.finite(rec$compute_elapsed_sec))
  expect_true(is.finite(rec$comm_ratio) || is.na(rec$comm_ratio))
  expect_true(is.numeric(rec$comm_calls) || is.integer(rec$comm_calls))
  expect_true(rec$comm_calls >= 1)
  expect_true(is.character(rec$comm_notes))
  expect_true(is.finite(rec$user_sec))
  expect_true(is.finite(rec$system_sec))
  expect_true(inherits(rec$timestamp_start, "POSIXct"))
  expect_true(inherits(rec$timestamp_end, "POSIXct"))
})

test_that("npRmpi profiling can be fully disabled via option", {
  begin.fun <- getFromNamespace(".npRmpi_profile_bootstrap_begin", "npRmpi")
  end.fun <- getFromNamespace(".npRmpi_profile_bootstrap_end", "npRmpi")
  clear.fun <- getFromNamespace(".npRmpi_profile_clear", "npRmpi")
  last.fun <- getFromNamespace(".npRmpi_profile_last", "npRmpi")

  old.enable <- getOption("npRmpi.profile.enable")
  old.level <- getOption("npRmpi.profile.level")
  on.exit(options(npRmpi.profile.enable = old.enable), add = TRUE)
  on.exit(options(npRmpi.profile.level = old.level), add = TRUE)
  options(npRmpi.profile.enable = FALSE, npRmpi.profile.level = "basic")
  clear.fun()

  ctx <- begin.fun("compute.bootstrap.errors.bandwidth", "inid", 3L, 25L, 17L)
  expect_null(ctx)

  out <- end.fun(list(baz = 3L), ctx)
  rec <- last.fun()
  expect_true(is.list(out))
  expect_null(rec)
  expect_false("timing.profile" %in% names(out))
})

test_that("summary timing formatter includes npRmpi bootstrap profile", {
  begin.fun <- getFromNamespace(".npRmpi_profile_bootstrap_begin", "npRmpi")
  add.fun <- getFromNamespace(".npRmpi_profile_add_comm_elapsed", "npRmpi")
  end.fun <- getFromNamespace(".npRmpi_profile_bootstrap_end", "npRmpi")
  clear.fun <- getFromNamespace(".npRmpi_profile_clear", "npRmpi")
  timing.fun <- getFromNamespace("genTimingStr", "npRmpi")

  old.enable <- getOption("npRmpi.profile.enable")
  old.level <- getOption("npRmpi.profile.level")
  old.summary <- getOption("npRmpi.profile.summary")
  on.exit(options(npRmpi.profile.enable = old.enable), add = TRUE)
  on.exit(options(npRmpi.profile.level = old.level), add = TRUE)
  on.exit(options(npRmpi.profile.summary = old.summary), add = TRUE)
  options(
    npRmpi.profile.enable = TRUE,
    npRmpi.profile.level = "detailed",
    npRmpi.profile.summary = TRUE
  )
  clear.fun()

  ctx <- begin.fun("compute.bootstrap.errors.rbandwidth", "inid", 9L, 100L, 25L)
  add.fun(0.001, "mpi.applyLB:inid")
  end.fun(list(), ctx)

  txt <- timing.fun(list(total.time = 1.0, optim.time = 0.2, fit.time = 0.8))
  expect_true(grepl("MPI Bootstrap Profile:", txt, fixed = TRUE))
  expect_true(grepl("comm_ratio=", txt, fixed = TRUE))
})
