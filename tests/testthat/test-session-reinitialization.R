# Pure startup-owner seam: no processes, transport mocks only. Real MPI lifecycle
# and raw-zero close proof belong in the external R36 qualification fixtures.
reinit_fixture <- function(world = 1L, pool = TRUE, owner = "spawn",
                           blocked = NULL, rank = 0L) {
  e <- new.env(parent = environment(npRmpi.init))
  e$events <- character()
  e$asNamespace <- function(package) e
  e$.npRmpi_profile_clear <- function() invisible(NULL)
  record <- function(x) e$events <- c(e$events, x)
  e$.npRmpi_abort_if_rmpi_attached <- function(...) NULL
  e$.npRmpi_fanout_assert_idle <- function(...) {
    if (identical(blocked, "busy")) stop("busy pool")
  }
  e$.npRmpi_lease_require_open <- function(...) {
    if (identical(blocked, "poisoned")) stop("poisoned pool")
  }
  e$mpi.comm.size <- function(comm) if (comm == 0L) world else if (pool) 3L else 0L
  e$.npRmpi_session_has_active_pool <- function(...) pool
  e$.npRmpi_worker_identity_admit <- function(...) {
    record("identity")
    if (identical(blocked, "identity")) stop("identity mismatch")
  }
  e$.npRmpi_session_reset_spmd_state <- function(...) record("local-reset")
  e$.npRmpi_session_reset_spmd_state_allranks <- function(...) record("all-reset")
  e$mpi.spawn.Rslaves <- function(...) record("spawn")
  e$mpi.bcast.cmd <- function(...) record("initialize-all")
  e$mpi.comm.dup <- function(...) record("dup")
  e$.npRmpi_session_ensure_comm <- function(...) 3L
  e$np.mpi.initialize <- function(...) record("initialize")
  e$mpi.barrier <- function(...) record("barrier")
  e$mpi.comm.rank <- function(...) rank
  e$.npRmpi_attach_next_session_id <- function(...) 72L
  e$mpi.bcast <- function(x, ...) { record("session-id"); 72L }
  e$.npRmpi_attach_register_exit_finalizer <- function(...) record("finalizer")
  e$.npRmpi_session_attach_worker_loop <- function(...) { record("worker-loop"); FALSE }
  fn <- npRmpi.init
  environment(fn) <- e
  e$init <- fn
  e
}

reinit_options <- function(owner = "spawn") {
  keys <- c("npRmpi.mpi.initialized", "npRmpi.pool.active", "npRmpi.profile.active",
            "npRmpi.attach.close.state", "npRmpi.attach.session.id",
            "npRmpi.autodispatch", "np.messages", "npRmpi.master.only",
            "npRmpi.autodispatch.verify.options", "npRmpi.autodispatch.option.sync",
            "npRmpi.autodispatch.option.snapshot")
  old <- setNames(lapply(keys, getOption), keys)
  options(npRmpi.mpi.initialized = TRUE, npRmpi.pool.active = owner != "fresh",
          npRmpi.profile.active = identical(owner, "profile"),
          npRmpi.attach.close.state = if (owner == "attach") "open" else "closed",
          npRmpi.attach.session.id = if (owner == "attach") 41L else NA_integer_,
          npRmpi.autodispatch = FALSE, np.messages = TRUE)
  old
}

test_that("reuse preserves actual attach ownership and resets all ranks", {
  old <- reinit_options("attach"); on.exit(options(old), add = TRUE)
  for (mode in c("spawn", "attach", "auto")) {
    e <- reinit_fixture(world = 3L, owner = "attach")
    expect_true(e$init(mode = mode, quiet = TRUE, np.messages = FALSE))
    expect_identical(e$events, c("identity", "all-reset"))
    expect_identical(getOption("npRmpi.attach.close.state"), "open")
    expect_identical(getOption("npRmpi.attach.session.id"), 41L)
    expect_true(getOption("npRmpi.autodispatch"))
    expect_false(getOption("np.messages"))
  }
})

test_that("spawn and profile adoption keep their existing lifecycle", {
  old <- reinit_options(); on.exit(options(old), add = TRUE)
  for (owner in c("spawn", "profile")) {
    reinit_options(owner)
    e <- reinit_fixture(world = if (owner == "spawn") 1L else 3L)
    expect_true(e$init(mode = "spawn", nslaves = 7L, quiet = TRUE))
    expect_identical(e$events, c("identity", "all-reset"))
    expect_identical(getOption("npRmpi.attach.close.state"), "closed")
    expect_identical(getOption("npRmpi.attach.session.id"), NA_integer_)
    expect_identical(getOption("npRmpi.profile.active"), owner == "profile")
  }
})

test_that("invalid transitions and existing guards precede mutation", {
  old <- reinit_options(); on.exit(options(old), add = TRUE)
  cases <- list(list("spawn", 1L, "attach", NULL, "pre-launched"),
                list("profile", 3L, "attach", NULL, "profile-managed"),
                list("profile", 3L, "auto", NULL, "profile-managed"),
                list("attach", 3L, "spawn", "busy", "busy"),
                list("attach", 3L, "spawn", "poisoned", "poisoned"),
                list("attach", 3L, "attach", "identity", "identity mismatch"))
  for (z in cases) {
    reinit_options(z[[1L]])
    before <- options()
    e <- reinit_fixture(world = z[[2L]], blocked = z[[4L]])
    expect_error(e$init(mode = z[[3L]], quiet = TRUE, np.messages = FALSE), z[[5L]])
    expect_identical(options(), before)
    expect_identical(e$events, if (identical(z[[4L]], "identity")) "identity" else character())
  }
})

test_that("fresh attach includes a pre-duplicated communicator and workers", {
  old <- reinit_options("fresh"); on.exit(options(old), add = TRUE)
  for (pool in c(FALSE, TRUE)) for (rank in 0:1) {
    reinit_options("fresh")
    e <- reinit_fixture(world = 3L, pool = pool, rank = rank)
    expect_identical(e$init(mode = "auto", quiet = TRUE), rank == 0L)
    expect_true(all(c("local-reset", "dup", "initialize", "barrier", "session-id") %in% e$events))
    expect_false("all-reset" %in% e$events)
    expect_identical(getOption("npRmpi.attach.session.id"), 72L)
    expect_identical(getOption("npRmpi.attach.close.state"), "open")
    expect_true(if (rank == 0L) "identity" %in% e$events else "worker-loop" %in% e$events)
  }
})

test_that("fresh spawn retains startup and resets stale attach metadata", {
  old <- reinit_options("fresh"); on.exit(options(old), add = TRUE)
  options(npRmpi.attach.session.id = 99L, npRmpi.attach.close.state = "open")
  e <- reinit_fixture(pool = FALSE)
  expect_true(e$init(quiet = TRUE))
  expect_identical(e$events, c("local-reset", "spawn", "initialize-all"))
  expect_identical(getOption("npRmpi.attach.session.id"), NA_integer_)
  expect_identical(getOption("npRmpi.attach.close.state"), "closed")
})


test_that("attach checks the canonical communicator before reuse", {
  old <- reinit_options("attach"); on.exit(options(old), add = TRUE)
  e <- reinit_fixture(world = 3L)
  e$.npRmpi_fanout_assert_idle <- function(comm, ...) {
    if (comm == 1L) stop("actual owner busy")
  }
  before <- options()
  expect_error(e$init(mode = "attach", comm = 2L, quiet = TRUE), "actual owner busy")
  expect_identical(options(), before)
  expect_identical(e$events, character())
})
