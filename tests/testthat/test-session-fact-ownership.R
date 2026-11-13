test_that("saved options cannot rewind package-owned lifecycle facts", {
  ns<-asNamespace("npRmpi")
  state<-get(".npRmpi_spmd_state",ns);oldfacts<-state$lifecycle
  oldopts<-options();on.exit({state$lifecycle<-oldfacts;options(oldopts)},add=TRUE)
  getfact<-get(".npRmpi_session_fact",ns);setfacts<-get(".npRmpi_session_facts",ns)
  before<-options()
  setfacts(npRmpi.mpi.initialized=TRUE,npRmpi.pool.active=TRUE,
    npRmpi.profile.active=FALSE,npRmpi.attach.close.state="open",
    npRmpi.attach.session.id=17L,npRmpi.attach.session.counter=17L)
  established<-state$lifecycle
  options(before);options(npRmpi.autodispatch=FALSE)
  expect_identical(state$lifecycle,established)
  expect_identical(getfact("npRmpi.attach.session.id"),17L)
  expect_identical(getfact("npRmpi.attach.close.state"),"open")
  expect_false(getOption("npRmpi.autodispatch"))
  expect_true(get(".npRmpi_master_local_entry_needed",ns)())
  get(".npRmpi_session_comm_replaced",ns)(2L)
  expect_identical(state$lifecycle,established)
  get(".npRmpi_session_comm_replaced",ns)(1L)
  expect_false(getfact("npRmpi.pool.active"))
  expect_identical(getfact("npRmpi.attach.close.state"),"closed")
  expect_identical(getfact("npRmpi.attach.session.id"),NA_integer_)
  expect_identical(getfact("npRmpi.attach.session.counter"),17L)
  options(established)
  expect_false(getfact("npRmpi.pool.active"))
})

test_that("lifecycle transitions mirror observations but retain preferences", {
  ns<-asNamespace("npRmpi");state<-get(".npRmpi_spmd_state",ns)
  oldfacts<-state$lifecycle;oldopts<-options()
  on.exit({state$lifecycle<-oldfacts;options(oldopts)},add=TRUE)
  setfacts<-get(".npRmpi_session_facts",ns);getfact<-get(".npRmpi_session_fact",ns)
  options(npRmpi.autodispatch=FALSE)
  setfacts(npRmpi.mpi.initialized=TRUE,npRmpi.pool.active=TRUE)
  expect_true(getOption("npRmpi.pool.active"))
  expect_false(getOption("npRmpi.autodispatch"))
  get(".npRmpi_sync_mpi_inactive",ns)()
  expect_false(getfact("npRmpi.mpi.initialized"))
  expect_false(getfact("npRmpi.pool.active"))
  options(npRmpi.mpi.initialized=TRUE,npRmpi.pool.active=TRUE)
  expect_false(getfact("npRmpi.mpi.initialized"))
  expect_false(getfact("npRmpi.pool.active"))
})


test_that("successful communicator mutations invalidate only their replaced loop slot", {
  ns <- asNamespace("npRmpi")
  state <- get(".npRmpi_spmd_state", ns)
  oldfacts <- state$lifecycle
  oldopts <- options()
  on.exit({state$lifecycle <- oldfacts; options(oldopts)}, add = TRUE)
  setfacts <- get(".npRmpi_session_facts", ns)
  getfact <- get(".npRmpi_session_fact", ns)
  e <- new.env(parent = ns)
  native.result <- 1L
  e$.Call <- function(...) native.result
  e$is.loaded <- function(...) TRUE
  wrappers <- c("mpi.comm.dup", "mpi.comm.split", "mpi.comm.free",
                "mpi.comm.disconnect", "mpi.comm.spawn",
                "mpi.comm.get.parent", "mpi.intercomm.merge")
  for (name in wrappers) {
    fun <- get(name, ns)
    environment(fun) <- e
    assign(name, fun, e)
  }
  e$mpi.comm.size <- function(...) 2L
  calls <- function(slot) list(
    list(comm = 0L, newcomm = slot),
    list(comm = 0L, key = 0L, newcomm = slot),
    list(comm = slot), list(comm = slot),
    list(slave = "unused-test-only", nslaves = 2L, intercomm = slot),
    list(comm = slot), list(intercomm = 2L, comm = slot))
  for (slot in c(1L, 2L)) for (i in seq_along(wrappers)) {
    for (success in c(FALSE, TRUE)) {
      setfacts(npRmpi.pool.active = TRUE, npRmpi.attach.close.state = "open",
               npRmpi.attach.session.id = 31L)
      native.result <- if (!success) 0L else if (i == 5L) 2L else 1L
      out <- do.call(get(wrappers[i], e), calls(slot)[[i]])
      expect_identical(out, native.result)
      invalidated <- slot == 1L && success
      expect_identical(getfact("npRmpi.pool.active"), !invalidated)
      expect_identical(getfact("npRmpi.attach.close.state"),
                       if (invalidated) "closed" else "open")
      expect_identical(getfact("npRmpi.attach.session.id"),
                       if (invalidated) NA_integer_ else 31L)
    }
  }
})
