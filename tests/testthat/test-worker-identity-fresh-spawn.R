test_that("fresh spawn admits worker identity after merge and before public use", {
  e <- new.env(parent=environment(npRmpi.init)); events <- character()
  record <- function(x) events <<- c(events,x)
  e$is.loaded <- function(...) TRUE
  e$.Platform <- .Platform; e$.Platform$OS <- "unix"
  e$.npRmpi_lease_require_open <- function(...) NULL
  e$.npRmpi_session_runtime_state <- function(...) list(initialized=TRUE,finalized=FALSE)
  e$mpi.comm.size <- function(...) 0L
  e$mpi.comm.spawn <- function(...) {record("spawn");1L}
  e$mpi.intercomm.merge <- function(...) {record("merge");TRUE}
  e$mpi.comm.set.errhandler <- function(...) record("handler")
  e$mpi.comm.disconnect <- function(...) record("disconnect")
  e$mpi.bcast <- function(...) record("bcast")
  e$.npRmpi_sync_runtime_options_to_slaves <- function(comm,force) {
    expect_true(force);record("sync")
  }
  fail <- FALSE
  e$.npRmpi_worker_identity_admit <- function(comm) {
    record("identity");if(fail) stop("injected worker mismatch")
  }
  e$mpi.hostinfo <- function(...) record("hostinfo")
  fn <- getFromNamespace("mpi.spawn.Rslaves","npRmpi"); environment(fn)<-e
  fn(nslaves=1L,quiet=FALSE)
  expect_identical(events,c("spawn","merge","handler","disconnect","bcast","bcast","sync","identity","hostinfo"))
  events<-character();fail<-TRUE
  expect_error(fn(nslaves=1L,quiet=FALSE),"injected worker mismatch")
  expect_identical(tail(events,1L),"identity")
  expect_false("hostinfo" %in% events)
})
