# Pure boundary contracts: no MPI pool, launch, or default-lane search matrix.
service_shutdown_fixture <- function() {
  e <- new.env(parent = asNamespace("npRmpi"))
  s <- new.env(parent = emptyenv())
  s$status <- "open"
  s$events <- character()
  s$clock <- 0
  s$seen <- integer()
  s$missing <- integer()
  s$failed <- integer()
  s$wrong <- integer()
  e$.npRmpi_lease_state <- s
  for (name in c(".npscoefbw_nomad_pool_stop", ".npscoefbw_nomad_service_owner",
                 ".npRmpi_lease_require_open", ".npRmpi_lease_poison")) {
    fn <- get(name, asNamespace("npRmpi"))
    environment(fn) <- e
    assign(name, fn, e)
  }
  e$.npRmpi_shutdown_notice <- function(text) invisible(NULL)
  e$.npRmpi_protocol_rank_tag <- function(base, rank, ...) 100L + rank
  e$.npRmpi_attach_close_ack_timeout <- function() .1
  e$Sys.time <- function() { s$clock <- s$clock + .025; s$clock }
  e$Sys.sleep <- function(...) NULL
  e$mpi.send.Robj <- function(obj, dest, ...) {
    s$events <- c(s$events, paste0("send:", dest))
    if (dest %in% s$failed) stop("send failed")
    1L
  }
  e$mpi.iprobe <- function(source, ...) !source %in% c(s$missing, s$seen)
  e$mpi.recv.Robj <- function(source, ...) {
    s$seen <- c(s$seen, source)
    list(type = "stopped", rank = if (source %in% s$wrong) 99L else source)
  }
  list(e = e, s = s, pool = list(nslaves = 3L, comm = 1L, req.base = 10L))
}

test_that("service stop requires valid acknowledgements from every started rank", {
  z <- service_shutdown_fixture()
  expect_null(z$e$.npscoefbw_nomad_pool_stop(z$pool))
  expect_identical(z$s$seen, 1:3)
  expect_identical(z$s$status, "open")
  for (fault in c("failed", "missing", "wrong")) {
    z <- service_shutdown_fixture()
    z$s[[fault]] <- 2L
    ex <- tryCatch({ z$e$.npscoefbw_nomad_pool_stop(z$pool); NULL },
                   error = identity)
    expect_s3_class(ex, "npRmpi_service_stop_error")
    expect_identical(ex$outcome$expected, 1:3)
    expect_true(2L %in% ex$outcome$missing)
    expect_identical(z$s$status, "poisoned")
  }
})

test_that("service owner never retires or retries an uncertain generation", {
  for (fault in c("start", "stop", "retire", "none")) {
    z <- service_shutdown_fixture()
    z$e$.npscoefbw_nomad_pool_start <- function(...) {
      z$s$events <- c(z$s$events, "start")
      if (fault == "start") stop("activation failed")
      z$pool
    }
    z$e$.npscoefbw_nomad_pool_stop <- function(...) {
      z$s$events <- c(z$s$events, "stop")
      if (fault == "stop") stop("stop failed")
      invisible(NULL)
    }
    z$e$.npscoefbw_nomad_context_cleanup <- function(...) {
      z$s$events <- c(z$s$events, "retire")
      if (fault == "retire") stop("retire failed")
      invisible(NULL)
    }
    owner <- z$e$.npscoefbw_nomad_service_owner(list())
    ex <- tryCatch({ owner$start(); owner$cleanup(); NULL }, error = identity)
    before <- z$s$events
    owner$cleanup(primary = ex)
    expect_identical(z$s$events, before)
    if (fault == "none") {
      expect_null(ex)
      expect_identical(before, c("start", "stop", "retire"))
      expect_identical(z$s$status, "open")
    } else {
      expect_s3_class(ex, "error")
      expect_identical(z$s$status, "poisoned")
      expect_identical(sum(before == "retire"), as.integer(fault == "retire"))
    }
  }
})

test_that("failed stop probes and receives are not retried", {
  for (phase in c("probe", "receive")) {
    z <- service_shutdown_fixture()
    count <- 0L
    if (phase == "probe") {
      z$e$mpi.iprobe <- function(source, ...) {
        if (source == 2L) { count <<- count + 1L; stop("probe failed") }
        TRUE
      }
    } else {
      z$e$mpi.recv.Robj <- function(source, ...) {
        if (source == 2L) { count <<- count + 1L; stop("receive failed") }
        list(type = "stopped", rank = source)
      }
    }
    expect_error(z$e$.npscoefbw_nomad_pool_stop(z$pool),
                 class = "npRmpi_service_stop_error")
    expect_identical(count, 1L)
    expect_identical(z$s$status, "poisoned")
  }
})

test_that("shutdown notice cannot invoke a throwing outer message handler", {
  expect_no_error(withCallingHandlers(
    capture.output(.npRmpi_shutdown_notice("shutdown diagnostic"), type = "message"),
    message = function(m) stop("outer message handler")))
})
