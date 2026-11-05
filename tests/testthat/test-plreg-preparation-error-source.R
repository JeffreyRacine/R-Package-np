test_that("preparation failure exchange retains source and receiver values", {
  payload <- get(".npplregbw_build_plbandwidth_bcast_payload", asNamespace("npRmpi"))
  failure <- structure(list(condition = simpleError("captured preparation")),
                       class = "npRmpi_local_failure")
  for (source in c(0L, 1L, 3L)) for (participant in 0:3) {
    scope <- new.env(parent = environment(payload))
    scope$mpi.comm.rank <- function(...) participant
    scope$mpi.comm.size <- function(...) 4L
    scope$mpi.allreduce <- function(...) source
    scope$exchanges <- 0L
    scope$mpi.bcast.Robj <- function(obj, rank, comm) {
      scope$exchanges <- scope$exchanges + 1L
      if (participant == rank) NULL else failure
    }
    scope$.npplregbw_build_plbandwidth <- function(...) {
      if (participant == source) failure else list(ok = TRUE)
    }
    local.payload <- payload
    environment(local.payload) <- scope
    expect_identical(local.payload(NULL, NULL, NULL, NULL, NULL, NULL), failure)
    expect_identical(scope$exchanges, 1L)
  }
  scope$mpi.allreduce <- function(...) 4L
  scope$.npplregbw_build_plbandwidth <- function(...) list(ok = TRUE)
  scope$exchanges <- 0L
  expect_identical(local.payload(NULL, NULL, NULL, NULL, NULL, NULL), list(ok = TRUE))
  expect_identical(scope$exchanges, 0L)
})

test_that("selected preparation conditions survive actual ranks and pool reuse", {
  skip_on_cran()
  skip_if(isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
  if (!spawn_mpi_slaves()) skip("MPI pool unavailable")
  withr::defer(close_mpi_slaves(force = TRUE))
  withr::local_options(np.messages = FALSE)
  ns <- asNamespace("npRmpi")
  command <- get(".npRmpi_bcast_cmd_expr", ns)
  remote <- get("mpi.remote.exec", ns)
  last <- get("mpi.comm.size", ns)(1L) - 1L
  pids <- remote(Sys.getpid(), simplify = FALSE)
  diagnostics <- get(".npRmpi_lease_diagnostics", ns)
  command(quote({
    options(np.test.preparation.fault = integer(), np.test.preparation.returns = list())
    trace(".npplregbw_build_plbandwidth", where = asNamespace("npRmpi"),
      print = FALSE, tracer = quote({
        rank <- as.integer(mpi.comm.rank(1L))
        if (rank %in% getOption("np.test.preparation.fault"))
          stop(structure(list(message = paste("preparation source", rank),
            call = quote(preparation_fixture()), source_rank = rank, marker = "P3"),
            class = c("np_test_preparation_error", "error", "condition")))
      }))
    trace(".npplregbw_build_plbandwidth_bcast_payload", where = asNamespace("npRmpi"),
      print = FALSE, exit = quote({
        ans <- returnValue()
        if (inherits(ans, "npRmpi_local_failure")) {
          seen <- getOption("np.test.preparation.returns")
          seen[[length(seen) + 1L]] <- ans$condition
          options(np.test.preparation.returns = seen)
        }
      }))
  }), caller.execute = TRUE)
  withr::defer(command(quote({
    untrace(".npplregbw_build_plbandwidth", where = asNamespace("npRmpi"))
    untrace(".npplregbw_build_plbandwidth_bcast_payload", where = asNamespace("npRmpi"))
    options(np.test.preparation.fault = NULL, np.test.preparation.returns = NULL)
  }), caller.execute = TRUE))
  set.seed(703)
  n <- 24L
  z <- data.frame(z = runif(n))
  y <- sin(3 * z$z) + seq_len(n) / n
  f <- factor(rep(c("30", "10", "50"), 8L), levels = c("50", "10", "30", "70"))
  args <- list(xdat = data.frame(x = f), ydat = y, zdat = z,
    regtype = "lp", degree.select = "coordinate", degree.min = 0L,
    degree.max = 1L, search.engine = "cell", nmulti = 1L,
    itmax = 2L, powell.remin = FALSE)
  for (ranks in list(0L, as.integer(last), seq.int(0L, last))) {
    command(bquote(options(np.test.preparation.fault = .(ranks))), caller.execute = TRUE)
    before <- diagnostics()
    err <- tryCatch(do.call(npplregbw, args), error = identity)
    expect_identical(diagnostics(), before)
    cause <- if (inherits(err, "npRmpi_coordinated_error")) err$cause else err
    expect_s3_class(cause, "np_test_preparation_error")
    expect_identical(cause$source_rank, min(ranks))
    expect_identical(cause$marker, "P3")
    expect_identical(conditionCall(cause), quote(preparation_fixture()))
    expect_identical(conditionMessage(cause), paste("preparation source", min(ranks)))
    command(quote(options(np.test.preparation.fault = integer())), caller.execute = TRUE)
    b <- do.call(npplregbw, args)
    expect_true(all(is.finite(fitted(npplreg(b, se = FALSE)))))
    expect_identical(remote(Sys.getpid(), simplify = FALSE), pids)
  }
  seen <- c(list(getOption("np.test.preparation.returns")),
    remote(getOption("np.test.preparation.returns"), simplify = FALSE))
  expect_true(all(vapply(seen, length, integer(1L)) == 3L))
  expect_true(all(vapply(seen, function(x)
    identical(vapply(x, `[[`, integer(1L), "source_rank"), c(0L, as.integer(last), 0L)),
    logical(1L))))
})
