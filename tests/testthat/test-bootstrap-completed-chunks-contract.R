test_that("collective progress credits completed batches, including empty tails", {
  ns <- asNamespace("npRmpi")
  fun <- get(".npRmpi_bootstrap_collective_apply", ns)
  for (rank in 0:3) for (B in c(1L, 9L, 17L, 399L)) {
    env <- new.env(parent = ns)
    env$mpi.comm.size <- function(...) 4L
    env$mpi.comm.rank <- local({ rr <- rank; function(...) rr })
    env$mpi.allreduce <- function(x, ...) x
    environment(fun) <- env
    seen <- new.env(parent = emptyenv())
    seen$completed <- integer()
    seen$progress <- integer()
    out <- fun(B, 8L, function(ids, pos) {
      seen$completed <- c(seen$completed, ids)
      as.double(ids)
    }, progress = function(done) {
      expect_lte(length(seen$completed), done)
      seen$progress <- c(seen$progress, done)
    })
    ids <- seq_len(B)
    expected <- ids[(ids - 1L) %% 4L == rank]
    expect_identical(as.integer(unlist(out)), expected)
    if (rank == 0L) {
      expect_identical(tail(seen$progress, 1L), B)
      expect_true(all(diff(seen$progress) > 0))
      if (B == 399L) expect_lt(seen$progress[[1]], B / 4)
    } else expect_length(seen$progress, 0L)
  }
})

test_that("a failed batch cannot be reported complete", {
  ns <- asNamespace("npRmpi")
  fun <- get(".npRmpi_bootstrap_collective_apply", ns)
  env <- new.env(parent = ns)
  env$mpi.comm.size <- function(...) 4L
  env$mpi.comm.rank <- function(...) 0L
  env$mpi.allreduce <- function(x, ...) x
  environment(fun) <- env
  seen <- new.env(parent = emptyenv())
  seen$n <- 0L
  expect_error(fun(399L, 8L, function(...) stop("sentinel failure"),
    progress = function(done) seen$n <- seen$n + 1L), "sentinel failure")
  expect_identical(seen$n, 0L)
})
