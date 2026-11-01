test_that("object all-gather preserves rank payloads without forcing full collection", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  withr::defer(close_mpi_slaves(force = TRUE))
  ns <- asNamespace("npRmpi")
  command <- get(".npRmpi_bcast_cmd_expr", ns)
  # Use the real wrapper and native collectives on every rank. Only explicit
  # gc() is intercepted; R's ordinary allocation-triggered GC remains enabled.
  result <- command(quote(local({
    fn <- get("mpi.allgather.Robj", asNamespace("npRmpi"))
    env <- new.env(parent = environment(fn)); environment(fn) <- env
    env$gc <- function(...) stop("object all-gather forced a full GC")
    rank <- mpi.comm.rank(1L); size <- mpi.comm.size(1L)
    make <- function(i) list(rank = i, empty = raw(0),
      bytes = as.raw(rep(0:255, length.out = 257L + i)),
      category = factor(c("b", "a"), levels = c("b", "a", "unused")),
      nested = list(NULL, seq_len(i + 1L)))
    expected <- lapply(0:(size - 1L), make)
    for (iteration in 1:3) stopifnot(identical(fn(list(make(rank))), expected))
    large <- as.raw(rep(0:255, length.out = 65536L + rank))
    gathered <- fn(list(large))
    stopifnot(identical(lengths(gathered), 65536L + 0:(size - 1L)))
    stopifnot(all(vapply(seq_along(gathered), function(j)
      identical(gathered[[j]], as.raw(rep(0:255, length.out = 65535L + j))), TRUE)))
    list(rank = rank, size = size, complete = TRUE)
  })), caller.execute = TRUE)
  expect_true(result$complete)
  expect_identical(result$rank, 0L)
  expect_gte(result$size, 2L)
  split <- get(".npRmpi_split_raw_by_counts", ns)
  expect_error(split(raw(2), c(1L, 2L), "test"), "payload length")
  alive <- command(quote(npRmpi:::mpi.allgather.Robj(npRmpi:::mpi.comm.rank(1L))),
                   caller.execute = TRUE)
  expect_identical(alive, 0:(result$size - 1L))
})
