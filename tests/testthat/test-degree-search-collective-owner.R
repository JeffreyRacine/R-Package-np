test_that("automatic degree search distinguishes command ownership from pool flags", {
  f <- .npRmpi_degree_search_collective
  e <- new.env(parent = environment(f))
  environment(f) <- e
  grid <- expand.grid(local = c(FALSE, TRUE), automatic = c(FALSE, TRUE),
                      manual = c(FALSE, TRUE), pool = c(FALSE, TRUE),
                      size = c(1L, 2L, 3L))
  old <- options(npRmpi.local.regression.mode = FALSE)
  on.exit(options(old), add = TRUE)
  for (i in seq_len(nrow(grid))) {
    z <- grid[i, ]
    e$queries <- 0L
    e$pool.queries <- 0L
    e$.npRmpi_autodispatch_in_context <- function() z$automatic
    e$.npRmpi_autodispatch_called_from_bcast <- function() z$manual
    e$.npRmpi_has_active_slave_pool <- function(comm) {
      expect_identical(comm, 7L)
      e$pool.queries <- e$pool.queries + 1L
      z$pool
    }
    e$mpi.comm.size <- function(comm) {
      expect_identical(comm, 7L)
      e$queries <- e$queries + 1L
      z$size
    }
    options(npRmpi.local.regression.mode = z$local)
    command <- z$automatic || z$manual
    expect_identical(f(7L), !z$local && if (command) z$size >= 2L else z$pool)
    expect_identical(e$queries, as.integer(!z$local && command))
    expect_identical(e$pool.queries, as.integer(!z$local && !command))
  }
})

test_that("a failed communicator query is not hidden by local fallback", {
  f <- .npRmpi_degree_search_collective
  e <- new.env(parent = environment(f)); environment(f) <- e
  e$.npRmpi_autodispatch_in_context <- function() TRUE
  e$mpi.comm.size <- function(...) stop("test communicator unavailable")
  old <- options(npRmpi.local.regression.mode = FALSE)
  on.exit(options(old), add = TRUE)
  expect_error(f(), "test communicator unavailable")
  options(npRmpi.local.regression.mode = TRUE)
  expect_identical(f(), FALSE)
})
