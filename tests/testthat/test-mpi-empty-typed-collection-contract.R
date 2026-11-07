test_that("typed MPI collection retains zero-length worker slices", {
  ns <- asNamespace("npRmpi")
  for (type in c("integer", "double", "raw"))
    for (lens in list(c(0L, 0L), c(0L, 2L), c(2L, 0L), c(1L, 3L), c(2L, 2L)))
      for (simplify in c(FALSE, TRUE)) {
        payload <- lapply(seq_along(lens), function(i) switch(type,
          integer = rep.int(as.integer(i), lens[[i]]),
          double = rep.int(as.double(i), lens[[i]]),
          raw = as.raw(rep.int(i, lens[[i]]))))
        e <- new.env(parent = ns)
        f <- get("mpi.remote.exec", ns); environment(f) <- e
        e$mpi.comm.size <- function(...) 3L
        e$runif <- function(...) 0.5
        e$mpi.bcast.cmd <- e$.npRmpi_bcast_prepared <-
          e$.npRmpi_transport_trace <- function(...) invisible(NULL)
        code <- switch(type, integer = 1L, double = 2L, raw = 4L)
        e$mpi.allgather <- function(...) c(0L, 0L, code, lens[[1L]], code, lens[[2L]])
        e$mpi.gather <- e$mpi.gatherv <- function(x, ...) c(x, do.call(c, payload))
        expected <- if (simplify && lens[[1L]] == lens[[2L]]) {
          data.frame(matrix(do.call(c, payload), nrow = lens[[1L]], ncol = 2L))
        } else setNames(payload, paste0("slave", 1:2))
        got <- tryCatch(f(NULL, simplify = simplify), error = identity)
        expect_identical(got, expected,
                         info = paste(type, paste(lens, collapse = ","), simplify))
      }
})
