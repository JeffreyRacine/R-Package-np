test_that("MPI collectors retain NULL values as indexed slots", {
  ns <- asNamespace("npRmpi")
  fixtures <- list(
    all_null = rep(list(NULL), 3L),
    first_null = list(NULL, list(id = 2L), list(id = 3L)),
    middle_null = list(list(id = 1L), NULL, list(id = 3L)),
    last_null = list(list(id = 1L), list(id = 2L), NULL),
    nested = list(list(NULL), simpleError("condition as value"), list(a = 3L)),
    ordinary = list(list(id = 1L), list(id = 2L), list(id = 3L)))
  for (payload in fixtures) for (order in list(1:3, 3:1, c(2L, 1L, 3L))) {
    for (owner in c("remote", "scatter", "dynamic")) {
      e <- new.env(parent = ns)
      name <- if (owner == "remote") "mpi.remote.exec" else ".npRmpi_fanout_apply"
      f <- get(name, ns); environment(f) <- e
      e$mpi.comm.size <- function(...) 4L
      e$mpi.bcast.cmd <- e$.npRmpi_bcast_prepared <- function(...) invisible(NULL)
      e$.npRmpi_transport_trace <- function(...) invisible(NULL)
      e$runif <- function(...) 0.5
      cursor <- 0L
      if (owner == "remote") {
        e$mpi.allgather <- function(...) c(0L, 0L, rep(c(-1L, 0L), 3L))
        e$mpi.any.source <- function(...) -1L
        e$mpi.recv.Robj <- function(...) {cursor <<- cursor + 1L; payload[[order[[cursor]]]]}
        e$mpi.get.sourcetag <- function(...) c(order[[cursor]], 1L)
        got <- f(NULL, simplify = FALSE)
        expected <- setNames(payload, paste0("slave", 1:3))
      } else {
        tx <- new.env(parent = emptyenv()); tx$rank <- rep("busy", 3L)
        tx$assigned <- vector("list", 3L)
        e$.npRmpi_fanout_new <- function(...) tx
        e$.npRmpi_fanout_header <- function(...) list()
        e$.npRmpi_scatter_prepare <- function(...) list()
        e$.npRmpi_scatter_prepared <- e$.npRmpi_fanout_send <-
          e$.npRmpi_fanout_native <- function(...) invisible(NULL)
        e$.npRmpi_fanout_run <- function(tx, expr) force(expr)
        e$.npRmpi_fanout_receive <- function(...) {
          cursor <<- cursor + 1L
          if (cursor == 3L) tx$rank <- rep("closed", 3L)
          id <- order[[cursor]]
          list(kind = "result", tasks = id, value = list(payload[[id]]))
        }
        got <- f(as.list(1:3), identity, list(), comm = 1L, dynamic = owner == "dynamic")
        expected <- payload
      }
      expect_identical(got, expected, info = paste(owner, paste(order, collapse = ",")))
    }
  }
})
