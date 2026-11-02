test_that("command errors retain identity and deferred interrupt evidence", {
  scope <- getFromNamespace(".npRmpi_with_command_interrupt_scope", "npRmpi")
  env <- new.env(parent = environment(scope))
  environment(scope) <- env
  events <- integer()
  env$.npRmpi_mpi_interrupt_scope <- function(action) {
    events <<- c(events, action)
    TRUE
  }
  cause <- structure(list(message = "failure", call = quote(owner(x)), code = 17L),
                     class = c("test_error", "error", "condition"), marker = "kept")
  got <- tryCatch(scope(stop(cause)), error = identity, interrupt = identity)
  expect_s3_class(attr(got, ".npRmpi.secondary.interrupt", exact = TRUE), "interrupt")
  attr(got, ".npRmpi.secondary.interrupt") <- NULL
  expect_identical(got, cause)
  expect_identical(events, c(1L, 2L))

  interrupted <- structure(list(message = "specific interrupt", call = quote(owner())),
                           class = c("interrupt", "condition"))
  expect_identical(tryCatch(scope(stop(interrupted)), interrupt = identity), interrupted)
  expect_s3_class(tryCatch(scope(1L), interrupt = identity), "interrupt")
})

test_that("nested command scopes preserve one primary error and close both levels", {
  scope <- getFromNamespace(".npRmpi_with_command_interrupt_scope", "npRmpi")
  env <- new.env(parent = environment(scope))
  environment(scope) <- env
  depth <- 0L
  events <- integer()
  env$.npRmpi_mpi_interrupt_scope <- function(action) {
    events <<- c(events, action)
    depth <<- depth + if (action == 1L) 1L else -1L
    if (action == 1L) TRUE else depth == 0L
  }
  cause <- simpleError("nested failure", call = quote(nested()))
  got <- tryCatch(scope(scope(stop(cause))), error = identity, interrupt = identity)
  expect_s3_class(attr(got, ".npRmpi.secondary.interrupt", exact = TRUE), "interrupt")
  attr(got, ".npRmpi.secondary.interrupt") <- NULL
  expect_identical(got, cause)
  expect_identical(events, c(1L, 1L, 2L, 2L))
  expect_identical(depth, 0L)
})
