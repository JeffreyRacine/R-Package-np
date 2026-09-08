test_that("index row failures use the existing gather and retain their cause", {
  cause <- structure(list(message = "local index fit failed",
                          call = quote(npreg(bws = b)), stage = "fit", id = 17L),
                     class = c("index_test_error", "error", "condition"))
  state <- new.env(parent = emptyenv())
  state$gathers <- 0L
  testthat::local_mocked_bindings(
    .npindex_spmd_row_task = function(...) list(active = TRUE, rows = 1:2),
    .npindex_spmd_active = function(...) TRUE,
    .npRmpi_autodispatch_in_context = function(...) FALSE,
    mpi.comm.rank = function(...) 0L,
    mpi.allgather.Robj = function(payload, ...) {
      state$gathers <- state$gathers + 1L
      payload
    },
    .package = "npRmpi"
  )
  got <- tryCatch(.npindex_spmd_eval_rows(2L, 1L, function(rows) stop(cause),
                                        "index test"), error = identity)
  expect_identical(got, cause)
  expect_identical(state$gathers, 1L)
  good <- .npindex_spmd_eval_rows(2L, 1L, function(rows) matrix(rows, ncol = 1L),
                                 "index test")
  expect_equal(good, matrix(1:2, ncol = 1L))
  expect_identical(state$gathers, 2L)
  expect_error(.npindex_spmd_eval_rows(2L, 1L,
    function(rows) matrix(1:3, ncol = 1L), "index test"), "malformed chunk")
  expect_identical(state$gathers, 3L)
})

test_that("index row completion does not classify successful condition values as exceptions", {
  value <- structure(matrix(1:2, ncol = 1L), class = c("error", "matrix", "array"))
  testthat::local_mocked_bindings(
    .npindex_spmd_row_task = function(...) list(active = TRUE, rows = 1:2),
    .npindex_spmd_active = function(...) TRUE,
    .npRmpi_autodispatch_in_context = function(...) FALSE,
    mpi.comm.rank = function(...) 0L,
    mpi.allgather.Robj = function(payload, ...) payload,
    .package = "npRmpi"
  )
  got <- .npindex_spmd_eval_rows(2L, 1L, function(rows) value, "index test")
  expect_equal(got, matrix(1:2, ncol = 1L))
})
