test_that("completed local failures preserve conditions without classifying values", {
  original <- structure(list(message = "pilot failure", call = quote(variance()), token = 31L),
                        class = c("pilot_test_error", "error", "condition"))
  expect_identical(.npRmpi_capture_local_work(original), original)
  captured <- .npRmpi_capture_local_work(stop(original))
  expect_s3_class(captured, "npRmpi_local_failure")
  testthat::local_mocked_bindings(
    .npRmpi_autodispatch_in_context = function(...) FALSE,
    .package = "npRmpi"
  )
  expect_identical(tryCatch(.npRmpi_raise_completed_failure(captured), error = identity),
                   original)
})

test_that("the root pilot sends its failure before publishing an error", {
  old <- options(npRmpi.autodispatch.context = TRUE)
  on.exit(options(old), add = TRUE)
  sent <- NULL
  sends <- 0L
  original <- simpleError("pilot failed", call = quote(scale.fit()))
  testthat::local_mocked_bindings(
    mpi.comm.size = function(...) 2L,
    mpi.comm.rank = function(...) 0L,
    .npRmpi_with_local_regression = function(expr) force(expr),
    .nplsqreg_scale_pilot_fit = function(...) stop(original),
    mpi.bcast.Robj = function(obj, ...) { sent <<- obj; sends <<- sends + 1L },
    .package = "npRmpi"
  )
  got <- tryCatch(.nplsqreg_scale_pilot(data.frame(x = 1:4), 1:4,
       dots = list(bwtype = "generalized_nn"), regtype.pilot = "ll"),
       error = identity)
  expect_identical(sends, 1L)
  expect_identical(sent[["condition", exact = TRUE]], original)
  expect_s3_class(got, "npRmpi_coordinated_error")
  expect_identical(got[["cause", exact = TRUE]], original)
})
