test_that("broadcast sync helpers snapshot global objects via mget", {
  robj.body <- paste(deparse(body(mpi.bcast.Robj2slave), width.cutoff = 500L), collapse = " ")
  rfun.body <- paste(deparse(body(mpi.bcast.Rfun2slave), width.cutoff = 500L), collapse = " ")

  expect_match(robj.body, "mget\\(master\\.objects")
  expect_match(rfun.body, "mget\\(master\\.fun")
  expect_no_match(rfun.body, "lapply\\(lapply\\(")
})
