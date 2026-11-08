test_that("loaded worker identity binds its native source", {
  r <- .npRmpi_worker_identity_record()
  expect_true(isTRUE(r$ok))
  expect_identical(r$source$native, .Call(C_np_worker_native_identity))
  expect_identical(.npRmpi_worker_identity_problems(list(r, r)), c("", ""))
  alternative <- r
  alternative$library <- "/another/installation"
  alternative$platform <- "another-supported-platform"
  expect_identical(.npRmpi_worker_identity_problems(list(r, alternative)), c("", ""))
})

test_that("worker identity refuses mismatches and malformed records", {
  r <- .npRmpi_worker_identity_record()
  changed <- r
  changed$source$r <- paste(rep("0", 32), collapse = "")
  expect_true(nzchar(.npRmpi_worker_identity_problems(list(r, changed))[[2L]]))
  changed <- r
  changed$native <- paste(rep("0", 32), collapse = "")
  expect_true(nzchar(.npRmpi_worker_identity_problems(list(r, changed))[[2L]]))
  for (bad in list(NULL, 1L, list(), list(ok = FALSE), list(ok = TRUE)))
    expect_true(nzchar(.npRmpi_worker_identity_problems(list(r, bad))[[2L]]))
})
