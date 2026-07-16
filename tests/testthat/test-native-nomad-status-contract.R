test_that("native NOMAD status decoding preserves R condition semantics", {
  decode <- getFromNamespace(".np_nomad_native_status", "np")

  ok <- list(status = 0L, result_status = 0L, message = "")
  expect_invisible(decode(ok, "native test route"))

  expect_error(
    decode(
      list(status = 1L, result_status = 0L, message = "ordinary failure"),
      "native test route"
    ),
    "native test route failed.*ordinary failure"
  )
  expect_error(
    decode(list(status = 0L, result_status = NA_integer_), "native test route"),
    "native test route failed"
  )

  interrupted <- tryCatch(
    decode(
      list(status = 4L, result_status = 4L, message = "user interrupt"),
      "native test route"
    ),
    interrupt = identity
  )
  expect_s3_class(interrupted, "interrupt")
  expect_identical(conditionMessage(interrupted), "user interrupt")
})
