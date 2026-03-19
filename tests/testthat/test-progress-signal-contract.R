test_that("native progress signal accepts scalar payloads safely", {
  library(npRmpi)
  expect_silent(.Call(
    "C_np_progress_signal",
    "bandwidth_activity_step",
    "bandwidth",
    1L,
    2L,
    PACKAGE = "npRmpi"
  ))
})
