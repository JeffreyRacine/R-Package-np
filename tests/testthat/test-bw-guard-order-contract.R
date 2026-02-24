test_that("core bw wrappers validate arguments before MPI pool guard", {
  ns <- asNamespace("npRmpi")

  assert_guard_after_validation <- function(fname) {
    fun <- get(fname, envir = ns)
    txt <- paste(deparse(body(fun)), collapse = "\n")
    pos_validate <- regexpr("npValidate", txt, fixed = TRUE)[1]
    pos_guard <- regexpr(".npRmpi_require_active_slave_pool", txt, fixed = TRUE)[1]
    expect_true(pos_validate > 0L)
    expect_true(pos_guard > 0L)
    expect_true(pos_validate < pos_guard)
  }

  wrappers <- c(
    "npregbw.rbandwidth",
    "npcdensbw.conbandwidth",
    "npcdistbw.condbandwidth",
    "npindexbw.sibandwidth",
    "npscoefbw.scbandwidth",
    "npplregbw.default",
    "npplregbw.plbandwidth"
  )

  lapply(wrappers, assert_guard_after_validation)
})
