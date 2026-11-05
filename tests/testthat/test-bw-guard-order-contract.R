test_that("core bw wrappers validate arguments before MPI pool guard", {
  ns <- asNamespace("npRmpi")

  assert_guard_after_validation <- function(fname) {
    fun <- get(fname, envir = ns)
    if (identical(fname, "npregbw.rbandwidth")) {
      # Numeric-response forwarding preserves the public call; validation and
      # pool ownership live in its explicit implementation delegate.
      forward <- body(fun)[[2L]]
      expect_identical(forward[[1L]], as.name(".np_regression_forward"))
      expect_identical(forward[[2L]], as.name(".npregbw_bandwidth_impl"))
      fun <- get(".npregbw_bandwidth_impl", envir = ns, inherits = FALSE)
    }
    txt <- paste(deparse(body(fun)), collapse = "\n")
    pos_validate <- regexpr("npValidate", txt, fixed = TRUE)[1]
    pos_guard <- regexpr(".npRmpi_require_active_slave_pool", txt, fixed = TRUE)[1]
    expect_true(pos_validate > 0L, info = fname)
    expect_true(pos_guard > 0L, info = fname)
    expect_true(pos_validate < pos_guard, info = fname)
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


test_that("forwarded regression arguments fail before touching the MPI pool", {
  guard.calls <- 0L
  local_mocked_bindings(.npRmpi_require_active_slave_pool = function(...) {
    guard.calls <<- guard.calls + 1L
    stop("MPI guard reached")
  }, .package = "npRmpi")
  bw.wrapper <- get("npregbw.rbandwidth", asNamespace("npRmpi"), inherits = FALSE)
  expect_error(bw.wrapper(xdat = data.frame(x = 1:3), ydat = 1:3,
    bws = NULL, bandwidth.compute = NA), "bandwidth.compute", fixed = TRUE)
  expect_identical(guard.calls, 0L)
})
