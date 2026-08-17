test_that("rank-local regression context does not depend on coordinator-only pool state", {
  resolve_context <- getFromNamespace(
    ".npRmpi_rank_local_regression_context", "npRmpi"
  )
  classify <- function(pool = FALSE,
                       autodispatch = FALSE,
                       manual = FALSE,
                       called = FALSE,
                       local = FALSE) {
    old.local <- getOption("npRmpi.local.regression.mode", FALSE)
    options(npRmpi.local.regression.mode = local)
    on.exit(options(npRmpi.local.regression.mode = old.local), add = TRUE)
    testthat::with_mocked_bindings(
      resolve_context(),
      .npRmpi_has_active_slave_pool = function(comm = 1L) pool,
      .npRmpi_autodispatch_in_context = function() autodispatch,
      .npRmpi_manual_bcast_in_context = function() manual,
      .npRmpi_autodispatch_called_from_bcast = function() called,
      .package = "npRmpi"
    )
  }

  expect_false(classify())
  expect_true(classify(pool = TRUE))
  expect_true(classify(autodispatch = TRUE))
  expect_true(classify(manual = TRUE))
  expect_true(classify(called = TRUE))
  expect_true(classify(local = TRUE))
})
