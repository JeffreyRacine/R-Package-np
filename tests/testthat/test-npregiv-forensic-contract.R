iv_forensic_fixture <- function(n = 26L, seed = 20260721L) {
  set.seed(seed)
  w <- rnorm(n)
  v <- rnorm(n, sd = 0.22)
  z <- 0.45 * w + v
  y <- z^2 - 0.4 * v + rnorm(n, sd = 0.055)
  data.frame(y, z, w)
}

iv_forensic_skip_mpi_numerics <- function() {
  if ("npRmpi" %in% loadedNamespaces())
    skip("installed one- and multi-worker sentinels own npRmpi numerical coverage")
}

test_that("IV scalar controls fail early and deterministically", {
  number <- getFromNamespace(".np_iv_validate_number", "npRmpi")
  integer <- getFromNamespace(".np_iv_validate_integer", "npRmpi")
  flag <- getFromNamespace(".np_iv_validate_flag", "npRmpi")
  expect_error(
    number(c(0.25, 0.5), "constant", 0, 1, TRUE, TRUE),
    "constant must be one finite number", fixed = TRUE
  )
  expect_error(
    integer(2.5, "iterate.max", 2),
    "iterate.max must be one integer", fixed = TRUE
  )
  expect_error(
    flag(c(TRUE, FALSE), "stop.on.increase"),
    "stop.on.increase must be TRUE or FALSE", fixed = TRUE
  )
})

test_that("npregiv starting values and first-state replay are coherent", {
  iv_forensic_skip_mpi_numerics()
})

test_that("multivariate npregiv derivatives are coordinate-specific", {
  orders <- getFromNamespace(".np_iv_partial_orders", "npRmpi")
  expect_identical(orders(2L, 1L), diag(1L, 2L, 2L))
  expect_identical(orders(2L, 2L), diag(2L, 2L, 2L))
  iv_forensic_skip_mpi_numerics()
})

test_that("npregiv owns categorical kernels and bandwidth scaling arguments", {
  iv_forensic_skip_mpi_numerics()
})

test_that("npregivderiv evaluation is output-only and S3 accessors stay training-aligned", {
  iv_forensic_skip_mpi_numerics()
})
