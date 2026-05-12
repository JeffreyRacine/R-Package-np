test_that("mpi.cart.get parser validates output shape", {
  parse.fun <- getFromNamespace(".npRmpi_parse_cart_get", "npRmpi")

  parsed <- parse.fun(1:6, 2L)
  expect_identical(parsed$dims, 1:2)
  expect_identical(parsed$periods, 3:4)
  expect_identical(parsed$coords, 5:6)

  expect_error(parse.fun(1:5, 2L), "expected 6 values")
  expect_error(parse.fun(1:6, 0L), "single positive integer")
})

test_that("tail slave index helper stays in bounds for small sizes", {
  idx.fun <- getFromNamespace(".npRmpi_tail_slave_indices", "npRmpi")

  expect_identical(idx.fun(2L), 1L)
  expect_identical(idx.fun(3L), 1:2)
  expect_identical(idx.fun(9L), 7:8)
  expect_error(idx.fun(1L), "integer >= 2")
})

test_that("npRmpi.quit is harmless before MPI initialization", {
  old <- getOption("npRmpi.mpi.initialized")
  old.pool <- getOption("npRmpi.pool.active")
  on.exit(options(npRmpi.mpi.initialized = old,
                  npRmpi.pool.active = old.pool), add = TRUE)

  options(npRmpi.mpi.initialized = FALSE)
  expect_false(npRmpi.quit())
})

test_that("active slave pool guard fails closed before npRmpi.init", {
  has.pool <- getFromNamespace(".npRmpi_has_active_slave_pool", "npRmpi")
  old.pool <- getOption("npRmpi.pool.active")
  old.ctx <- getOption("npRmpi.autodispatch.context")
  on.exit(options(npRmpi.pool.active = old.pool), add = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)

  options(npRmpi.pool.active = FALSE)
  expect_false(has.pool())
  expect_error(
    getFromNamespace(".npRmpi_require_active_slave_pool", "npRmpi")(
      where = "abuse probe"
    ),
    "requires an active MPI slave pool"
  )

  options(npRmpi.autodispatch.context = TRUE)
  expect_silent(
    getFromNamespace(".npRmpi_require_active_slave_pool", "npRmpi")(
      where = "worker payload"
    )
  )
})
