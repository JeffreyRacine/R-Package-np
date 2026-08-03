test_that("exact local operators support multiple right-hand sides", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  set.seed(7321)
  train <- data.frame(x = runif(31), z = runif(31))
  eval <- data.frame(x = seq(0.08, 0.92, length.out = 7),
                     z = seq(0.91, 0.09, length.out = 7))
  rhs <- cbind(a = rnorm(nrow(train)), b = runif(nrow(train)))

  for (kernel in c("gaussian", "beta")) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bw <- npudistbw(
        dat = train,
        bws = if (identical(bwtype, "fixed")) c(0.19, 0.23) else c(5, 6),
        bwtype = bwtype,
        ckertype = kernel,
        ckerbound = if (identical(kernel, "beta")) "range" else "none",
        bandwidth.compute = FALSE
      )

      for (operator.name in c("normal", "integral")) {
        operator <- rep.int(operator.name, ncol(train))
        got <- npRmpi:::.np_exact_operator_apply(
          kbw = bw,
          txdat = train,
          exdat = eval,
          operator = operator,
          rhs = rhs,
          where = "multi-RHS contract test"
        )
        reference <- do.call(cbind, lapply(seq_len(ncol(rhs)), function(j) {
          npRmpi:::.np_exact_operator_apply(
            kbw = bw,
            txdat = train,
            exdat = eval,
            operator = operator,
            rhs = rhs[, j, drop = FALSE],
            where = "single-RHS contract oracle"
          )
        }))

        expect_equal(unname(got), unname(reference), tolerance = 2e-13)
        expect_identical(dim(got), c(nrow(eval), ncol(rhs)))

        if (identical(bwtype, "fixed") || identical(kernel, "beta")) {
          direct <- npRmpi:::.np_direct_operator_apply(
            kbw = bw,
            txdat = train,
            exdat = eval,
            operator = operator,
            rhs = rhs,
            where = "direct multi-RHS contract oracle"
          )
          expect_equal(unname(got), unname(direct), tolerance = 2e-13)
        }
      }
    }
  }
})

test_that("native output extents fail before integer overflow", {
  expect_identical(
    npRmpi:::.np_native_output_extent(7L, 3L, where = "test extent"),
    21L
  )
  expect_error(
    npRmpi:::.np_native_output_extent(.Machine$integer.max, 2L,
                                  where = "test extent"),
    "exceeds the native output-size capacity",
    fixed = TRUE
  )
})
