test_that("native kernel sums reject incomplete internal options before reads", {
  # No MPI or numeric owner is entered: the option ABI is validated first.
  count <- 28L
  for (entry in c("C_np_kernelsum", "C_np_kernelsum_power12")) {
    for (n in c(0L, count - 3L, count - 2L, count - 1L)) {
      args <- c(list(entry), rep(list(double()), 12L),
                list(integer(n), 1, 0L, 0L, 0L, double(), double(),
                     PACKAGE = "npRmpi"))
      expect_error(do.call(.Call, args),
                   "invalid internal option vector", fixed = TRUE)
    }
  }
})
