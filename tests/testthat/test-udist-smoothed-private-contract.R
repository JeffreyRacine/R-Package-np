test_that("private CDF percentile support follows the complete corrected row", {
  bound <- .np_udist_smoothed_support()
  z <- c(-8, -sqrt(3), -.4, 0, .7, sqrt(3), 8)
  row <- .np_udist_bc_continuous(data.frame(x = 0), data.frame(x = z), 1)$corrected[, 1]
  expect_equal(row, pnorm(z) + z * dnorm(z) / 2, tolerance = 1e-14)
  expect_true(all(bound[1] <= row & row <= bound[2]))
  expect_equal(row[c(2, 6)], unname(bound), tolerance = 1e-14)
  expect_lt(bound[1], 0)
  expect_gt(bound[2], 1)
  expect_error(.np_udist_smoothed_stream_check(list(1:7)), "CMRG")
})
