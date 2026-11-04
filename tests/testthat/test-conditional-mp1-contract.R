test_that("MP1 is a separate explicit construction", {
  expect_identical(np_conditional_bias_control(), np_conditional_bias_control("analytic"))
  expect_identical(np_conditional_bias_control("jackknife-mp1")$construction, "jackknife-mp1")
  expect_error(np_conditional_bias_control("full-t"))
})

test_that("MP1 signed inversion uses opposite endpoint ranks", {
  T <- matrix(c(-4, 1, 0, 2, 3, 5, 1, -2), nrow=2)
  a <- .np_conditional_mp1_reduce(c(7, 9), c(2, 3), T, .1)
  expect_equal(unname(a$bands$pointwise), matrix(c(1, -6, 15, 15), 2))
  expect_equal(unname(a$bands$simultaneous), matrix(c(-3, -6, 15, 21), 2))
  expect_error(.np_conditional_mp1_reduce(c(7,9),c(0,3),T,.1), "invalid MP1")
})
