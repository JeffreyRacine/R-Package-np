test_that("conditional bandwidth selectors reject unused named dots", {
  set.seed(20260511)
  d <- data.frame(x = runif(30), y = rnorm(30))

  expect_error(
    npcdensbw(y ~ x, data = d, ckertype = "banana", nmulti = 1L),
    "unused argument"
  )
  expect_error(
    npcdistbw(y ~ x, data = d, ckertype = "banana", nmulti = 1L),
    "unused argument"
  )
})
