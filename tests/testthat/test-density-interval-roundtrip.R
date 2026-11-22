test_that("published density tail intervals replay their own endpoints", {
  withr::local_seed(3405)
  withr::local_options(list(np.messages = FALSE, np.plot.progress = FALSE))
  d <- data.frame(x = runif(92))
  b <- npudensbw(~x, data = d, bws = .16, bandwidth.compute = FALSE)
  draw <- function(interval) suppressWarnings(plot(b, errors = "bootstrap", center = "bias-corrected",
    band = "pointwise", interval = interval, B = 59L, random.seed = 33L, neval = 9L, output = "data"))
  a <- draw("separate-tails")
  reported <- a$inference[[1L]]$components$pointwise$interval
  expect_identical(reported, "separate-tails")
  expect_identical(draw(reported)$panels, a$panels)
  expect_false(identical(draw("asymmetric-sd")$panels[[1L]]$intervals, a$panels[[1L]]$intervals))
})
