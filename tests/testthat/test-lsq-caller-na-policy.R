test_that("LSQ formula NA policies belong to their caller", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 24L
  data <- data.frame(x = seq(.05, .95, length.out = n))
  data$y <- sin(3 * data$x)
  data$y[3L] <- NA_real_
  for (route in c("bandwidth", "fit")) {
    calls <- new.env(parent = emptyenv())
    calls$local <- calls$foreign <- calls$data <- 0L
    foreign <- new.env(parent = baseenv())
    foreign$policy <- function(x) {
      calls$foreign <- calls$foreign + 1L
      stats::na.omit(x)
    }
    formula <- as.formula("y ~ x", env = foreign)
    invoke <- function() {
      policy <- function(x) {
        calls$local <- calls$local + 1L
        stats::na.exclude(x)
      }
      if (route == "bandwidth")
        nplsqregbw(bws = formula, data = {calls$data <- calls$data + 1L; data},
          na.action = policy, bw = .2, scale = rep(1, n), delta = .5,
          regtype = "ll", bandwidth.compute = FALSE)
      else
        nplsqreg(bws = formula, data = {calls$data <- calls$data + 1L; data},
          na.action = policy, bw = .2, scale = rep(1, n), delta = .5,
          regtype = "ll", bandwidth.compute = FALSE)
    }
    result <- invoke()
    expect_identical(calls$local, 1L)
    expect_identical(calls$foreign, 0L)
    expect_identical(calls$data, 1L)
    if (route == "fit") {
      expect_length(fitted(result), n)
      expect_true(is.na(fitted(result)[3L]))
    } else {
      expect_identical(as.integer(result$rows.omit), 3L)
    }
  }
})
