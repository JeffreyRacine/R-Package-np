test_that("conditional mode vectors retain names for subsequent prediction", {
  withr::local_options(np.messages = FALSE)
  set.seed(3899)
  x <- runif(80)
  y <- factor(ifelse(x + rnorm(80, sd = .3) > .5, "yes", "no"))
  E <- data.frame(x = c(.2, .5, .8))
  for (tree in c(FALSE, TRUE)) {
    withr::local_options(np.tree = tree)
    reference <- npconmode(bws = c(.25, .1), txdat = data.frame(x = x),
      tydat = data.frame(y = y), exdat = E, probabilities = TRUE, gradients = TRUE)
    wrap <- function(...) npconmode(...)
    for (fit in list(
      npconmode(bws = c(.25, .1), txdat = x, tydat = y, exdat = E,
                probabilities = TRUE, gradients = TRUE),
      npconmode(c(.25, .1), x, y, exdat = E,
                probabilities = TRUE, gradients = TRUE),
      wrap(bws = c(.25, .1), txdat = x, tydat = y, exdat = E,
           probabilities = TRUE, gradients = TRUE))) {
      expect_identical(fit$bws$xnames, "x")
      expect_identical(fit$bws$ynames, "y")
      expect_identical(fitted(fit), fitted(reference))
      expect_identical(fit$probabilities, reference$probabilities)
      expect_identical(gradients(fit), gradients(reference))
      expect_identical(predict(fit, newdata = E), predict(fit, exdat = E))
      expect_identical(predict(unserialize(serialize(fit, NULL)), newdata = E),
                       predict(fit, exdat = E))
      expect_identical(predict(fit, newdata = data.frame(wrong = 1), exdat = E),
                       predict(fit, exdat = E))
      expect_error(predict(fit, newdata = data.frame(wrong = 1)), "x")
    }
  }
})
