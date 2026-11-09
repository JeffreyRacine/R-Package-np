test_that("LSQ formula scale-factor controls are not auxiliary scale vectors", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 24L
  data <- data.frame(x = seq(.05, .95, length.out = n))
  data$y <- sin(3 * data$x) + cos(9 * data$x) / 7
  for (control in list(list(scale.factor.search.lower = .05),
                       list(scale.factor.init = .5))) {
    args <- c(list(bws = y ~ x, data = data, bw = .2, delta = .5,
                   regtype = "ll", regtype.pilot = "ll",
                   bandwidth.compute = FALSE,
                   pilot.args = list(nmulti = 1L, itmax = 10L,
                                     powell.remin = FALSE)), control)
    for (fun in list(nplsqregbw, nplsqreg)) {
      set.seed(882)
      actual <- do.call(fun, args)
      args.null <- args
      args.null["scale"] <- list(NULL)
      set.seed(882)
      expected <- do.call(fun, args.null)
      if (inherits(actual, "lsqregressionbandwidth")) {
        expect_identical(actual$bw, expected$bw)
        expect_identical(actual$objective, expected$objective)
        expect_identical(actual$scale, expected$scale)
      } else expect_identical(fitted(actual), fitted(expected))
    }
  }
})
