test_that("external smooth-coefficient residuals are training residuals independent of se", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(190908)
  d <- data.frame(x = runif(30), z = runif(30))
  d$y <- (1 + d$z)*d$x + rnorm(30, sd = .2)
  e <- d[1:5, ]
  for (reg in c("lc", "ll", "lp")) for (bt in c("fixed", "generalized_nn", "adaptive_nn")) {
    hs <- if (bt == "fixed") c(.4, 1e8) else 20
    for (h in hs) {
      args <- list(formula = y ~ x | z, data = d, bws = h, bwtype = bt,
        regtype = reg, bandwidth.compute = FALSE)
      if (reg == "lp") args$degree <- 2
      bw <- do.call(npscoefbw, args)
      train <- npscoef(bw, iterate = FALSE)
      expected <- d$y - fitted(train)
      for (native in c(FALSE, TRUE)) {
        a <- list(bws = bw, iterate = FALSE, betas = TRUE, residuals = TRUE)
        if (native) {
          a <- c(a, list(txdat = d["x"], tzdat = d["z"], tydat = d$y,
                        exdat = e["x"], ezdat = e["z"], eydat = e$y))
        } else a <- c(a, list(newdata = e, y.eval = TRUE))
        off <- do.call(npscoef, c(a, list(se = FALSE)))
        on <- do.call(npscoef, c(a, list(se = TRUE)))
        expect_type(residuals(off), "double")
        expect_equal(as.vector(residuals(off)), as.vector(expected), tolerance = 1e-11)
        expect_equal(residuals(off), residuals(on), tolerance = 1e-12)
        expect_identical(fitted(off), fitted(on))
        expect_identical(off$beta, on$beta)
        expect_null(off$merr)
        expect_null(off$gerr)
      }
    }
  }
})
