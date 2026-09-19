test_that("smooth-coefficient row outputs restore the evaluation sample together", {
  old <- options(np.messages = FALSE, na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  set.seed(190904)
  d <- data.frame(x = runif(36), z = runif(36))
  d$y <- (1 + d$z)*d$x + rnorm(36, sd = .2)
  d$x[c(3, 9)] <- NA
  e <- d[c(1, 2, 4, 5, 6), ]
  e$z[c(2, 4)] <- NA
  for (reg in c("lc", "ll", "lp")) for (bt in c("fixed", "generalized_nn", "adaptive_nn")) {
    args <- list(formula = y ~ x | z, data = d,
      bws = if (bt == "fixed") .65 else 20, bwtype = bt,
      bandwidth.compute = FALSE, regtype = reg)
    if (reg == "lp") args$degree <- 2
    bw <- do.call(npscoefbw, args)
    for (external in c(FALSE, TRUE)) for (ss in c(FALSE, TRUE)) {
      fit.args <- list(bws = bw, betas = TRUE, se = ss, iterate = FALSE)
      if (external) fit.args$newdata <- e
      a <- do.call(npscoef, fit.args)
      omitted <- if (external) c(2L, 4L) else c(3L, 9L)
      nr <- if (external) 5L else 36L
      expect_identical(nrow(coef(a)), nr)
      expect_identical(nrow(a$grad), nr)
      expect_identical(which(is.na(coef(a)[, 1])), omitted)
      expect_identical(which(is.na(a$grad[, 1])), omitted)
      if (ss) {
        expect_identical(nrow(a$gerr), nr)
        expect_identical(which(is.na(a$gerr[, 1])), omitted)
      } else expect_null(a$gerr)
      expect_equal(unname(coef(a)[, -1, drop = FALSE]), unname(a$grad))
    }
    one <- npscoef(bw, newdata = e[1, , drop = FALSE], betas = TRUE, se = TRUE)
    expect_identical(dim(coef(one)), c(1L, 2L))
    expect_identical(dim(one$gerr), c(1L, 1L))
    none <- npscoef(bw, newdata = e, betas = FALSE)
    expect_identical(none$beta, NA)
  }
})
