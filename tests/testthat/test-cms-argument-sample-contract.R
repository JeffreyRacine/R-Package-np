test_that("CMS manual bandwidths belong to selection, not duplicated kernel dots", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(319); d <- data.frame(x = rnorm(30)); d$y <- d$x^2 + rnorm(30)
  models <- list(cms = lm(y ~ x, data = d, x = TRUE, y = TRUE),
                 qcms = quantreg::rq(y ~ x, data = d, model = TRUE))
  for (kind in names(models)) {
    fun <- get(paste0("np", kind, "test"), asNamespace("npRmpi"))
    model <- models[[kind]]
    for (distribution in c("asymptotic", "bootstrap")) {
      for (kernel in c("gaussian", "epanechnikov")) {
        args <- list(formula = y ~ x, data = d, model = model, B = 9,
                     distribution = distribution, bws = .4,
                     bandwidth.compute = FALSE, ckertype = kernel)
        result <- do.call(fun, args)
        score <- as.numeric(residuals(model))
        if (kind == "qcms") score <- as.numeric(score <= 0) - .5
        z <- outer(d$x, d$x, "-")/.4
        K <- if (kernel == "gaussian") dnorm(z)/.4 else
          .75*(1-z^2/5)*(abs(z) <= sqrt(5))/(sqrt(5)*.4)
        diag(K) <- 0
        expect_equal(result$In, sum(outer(score, score)*K)/30^2, tolerance = 2e-12)
        expect_equal(result$Omega.hat, 2*.4*sum(outer(score^2, score^2)*K^2)/30^2,
                     tolerance = 2e-12)
        expect_true(is.finite(result$P))
      }
    }
  }
})
