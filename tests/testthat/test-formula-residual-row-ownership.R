test_that("formula residuals retain training omissions with external evaluation", {
  old <- options(np.messages = FALSE, na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  set.seed(190905)
  d <- data.frame(x = runif(36, -1, 1), z = runif(36, -1, 1))
  d$y <- sin(d$x) + .2*d$z + rnorm(36, sd = .2)
  d$x[c(3, 9)] <- NA
  e <- d[c(1, 2, 4, 5, 6), ]
  e$z[c(2, 4)] <- NA
  for (policy in c("na.exclude", "na.omit")) {
    options(na.action = policy)
    for (family in c("npreg", "npindex", "npscoef")) {
      f <- if (family == "npscoef") y ~ x | z else y ~ x + z
      bw <- do.call(get(paste0(family, "bw")), list(formula = f,
        data = d, bws = switch(family, npreg = c(.65, .65),
          npindex = c(1, .4, .65), npscoef = .65),
        bandwidth.compute = FALSE, regtype = "lc", na.action = policy))
      train <- get(family)(bws = bw, residuals = TRUE, se = TRUE)
      for (explicit in c(FALSE, TRUE)) {
        args <- list(bws = bw, newdata = e, residuals = TRUE, se = TRUE)
        if (explicit) args$na.action <- get(policy, asNamespace("stats"))
        fit <- do.call(get(family), args)
        expect_equal(unname(residuals(fit)), unname(residuals(train)), tolerance = 1e-12)
        expect_identical(unname(which(is.na(residuals(fit)))),
                         if (policy == "na.exclude") c(3L, 9L) else integer())
        expect_length(fitted(fit), if (policy == "na.exclude") 5L else 3L)
      }
    }
  }
})
