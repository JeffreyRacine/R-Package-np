test_that("CMS inference uses the compact model sample with na.exclude", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(319); d <- data.frame(x = rnorm(30)); d$y <- d$x^2 + rnorm(30)
  d$x[c(4, 19)] <- NA
  for (kind in c("cms", "qcms", "glm")) {
    make <- function(policy) switch(kind,
      cms = lm(y ~ x, data = d, x = TRUE, y = TRUE, na.action = policy),
      glm = glm(y ~ x, data = d, x = TRUE, y = TRUE, na.action = policy),
      qcms = quantreg::rq(y ~ x, data = d, model = TRUE, na.action = policy))
    omit <- make(na.omit); exclude <- make(na.exclude)
    immutable <- c("x", "y", "residuals", "fitted.values", "na.action", "coefficients", "call")
    before <- unserialize(serialize(exclude[immutable], NULL))
    fun <- get(if (kind == "qcms") "npqcmstest" else "npcmstest", asNamespace("np"))
    for (distribution in c("asymptotic", "bootstrap")) {
      for (route in c("formula", "native")) {
        inputs <- if (route == "formula") list(formula = y ~ x, data = d) else
          list(xdat = d$x, ydat = d$y)
        args <- c(inputs, list(distribution = distribution, B = 9,
          bws = .4, bandwidth.compute = FALSE))
        a <- do.call(fun, c(args, list(model = omit)))
        b <- do.call(fun, c(args, list(model = exclude)))
        fields <- setdiff(names(a), c("pcall", "bws", "timing.profile"))
        expect_identical(a[fields], b[fields])
      }
    }
    expect_identical(exclude[immutable], before)
    expect_error(do.call(fun, list(xdat = rev(d$x), ydat = rev(d$y), model = exclude,
      bws = .4, bandwidth.compute = FALSE)), "same complete observations")
  }
})

test_that("CMS validates retained components rather than call spelling", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(527); d <- data.frame(x = rnorm(24)); d$y <- d$x + rnorm(24)
  keep <- TRUE
  models <- list(
    cms = list(lm(y ~ x, d, x = TRUE, y = TRUE), lm(y ~ x, d, x = keep, y = keep)),
    qcms = list(quantreg::rq(y ~ x, data = d, model = TRUE),
                quantreg::rq(y ~ x, data = d, model = keep)))
  for (kind in names(models)) {
    fun <- get(paste0("np", kind, "test"), asNamespace("np"))
    invoke <- function(model) do.call(fun, list(formula = y ~ x, data = d,
      model = model, bws = .4, bandwidth.compute = FALSE, distribution = "asymptotic"))
    a <- invoke(models[[kind]][[1L]]); b <- invoke(models[[kind]][[2L]])
    fields <- setdiff(names(a), c("pcall", "bws", "timing.profile"))
    expect_identical(a[fields], b[fields])
    broken <- models[[kind]][[1L]]; broken$x <- NULL
    expect_error(invoke(broken), "must retain")
  }
})

test_that("CMS manual bandwidths belong to selection, not duplicated kernel dots", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(319); d <- data.frame(x = rnorm(30)); d$y <- d$x^2 + rnorm(30)
  models <- list(cms = lm(y ~ x, data = d, x = TRUE, y = TRUE),
                 qcms = quantreg::rq(y ~ x, data = d, model = TRUE))
  for (kind in names(models)) {
    fun <- get(paste0("np", kind, "test"), asNamespace("np"))
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

test_that("CMS mixed-kernel contractions retain selected categorical metadata", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(645)
  d <- data.frame(x = rnorm(30), u = factor(rep(c("a","b","c"),10)),
                  o = ordered(rep(1:3,10)), y = rnorm(30))
  model <- lm(y~x+u+o, data=d, x=TRUE, y=TRUE)
  score <- residuals(model)
  K <- dnorm(outer(d$x,d$x,"-")/.4)/.4 *
    (ifelse(outer(as.character(d$u),as.character(d$u),"=="),1,.2)/(1+2*.2)) *
    .5^abs(outer(as.integer(d$o),as.integer(d$o),"-"))
  diag(K) <- 0
  for (weighted in c(FALSE,TRUE)) for (pivot in c(FALSE,TRUE)) {
    result <- npcmstest(y~x+u+o, data=d, model=model, bws=c(.4,.2,.5),
      bandwidth.compute=FALSE, ukertype="liracine", okertype="liracine",
      B=9, density.weighted=weighted, pivot=pivot)
    fhat <- if(weighted) rep(1,30) else rowSums(K)/30
    expect_equal(result$In,sum(score*(K%*%score)/fhat)/30^2,tolerance=2e-12)
    if(pivot) expect_equal(result$Omega.hat,
      2*.4*sum(score^2*((K^2)%*%(score^2))/fhat^2)/30^2,tolerance=2e-12)
  }
})
