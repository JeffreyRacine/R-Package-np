test_that("CMS refits retain weights, offsets and raw rq design", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  set.seed(520)
  d <- data.frame(x = rnorm(36), z = rnorm(36))
  d$y <- 1 + d$x + .4*d$z^2 + rnorm(36)
  d$w <- exp(d$x); d$o <- .7*d$z^2
  d$positive <- exp(.2 + .1*d$x) + runif(36, .01, .08)
  zero.weights <- d$w; zero.weights[c(3, 17)] <- 0
  models <- list(
    lm = lm(y ~ x, d, x = TRUE, y = TRUE),
    weighted = lm(y ~ x, d, weights = w, x = TRUE, y = TRUE),
    offset = lm(y ~ x + offset(o), d, x = TRUE, y = TRUE),
    glm = glm(positive ~ x + offset(o/10), d, weights = w,
              family = gaussian("log"), x = TRUE, y = TRUE),
    rq = quantreg::rq(y ~ x, data = d, tau = .5, model = TRUE),
    weighted.rq = quantreg::rq(y ~ x, data = d, weights = w, tau = .5, model = TRUE),
    zero.rq = quantreg::rq(y ~ x, data = d, weights = zero.weights, tau = .5, model = TRUE))
  K <- dnorm(outer(d$x, d$x, "-")/.5)/.5; diag(K) <- 0
  statistic <- function(e, pivot) {
    I <- sum(e * (K %*% e))/nrow(d)^2
    if (!pivot) return(I)
    omega <- 2*.5*sum(e^2 * ((K^2) %*% e^2))/nrow(d)^2
    nrow(d)*sqrt(.5)*I/sqrt(omega)
  }
  for (kind in names(models)) {
    m <- models[[kind]]; original <- m
    quantile <- inherits(m, "rq")
    raw.y <- model.response(m$model)
    raw.x <- model.matrix(m$terms, m$model, contrasts.arg = m$contrasts)
    e <- as.numeric(residuals(m, type = "response"))
    for (method in c("iid", "wild", "wild-rademacher")) {
      set.seed(42)
      oracle <- replicate(9, {
        if (method == "iid") {
          ys <- fitted(m) + e[sample.int(length(e), replace = TRUE)]
        } else {
          a <- if (method == "wild") -.6180339887499 else -1
          b <- if (method == "wild") 1.6180339887499 else 1
          p <- if (method == "wild") .72360679774998 else .5
          multiplier <- rep(b, length(e)); multiplier[runif(length(e)) <= p] <- a
          center <- if (quantile) mean(e) else 0
          ys <- fitted(m) + (e-center)*multiplier + center
        }
        if (quantile) {
          refit <- suppressWarnings(quantreg::rq(ys ~ raw.x - 1, tau = .5,
                                                  weights = m$weights))
        } else {
          refit <- glm(ys ~ raw.x - 1,
            weights = if (is.null(m$family)) m$weights else m$prior.weights,
            offset = m$offset, family = if (is.null(m$family)) gaussian() else m$family)
        }
        score <- as.numeric(residuals(refit, type = "response"))
        if (quantile) score <- as.numeric(score <= 0)-.5
        c(raw = statistic(score, FALSE), pivot = statistic(score, TRUE))
      })
      for (pivot in c(FALSE, TRUE)) {
        args <- list(xdat = d["x"], ydat = raw.y, model = m, B = 9,
          bws = .5, bandwidth.compute = FALSE, pivot = pivot, boot.method = method)
        result <- if (quantile) do.call(npqcmstest, args) else do.call(npcmstest, args)
        field <- if (pivot) "Jn.bootstrap" else "In.bootstrap"
        expected <- sort(oracle[if (pivot) "pivot" else "raw", ])
        expect_equal(result[[field]], expected, tolerance = 2e-12,
                     info = paste(kind, method, pivot))
        expect_identical(m, original)
      }
    }
  }
})

test_that("CMS omission metadata and quantile ownership are explicit", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  set.seed(911); d <- data.frame(x = rnorm(25), y = rnorm(25))
  d$x[c(2,8)] <- NA
  for (quantile in c(FALSE, TRUE)) {
    m <- if (quantile) quantreg::rq(y ~ x, data = d, na.action = na.exclude) else
      lm(y ~ x, d, x = TRUE, y = TRUE, na.action = na.exclude)
    args <- list(formula = y ~ x, data = d, model = m, bws = .5,
                 bandwidth.compute = FALSE, B = 9)
    result <- if (quantile) do.call(npqcmstest, args) else do.call(npcmstest, args)
    expect_identical(as.integer(result$na.index), c(2L,8L))
  }
  m <- quantreg::rq(y ~ x, data = d, tau = .8)
  seed <- .Random.seed
  expect_error(npqcmstest(y ~ x, d, model = m, B = 9), "tau must match")
  expect_identical(.Random.seed, seed)
  expect_error(npqcmstest(y ~ x, d, model = m, tau = c(.2,.8), B = 9),
               "finite scalar")
  expect_error(npqcmstest(y ~ x, d, model = m, tau = NA_real_, B = 9),
               "finite scalar")
  expect_silent(npqcmstest(y ~ x, d, model = m, tau = .8, B = 9,
                           bws = .5, bandwidth.compute = FALSE))
})
