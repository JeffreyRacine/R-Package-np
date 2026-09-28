test_that("sparse kernel powers and deletion follow the initialized support", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.02, .98, length.out = 40L))
  y <- sin(3*x$x)
  h <- .08
  u <- outer(x$x, x$x, "-")/h
  K <- ifelse(abs(u) < sqrt(5), 3/(4*sqrt(5))*(1-u^2/5), 0)
  for (tree in c(FALSE, TRUE)) for (power in 1:3) for (loo in c(FALSE, TRUE)) {
    options(np.tree = tree)
    W <- K^power
    if (loo) diag(W) <- 0
    got <- npksum(txdat = x, tydat = y, bws = h, bwscaling = FALSE,
                  ckertype = "epanechnikov", ckerorder = 2L,
                  kernel.pow = power, leave.one.out = loo)$ksum
    expect_equal(as.numeric(got), as.numeric(crossprod(y, W)), tolerance = 1e-12)
  }
})

test_that("conditional LP uncertainty retains sparse and dense agreement", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(3)
  x <- data.frame(x = runif(40L))
  y <- data.frame(y = x$x + rnorm(40L, sd = .2))
  for (cdf in c(FALSE, TRUE)) for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    make <- if (cdf) npcdistbw else npcdensbw
    fit <- if (cdf) npcdist else npcdens
    b <- make(xdat = x, ydat = y, bws = if (type == "fixed") c(.2, .2) else c(8, 8),
              bwtype = type, bwscaling = FALSE, regtype = "lp", degree = 1L,
              cxkertype = "epanechnikov", cykertype = "epanechnikov",
              bandwidth.compute = FALSE)
    for (external in c(FALSE, TRUE)) {
      a <- list(bws = b, txdat = x, tydat = y, gradients = TRUE,
                gradient.order = 1L, se = TRUE)
      if (external) {
        a$exdat <- data.frame(x = c(.25, .5, .75))
        a$eydat <- data.frame(y = c(.2, .6, .8))
      }
      options(np.tree = FALSE)
      dense <- do.call(fit, a)
      options(np.tree = TRUE)
      sparse <- do.call(fit, a)
      for (extract in list(fitted, se, gradients, function(z) gradients(z, se = TRUE))) {
        want <- extract(dense)
        got <- extract(sparse)
        expect_true(all(is.finite(want)) && all(is.finite(got)))
        expect_equal(got, want, tolerance = 2e-11)
      }
    }
  }
})
