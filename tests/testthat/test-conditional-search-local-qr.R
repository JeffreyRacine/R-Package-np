test_that("conditional LP search agrees with directly deleted weighted QR", {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = TRUE)
  on.exit(options(old), add = TRUE)
  set.seed(51005)
  n <- 80L
  x <- runif(n); y <- sin(2*x) + rnorm(n, sd = .3)
  X <- data.frame(x = x); Y <- data.frame(y = y)
  grid <- data.frame(y = seq(-.2, 1.5, length.out = 9L))
  ns <- asNamespace(getNamespaceName(environment(npcdensbw)))
  hy <- .35
  for (kernel in c("gaussian", "epanechnikov"))
    for (degree in 1:2) for (hx in c(.2, 100)) {
      density <- numeric(n); loss <- numeric(n); cdf <- numeric(n)
      for (i in seq_len(n)) {
        delta <- x[-i] - x[i]; z <- delta/hx
        w <- if (kernel == "gaussian") dnorm(z) else
          pmax(0, 3/(4*sqrt(5))*(1-z*z/5))
        # The retained all-large selector admits a constant x kernel.
        if (hx == 100) w[] <- 1
        design <- outer(delta, 0:degree, `^`)
        # Solve the deleted weighted design independently in R. Its first
        # coefficient is the prediction at the held-out x, in original units.
        row <- qr.solve(sqrt(w)*design, diag(sqrt(w)))[1L, ]
        yy <- y[-i]
        density[i] <- sum(row * dnorm((y[i]-yy)/hy)/hy)
        convolution <- dnorm(outer(yy, yy, `-`)/(sqrt(2)*hy))/(sqrt(2)*hy)
        loss[i] <- 2*density[i] - drop(row %*% convolution %*% row)
        fy <- vapply(grid$y, function(g) sum(row*pnorm((g-yy)/hy)), 0)
        cdf[i] <- mean(((y[i] <= grid$y)-fy)^2)
      }
      log.score <- ifelse(density > .Machine$double.xmin, log(pmax(density, .Machine$double.xmin)),
        ifelse(density < -.Machine$double.xmin,
          2*log(.Machine$double.xmin)-log(pmax(-density, .Machine$double.xmin)),
          log(.Machine$double.xmin)))
      args <- list(xdat = X, ydat = Y, bws = c(hy, hx),
        bandwidth.compute = FALSE, regtype = "lp", degree = degree,
        bernstein.basis = FALSE, cxkertype = kernel, cykertype = "gaussian")
      for (tree in c(FALSE, TRUE)) {
        options(np.tree = tree)
        for (method in c("cv.ml", "cv.ls")) {
          b <- do.call(npcdensbw, c(args, list(bwmethod = method)))
          value <- get(".npcdensbw_eval_only", ns)(X, Y, b, invalid.penalty = "dbmax")$objective
          expect_equal(value, if (method == "cv.ml") sum(log.score) else mean(loss),
                       tolerance = if (method == "cv.ls") 5e-10 else 1e-10)
        }
        b <- do.call(npcdistbw, args)
        value <- get(".npcdistbw_eval_only", ns)(X, Y, bws = b, gydat = grid,
                                                invalid.penalty = "dbmax")$objective
        expect_equal(value, mean(cdf), tolerance = 1e-10)
      }
    }
})
