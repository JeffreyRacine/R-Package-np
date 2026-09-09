test_that("conditional beta response allows a zero explanatory derivative total", {
  x <- data.frame(x = seq(.05, .95, length.out = 12L))
  y <- data.frame(y = c(.12, .65, .31, .82, .24, .57, .43, .91, .18, .73, .36, .52))
  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    estimator <- if (cdf) npcdist else npcdens
    b <- constructor(xdat = x, ydat = y, bws = c(.15, 2),
      bandwidth.compute = FALSE, regtype = "lc", cxkertype = "uniform",
      cykertype = "beta", cykerbound = "fixed", cykerlb = 0, cykerub = 1)
    fit <- estimator(bws = b, txdat = x, tydat = y, exdat = data.frame(x = .5),
      eydat = data.frame(y = .43), gradients = TRUE)
    z <- as.vector(npksum(bws = .15, txdat = y, exdat = data.frame(y = .43),
      ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      operator = if (cdf) "integral" else "normal", return.kernel.weights = TRUE)$kw)
    expect_equal(as.vector(fitted(fit)), mean(z), tolerance = 2e-12)
    expect_identical(as.vector(gradients(fit)), 0)
    expect_identical(as.vector(fit$congerr), 0)
    expect_error(suppressWarnings(estimator(bws = b, txdat = x, tydat = y,
      exdat = data.frame(x = 10), eydat = data.frame(y = .43), gradients = TRUE)),
      "canonical beta response-unit restoration failed")
  }
})

test_that("scalar conditional beta influence errors use sample covariance scaling", {
  for (n in c(12L, 25L)) for (beta.x in c(FALSE, TRUE)) {
    x <- data.frame(x = seq(.05, .95, length.out = n))
    y <- data.frame(y = .03 + .94*((seq_len(n)*7L) %% (n+1L))/(n+1L))
    ex <- .37; ey <- .43; hx <- .2; hy <- .15
    side <- function(t, e, h, beta) {
      if (!beta) return(dnorm((e-t)/h)/h)
      as.vector(npksum(bws=h, txdat=data.frame(v=t), exdat=data.frame(v=e),
        ckertype="beta", ckerbound="fixed", ckerlb=0, ckerub=1,
        return.kernel.weights=TRUE)$kw)
    }
    a <- list(xdat=x, ydat=y, bws=c(hy,hx), bandwidth.compute=FALSE,
      regtype="lc", cxkertype=if(beta.x) "beta" else "gaussian",
      cykertype=if(beta.x) "gaussian" else "beta")
    a <- c(a, if(beta.x) list(cxkerbound="fixed",cxkerlb=0,cxkerub=1)
      else list(cykerbound="fixed",cykerlb=0,cykerub=1))
    b <- do.call(npcdensbw, a)
    fit <- npcdens(bws=b, txdat=x, tydat=y, exdat=data.frame(x=ex),
      eydat=data.frame(y=ey), gradients=TRUE)
    w <- side(x$x, ex, hx, beta.x); z <- side(y$y, ey, hy, !beta.x)
    alpha <- w/sum(w); m <- sum(alpha*z)
    step <- 2e-6
    wp <- side(x$x, ex+step, hx, beta.x); wm <- side(x$x, ex-step, hx, beta.x)
    ap <- (wp/sum(wp)-wm/sum(wm))/(2*step)
    g <- sum(ap*z)
    u <- alpha*(z-m); v <- ap*(z-m)-alpha*g
    expect_equal(as.vector(fitted(fit)), m, tolerance=3e-10)
    expect_equal(as.vector(se(fit)), sqrt(n/(n-1)*sum(u*u)), tolerance=3e-10)
    expect_equal(as.vector(gradients(fit)), g, tolerance=3e-7)
    expect_equal(as.vector(fit$congerr), sqrt(n/(n-1)*sum(v*v)), tolerance=3e-7)
  }
})
