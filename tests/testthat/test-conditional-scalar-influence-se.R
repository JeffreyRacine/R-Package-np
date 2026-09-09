test_that("conditional beta response allows a zero explanatory derivative total", {
  x <- data.frame(x = seq(.05, .95, length.out = 12L))
  y <- data.frame(y = c(.12, .65, .31, .82, .24, .57, .43, .91, .18, .73, .36, .52))
  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    estimator <- if (cdf) npcdist else npcdens
    b <- constructor(xdat = x, ydat = y, bws = c(.15, 2),
      bandwidth.compute = FALSE, regtype = "lc", cxkertype = "uniform",
      cykertype = "beta", cykerbound = "fixed", cykerlb = 0, cykerub = 1)
    fit <- estimator(se = TRUE, bws = b, txdat = x, tydat = y, exdat = data.frame(x = .5),
      eydat = data.frame(y = .43), gradients = TRUE)
    z <- as.vector(npksum(bws = .15, txdat = y, exdat = data.frame(y = .43),
      ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1,
      operator = if (cdf) "integral" else "normal", return.kernel.weights = TRUE)$kw)
    expect_equal(as.vector(fitted(fit)), mean(z), tolerance = 2e-12)
    expect_identical(as.vector(gradients(fit)), 0)
    expect_identical(as.vector(fit$congerr), 0)
    # Invalid base normalization is checked in the serial counterpart.
    # Do not exercise the separately deferred native-unwind pool cleanup here.
  }
})
test_that("scalar conditional categorical errors pair same-sample influences", {
  n <- 25L
  for (beta.x in c(FALSE, TRUE)) for (ordered in c(FALSE, TRUE)) {
    group <- factor(rep(c("a","b","c"), length.out=n), ordered=ordered)
    x <- if(beta.x) data.frame(x=seq(.05,.95,length.out=n),u=group) else data.frame(u=group)
    y <- data.frame(y=.03+.94*((seq_len(n)*7L) %% (n+1L))/(n+1L))
    ex <- if(beta.x) data.frame(x=rep(.37,3L),u=factor(c("a","b","c"),levels=levels(group),ordered=ordered))
      else data.frame(u=factor(c("a","b","c"),levels=levels(group),ordered=ordered))
    ey <- data.frame(y=rep(.43,3L))
    lambda <- if(ordered) .25 else 0
    args <- list(xdat=x,ydat=y,bws=c(.15,if(beta.x) .2,lambda),
      bandwidth.compute=FALSE,regtype="lc",oxkertype="wangvanryzin",
      cxkertype=if(beta.x) "beta" else "gaussian",cykertype=if(beta.x) "gaussian" else "beta")
    args <- c(args,if(beta.x) list(cxkerbound="fixed",cxkerlb=0,cxkerub=1)
      else list(cykerbound="fixed",cykerlb=0,cykerub=1))
    bw <- do.call(npcdensbw,args)
    fit <- npcdens(se = TRUE, bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey,gradients=TRUE)
    beta.row <- function(train,point,h) as.vector(npksum(bws=h,
      txdat=data.frame(v=train),exdat=data.frame(v=point),ckertype="beta",
      ckerbound="fixed",ckerlb=0,ckerub=1,return.kernel.weights=TRUE)$kw)
    continuous <- if(beta.x) beta.row(x$x,.37,.2) else rep(1,n)
    z <- if(beta.x) dnorm((.43-y$y)/.15)/.15 else beta.row(y$y,.43,.15)
    endpoint <- function(level) {
      distance <- abs(as.integer(group)-level)
      category <- if(ordered) ifelse(distance==0,1-lambda,.5*(1-lambda)*lambda^distance)
        else ifelse(distance==0,1-lambda,lambda/2)
      w <- continuous*category
      a <- w/sum(w); m <- sum(a*z)
      list(m=m,u=a*(z-m))
    }
    for (level in 1:3) {
      alternate <- if(ordered) if(level==1L) 2L else level-1L else 1L
      a <- endpoint(level); b <- endpoint(alternate)
      expected <- sqrt(n/(n-1)*sum((a$u-b$u)^2))
      direction <- if(ordered && level==1L) -1 else 1
      expect_equal(fit$congrad[level,ncol(x)],direction*(a$m-b$m),tolerance=3e-10)
      expect_equal(fit$congerr[level,ncol(x)],expected,tolerance=3e-10)
    }
    if(!ordered) expect_equal(unname(fit$congerr[1L,ncol(x)]),0,tolerance=3e-10)
    else expect_equal(fit$congerr[1L,ncol(x)],fit$congerr[2L,ncol(x)],tolerance=3e-10)
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
    fit <- npcdens(se = TRUE, bws=b, txdat=x, tydat=y, exdat=data.frame(x=ex),
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
