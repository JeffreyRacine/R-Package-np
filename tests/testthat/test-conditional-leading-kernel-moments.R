test_that("conditional leading SE uses separate X and Y kernel moments", {
  n <- 41L
  x <- data.frame(x1 = seq(-.4, .4, length.out = n),
                  x2 = seq(-.2, .2, length.out = n))
  y <- data.frame(y = .4 * sin(seq_len(n)))
  ex <- data.frame(x1 = c(-.03, .02), x2 = c(.01, -.02))
  ey <- data.frame(y = c(.02, -.01))
  moments <- c(gaussian = 1/(2*sqrt(pi)), uniform = .5)
  for (cdf in c(FALSE, TRUE)) for (kx in names(moments))
    for (ky in names(moments)) {
      bwfun <- if(cdf) npcdistbw else npcdensbw
      fitfun <- if(cdf) npcdist else npcdens
      bw <- bwfun(xdat=x, ydat=y, bws=rep(.8,3),
                  bandwidth.compute=FALSE, regtype="lc",
                  cxkertype=kx, cykertype=ky)
      fit <- fitfun(bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey)
      wx <- vapply(seq_len(nrow(ex)), function(q) {
        products <- lapply(seq_len(ncol(x)), function(j) {
          u <- (x[[j]]-ex[[j]][q])/.8
          if(kx=="gaussian") dnorm(u) else .5*(abs(u)<1)
        })
        sum(Reduce(`*`,products))
      },numeric(1))
      value <- as.double(fitted(fit))
      variance <- if(cdf) value*(1-value)*moments[[kx]]^2/wx else
        value*moments[[kx]]^2*moments[[ky]]/(.8*wx)
      expect_equal(as.double(se(fit)),sqrt(variance),tolerance=2e-12)
    }
})

test_that("conditional analytic derivative SE uses derivative kernel energy", {
  x <- data.frame(x = seq(-.8, .8, length.out = 47L))
  y <- data.frame(y = .3*sin(seq_len(nrow(x))))
  ex <- data.frame(x = c(-.037, .021))
  ey <- data.frame(y = c(.01, -.02))
  for (family in c("gaussian", "epanechnikov")) for (order in c(2L,4L)) {
    radius <- if(family == "gaussian") Inf else sqrt(5)
    kernel <- if(family == "gaussian") {
      if(order == 2L) dnorm else function(u) dnorm(u)*(1.5-.5*u^2)
    } else {
      if(order == 2L) function(u) 3/(4*sqrt(5))*(1-u^2/5) else
        function(u) 15/(32*sqrt(5))*(3-2*u^2+.28*u^4)
    }
    deriv <- if(family == "gaussian") {
      if(order == 2L) function(u) -u*dnorm(u) else
        function(u) u*dnorm(u)*(-2.5+.5*u^2)
    } else {
      if(order == 2L) function(u) -3/(10*sqrt(5))*u else
        function(u) 15/(32*sqrt(5))*(-4*u+1.12*u^3)
    }
    ratio <- sqrt(integrate(function(u) deriv(u)^2,-radius,radius)$value /
      integrate(function(u) kernel(u)^2,-radius,radius)$value)
    for (cdf in c(FALSE,TRUE)) for (type in c("fixed","generalized_nn")) {
      bwfun <- if(cdf) npcdistbw else npcdensbw
      fitfun <- if(cdf) npcdist else npcdens
      bw <- bwfun(xdat=x,ydat=y,bws=if(type=="fixed") c(.8,.8) else c(15,15),
        bandwidth.compute=FALSE,bwtype=type,regtype="lc",
        cxkertype=family,cxkerorder=order)
      fit <- fitfun(bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey,gradients=TRUE)
      h <- if(type=="fixed") rep(.8,nrow(ex)) else
        vapply(ex$x,function(q) sort(abs(x$x-q))[15L],numeric(1))
      expect_equal(as.double(fit$congerr),as.double(se(fit))*ratio/h,
                   tolerance=3e-8)
    }
  }
})
