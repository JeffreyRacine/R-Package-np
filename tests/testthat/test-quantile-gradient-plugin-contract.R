test_that("quantile gradients retain the documented selected-density target", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920320)
  x <- data.frame(x=runif(40,.02,.98))
  y <- data.frame(y=pmin(.96,pmax(.03,.1+.65*x$x+rnorm(40,sd=.07))))
  ex <- data.frame(x=c(.2,.5,.8))
  for (kernel in c("gaussian","beta")) {
    for (bound in if(kernel=="beta") "fixed" else c("none","fixed")) {
      for (type in c("fixed","adaptive_nn")) {
        args <- list(txdat=x,tydat=y,exdat=ex,
          bws=if(type=="fixed") c(.25,.25) else c(10,12),bwtype=type,
          cykertype=kernel,cykerbound=bound,regtype="lc",se=FALSE)
        if(bound=="fixed") args <- c(args,list(cykerlb=0,cykerub=1))
        fit <- do.call(npqreg,c(args,list(tau=.4,gradients=TRUE)))
        ey <- data.frame(y=as.double(fitted(fit)))
        cdf <- npcdist(bws=fit$bws,txdat=x,tydat=y,exdat=ex,eydat=ey,
                      gradients=TRUE,se=FALSE)
        dens <- do.call(npcdens,c(args,list(eydat=ey)))
        expect_true(all(fitted(dens)>0))
        oracle <- -as.double(gradients(cdf))/as.double(fitted(dens))
        expect_equal(as.double(gradients(fit)),oracle,tolerance=1e-12)
      }
    }
  }
})
