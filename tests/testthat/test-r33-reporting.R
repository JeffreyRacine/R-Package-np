# Small metadata-only contracts for two sibling prediction methods.
test_that("smooth coefficient predictions retain known training MSE", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(445); x <- data.frame(x=runif(80)); z <- data.frame(z=runif(80)); y <- sin(x$x)+z$z+rnorm(80,sd=.1)
  b <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=.3,bandwidth.compute=FALSE)
  f <- npscoef(b,txdat=x,tydat=y,tzdat=z,se=TRUE)
  p <- predict(f,exdat=x[1:15,,drop=FALSE],ezdat=z[1:15,,drop=FALSE],eydat=y[1:15]+3,se.fit=TRUE)
  expect_equal(p$residual.scale,mean((y-fitted(f))^2))
  expect_identical(p$fit,predict(f,exdat=x[1:15,,drop=FALSE],ezdat=z[1:15,,drop=FALSE]))
})

test_that("LSQ predictions report training transformed-response MSE", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(448); x <- data.frame(x=runif(80)); y <- sin(x$x)+rnorm(80,sd=.1)
  b <- npregbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE)
  f <- nplsqreg(bws=b,tau=.5,txdat=x,tydat=y,se=TRUE)
  p <- predict(f,exdat=x[1:15,,drop=FALSE],se.fit=TRUE)
  expect_equal(p$residual.scale,f$fit$MSE)
})


test_that("partially linear plot grids do not fabricate goodness-of-fit scores", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(91); x<-data.frame(x=runif(80));z<-data.frame(z=runif(80))
  y<-2*x$x+sin(z$z)+rnorm(80,sd=.1)
  b<-npplregbw(xdat=x,zdat=z,ydat=y,bws=matrix(.3,2,1),bandwidth.compute=FALSE)
  f<-npplreg(bws=b,se=TRUE)
  for(common in c(TRUE,FALSE)) {
    p<-plot(f,plot.behavior="data",neval=8,errors="none",common.scale=common)
    for(one in p) expect_true(all(is.na(unlist(one[c("R2","MSE","MAE","MAPE","CORR")]))))
  }
})
