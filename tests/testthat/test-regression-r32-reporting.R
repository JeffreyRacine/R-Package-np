test_that("regression NN admission uses current complete sample size", {
  old<-options(np.messages=FALSE,np.extendednn=FALSE);on.exit(options(old),add=TRUE)
  set.seed(185);x<-data.frame(x=rnorm(90));y<-rnorm(90)
  for(type in c("generalized_nn","adaptive_nn")) {
    b<-npregbw(xdat=x,ydat=y,bws=85,bwtype=type,bandwidth.compute=FALSE)
    xn<-x;xn$x[1:10]<-NA
    expect_error(npregbw(xdat=xn,ydat=y,bws=b,bandwidth.compute=FALSE),"exceeds n-1")
    options(np.extendednn=TRUE)
    r<-npregbw(xdat=xn,ydat=y,bws=b,bandwidth.compute=FALSE)
    expect_equal(r$nobs,80L);expect_equal(r$bw,b$bw)
    options(np.extendednn=FALSE)
  }
})

test_that("predictor-only regression scores are unavailable without changing fits", {
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  set.seed(187);x<-data.frame(x=rnorm(90),z=rnorm(90));y<-sin(x$x)+.2*x$z
  for(reg in c("lc","ll","lp")) for(tree in c(FALSE,TRUE)) {
    options(np.tree=tree)
    args<-list(xdat=x,ydat=y,bws=c(.6,.7),regtype=reg,bandwidth.compute=FALSE)
    if(reg=="lp")args$degree<-c(2L,2L)
    b<-do.call(npregbw,args)
    f<-npreg(b,se=TRUE,gradients=TRUE)
    a<-npreg(b,exdat=x[1:20,],se=TRUE,gradients=TRUE)
    c<-npreg(b,exdat=x[1:20,],eydat=y[1:20],se=TRUE,gradients=TRUE)
    fields<-c("R2","MSE","MAE","MAPE","CORR","SIGN")
    expect_true(all(is.na(unlist(a[fields]))))
    expect_output(summary(a), "Regression Data")
    expect_output(print(a), "Regression Data")
    expect_identical(fitted(a),fitted(c))
    expect_identical(se(a),se(c))
    expect_identical(gradients(a),gradients(c))
    expect_equal(c$MSE,mean((y[1:20]-fitted(c))^2))
    expect_equal(f$MSE,mean((y-fitted(f))^2))
    expect_true(is.na(predict(f,exdat=x[1:20,],se.fit=TRUE)$residual.scale))
  }
})
