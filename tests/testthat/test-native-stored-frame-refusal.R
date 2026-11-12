test_that("stored conditional bandwidths refuse unused frame controls", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  d<-data.frame(x=seq(.02,.98,length.out=30));d$y<-sin(3*d$x)+cos(1:30)*.1
  for(fun in list(npcdensbw,npcdistbw)){
    b<-fun(xdat=d[1],ydat=d[2],bws=c(.3,.3),bandwidth.compute=FALSE)
    plain<-fun(xdat=d[1],ydat=d[2],bws=b,bandwidth.compute=FALSE)
    for (field in c("xbw", "ybw")) {
      expect_length(b[[field, exact = TRUE]], 1L)
      expect_length(plain[[field, exact = TRUE]], 1L)
      expect_equal(b[[field, exact = TRUE]], .3, tolerance = 0)
      expect_equal(plain[[field, exact = TRUE]], b[[field, exact = TRUE]], tolerance = 0)
    }
    calls<-0L;bad<-function(){calls<<-calls+1L;stop("must not force")}
    expect_error(fun(xdat=d[1],ydat=d[2],bws=b,subset=bad()),"subset=.*supported only")
    expect_error(fun(xdat=d[1],ydat=d[2],bws=b,na.action=bad()),"na.action.*formula interface")
    expect_identical(calls,0L)
    d0<-d;d0$y[5]<-NA
    calls<-0L;omit<-function(z){calls<<-calls+1L;na.omit(z)}
    f<-fun(y~x,data=d0,subset=x>.1,na.action=omit,bws=c(.3,.3),bandwidth.compute=FALSE)
    expect_identical(calls,1L);expect_equal(f$nobs,nrow(na.omit(subset(d0,x>.1))))
    set.seed(18);s<-fun(xdat=d[1],ydat=d[2],bws=b,nmulti=1L,itmax=15L)
    for (field in c("xbw", "ybw")) {
      value <- s[[field, exact = TRUE]]
      expect_length(value, 1L)
      expect_true(is.numeric(value) && all(is.finite(value)))
    }
  }
})

test_that("native model tests reject frame controls before forcing them", {
  calls<-0L;bad<-function(){calls<<-calls+1L;stop("must not force")}
  for(fun in list(npcmstest,npqcmstest)){
    expect_error(fun(xdat=1:10,ydat=1:10,subset=bad()),"subset=.*supported only")
    expect_error(fun(xdat=1:10,ydat=1:10,na.action=bad()),"na.action.*formula interface")
  }
  expect_identical(calls,0L)
})
