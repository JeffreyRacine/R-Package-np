test_that("retained training samples are never spliced into partial replacements", {
  set.seed(267); n <- 80L
  x <- data.frame(x=runif(n)); y <- sin(4*x$x)+rnorm(n,sd=.2)
  X <- data.frame(x=runif(n)); Y <- cos(4*X$x)+rnorm(n,sd=.2)
  b <- npregbw(xdat=x,ydat=y,bws=.08,bandwidth.compute=FALSE)
  expect_error(npreg(b,X), "training data 'tydat' missing")
  expect_error(npreg(b,tydat=Y), "training data 'txdat' missing")
  fit <- npreg(b,X,Y)
  oracle <- vapply(X$x,function(e){w<-dnorm((X$x-e)/.08);sum(w*Y)/sum(w)},0.0)
  expect_equal(fitted(fit),oracle,tolerance=1e-12)
  E <- X[1:11,,drop=FALSE]
  pred <- vapply(E$x,function(e){w<-dnorm((x$x-e)/.08);sum(w*y)/sum(w)},0.0)
  expect_equal(fitted(npreg(b,exdat=E)),pred,tolerance=1e-12)
  expect_equal(predict(npreg(b),se.fit=FALSE,E),pred,tolerance=1e-12)
  # Hat application replaces a right-hand side, not the smoother's design.
  expect_equal(as.vector(npreghat(b,y=Y,output='apply')),
    as.vector(npreghat(b,txdat=x,y=Y,output='apply')),tolerance=0)
  z <- data.frame(z=runif(n)); Z <- data.frame(z=runif(n))
  bs <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=.2,bandwidth.compute=FALSE)
  expect_equal(fitted(npscoef(bs,X,Y)),fitted(npscoef(bs,txdat=X,tydat=Y)),tolerance=0)
  W <- cbind(1,X$x)
  sw <- vapply(seq_len(n),function(i){w<-dnorm((X$x-X$x[i])/.2);sum(W[i,]*solve(crossprod(W,W*w),crossprod(W,w*Y)))},0.0)
  expect_equal(fitted(npscoef(bs,X,Y)),sw,tolerance=1e-10)
  bp <- npplregbw(xdat=x,ydat=y,zdat=z,bws=matrix(.2,2,1),bandwidth.compute=FALSE)
  expect_error(npplreg(bp,X,Y),"training data 'tzdat' missing")
  for (pair in list(list(npcdensbw,npcdens),list(npcdistbw,npcdist))) {
    bc <- pair[[1]](xdat=x,ydat=data.frame(y=y),bws=c(.2,.2),bandwidth.compute=FALSE)
    expect_error(pair[[2]](bc,X),"training data 'tydat' missing")
    expect_equal(fitted(pair[[2]](bc,X,data.frame(y=Y))),
      fitted(pair[[2]](bc,txdat=X,tydat=data.frame(y=Y))),tolerance=0)
  }
})

test_that("LSQ positional formula prediction preserves names and transformations", {
  set.seed(268); n <- 80L
  d <- data.frame(x=runif(n,.2,3),z=rnorm(n));d$y<-log(d$x)+.3*d$z+rnorm(n,sd=.2)
  for (tau in list(.5,c(.25,.75))) {
    b <- npregbw(y~log(x)+z,data=d,bws=c(.3,.4),bandwidth.compute=FALSE)
    f <- nplsqreg(b,tau=tau,scale=rep(.3,n),delta=.5,bandwidth.compute=FALSE)
    E <- d[1:11,c('z','x')];E$x <- E$x+.1
    expect_equal(predict(f,se.fit=TRUE,E),predict(f,newdata=E,se.fit=TRUE),tolerance=0)
    native <- data.frame('log(x)'=log(E$x),z=E$z,check.names=FALSE)
    expect_equal(predict(f,se.fit=TRUE,E),predict(f,se.fit=TRUE,exdat=native),tolerance=1e-12)
    expect_equal(predict(f,se.fit=TRUE,E,exdat=native),predict(f,se.fit=TRUE,exdat=native),tolerance=0)
    expect_equal(NROW(predict(f,se.fit=TRUE,E)$fit),11L)
    expect_error(predict(f,se.fit=TRUE,E['z']),"newdata must contain")
  }
})

test_that("LSQ training controls cannot silently reuse frozen pseudo responses", {
  set.seed(269);n<-80L;d<-data.frame(x=rnorm(n),z=runif(n));d$y<-sin(d$x)+d$z+rnorm(n,sd=.2)
  b<-npregbw(y~x+z,data=d,bws=c(.5,.5),bandwidth.compute=FALSE)
  f<-nplsqreg(b,scale=rep(.3,n),delta=.5,bandwidth.compute=FALSE)
  for(tau in list(.5,c(.25,.75))) {
    lb<-nplsqreg(b,tau=tau,scale=rep(.3,n),delta=.5,bandwidth.compute=FALSE)$bws
    expect_equal(fitted(nplsqreg(lb,data=NULL,subset=TRUE)),fitted(nplsqreg(lb)),tolerance=0)
    expect_equal(fitted(nplsqreg(lb,subset=NULL)),fitted(nplsqreg(lb)),tolerance=0)
    expect_error(nplsqreg(lb,data=d),"stored LSQ bandwidth object")
    expect_error(nplsqreg(lb,subset=21:60),"stored LSQ bandwidth object")
    expect_error(nplsqreg(lb,na.action=na.omit),"stored LSQ bandwidth object")
    expect_equal(fitted(nplsqreg(lb)),fitted(nplsqreg(lb,txdat=lb$xdat,tydat=lb$ydat)),tolerance=0)
  }
  args<-list(bws=b,scale=rep(.3,n),delta=.5,bandwidth.compute=FALSE)
  expect_equal(fitted(do.call(nplsqreg,c(args,list(subset=NULL)))),fitted(f),tolerance=0)
  expect_equal(fitted(do.call(nplsqreg,c(args,list(subset=TRUE)))),fitted(f),tolerance=0)
  counter<-new.env();counter$n<-0L
  g<-nplsqreg(b,data=d,subset={counter$n<-counter$n+1L;21:60},scale=rep(.3,40),delta=.5,bandwidth.compute=FALSE)
  expect_identical(counter$n,1L);expect_equal(g$bws$ydat,d$y[21:60],tolerance=0)
  g<-nplsqreg(b,data=d,subset=z>.5,scale=rep(.3,sum(d$z>.5)),delta=.5,bandwidth.compute=FALSE)
  expect_equal(g$bws$ydat,d$y[d$z>.5],tolerance=0)
  set.seed(5);idx<-sample(n,40)
  set.seed(5);g<-nplsqreg(b,data=d,subset=sample(n,40),scale=rep(.3,40),delta=.5,bandwidth.compute=FALSE)
  expect_equal(g$bws$ydat,d$y[idx],tolerance=0)
})

test_that("index explicit diagnostics preserve the training response scale", {
  set.seed(263);n<-90L;d<-data.frame(x=rnorm(n),z=rnorm(n));latent<-d$x+.5*d$z+rnorm(n,sd=.5)
  y<-factor(cut(latent,quantile(latent,0:4/4),include.lowest=TRUE,labels=FALSE)-1L,levels=0:3)
  for (ordered in c(FALSE,TRUE)) {
    yy<-if(ordered)ordered(y,levels=levels(y))else y
    b<-npindexbw(xdat=d,ydat=yy,bws=c(1,.5,.4),bandwidth.compute=FALSE)
    E<-d[1:15,];ey<-yy[1:15];numeric.ey<-as.numeric(as.character(ey))
    f<-npindex(b,exdat=E,eydat=ey,se=FALSE)
    training <- data.frame(y=yy,d)
    bf <- npindexbw(y~x+z,data=training,bws=c(1,.5,.4),bandwidth.compute=FALSE)
    eval <- data.frame(y=numeric.ey,E)
    auto <- npindex(bf,newdata=eval,se=FALSE)
    expect_equal(fitted(auto),fitted(f),tolerance=1e-12)
    if(ordered) expect_equal(auto$MSE,f$MSE,tolerance=1e-12) else expect_true(is.na(auto$MSE))
    oracle<-mean(((if(ordered)numeric.ey else as.integer(ey))-fitted(f))^2)
    expect_equal(f$MSE,oracle,tolerance=1e-12)
    if(ordered) {
      expect_equal(npindex(b,exdat=E,eydat=numeric.ey,y.eval=TRUE,se=FALSE)$MSE,oracle,tolerance=1e-12)
    } else {
      expect_error(npindex(b,exdat=E,eydat=numeric.ey,y.eval=TRUE,se=FALSE),"must be a factor")
      expect_error(npindex(b,exdat=E,eydat=as.character(ey),se=FALSE),"must be a numeric vector or factor")
    }
  }
})
