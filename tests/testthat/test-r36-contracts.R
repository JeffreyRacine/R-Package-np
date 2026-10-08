test_that("smooth coefficient refits retain the z=x convention they used", {
  set.seed(3601); n <- 90L
  x <- data.frame(x=runif(n)); z <- data.frame(z=runif(n))
  y <- sin(4*x$x)+rnorm(n,sd=.2)
  make <- function() {
    xl <- x; zl <- z; yl <- y
    npscoefbw(xdat=xl,ydat=yl,zdat=zl,bws=.3,bandwidth.compute=FALSE)
  }
  b <- make(); X <- data.frame(x=runif(n)); Y <- cos(3*X$x)+rnorm(n,sd=.2)
  W <- cbind(1,X$x)
  oracle <- function(E,Z=E) vapply(seq_len(nrow(E)),function(i) {
    w <- dnorm((X$x-Z[[1]][i])/.3)
    sum(c(1,E$x[i])*solve(crossprod(W,W*w),crossprod(W,w*Y)))
  },0.0)
  E <- X[1:11,,drop=FALSE]
  for (messages in c(FALSE, TRUE)) for (positional in c(FALSE,TRUE)) {
    withr::local_options(np.messages=messages)
    f <- if(positional) npscoef(b,X,Y) else npscoef(b,txdat=X,tydat=Y)
    expect_equal(fitted(f),oracle(X),tolerance=1e-11)
    expect_equal(fitted(npscoef(b,txdat=X,tydat=Y,tzdat=NULL,newdata=E)),oracle(E),tolerance=1e-11)
    expect_null(f$bws$zdati); expect_null(f$bws$znames); expect_null(f$bws$varnames$z)
    expect_equal(fitted(npscoef(f$bws)),oracle(X),tolerance=1e-11)
    expect_equal(predict(f,exdat=E),oracle(E),tolerance=1e-11)
    expect_equal(predict(f,exdat=E,ezdat=E),oracle(E),tolerance=1e-11)
    expect_equal(predict(f,newdata=E),oracle(E),tolerance=1e-11)
    expect_equal(predict(unserialize(serialize(f,NULL)),newdata=E),oracle(E),tolerance=1e-11)
  }
  # Explicit separate smoothing data retains its own role.
  f <- npscoef(b,txdat=X,tydat=Y,tzdat=z)
  expect_identical(f$bws$znames,b$znames)
  expect_equal(predict(f,exdat=X,ezdat=z),fitted(f),tolerance=1e-11)
  # A formula with a separate z continues to retain all its formula roles.
  d <- data.frame(x=x$x,z=z$z,y=y)
  bf <- npscoefbw(y~x|z,data=d,bws=.3,bandwidth.compute=FALSE)
  ff <- npscoef(bf)
  expect_equal(predict(ff,newdata=d),fitted(ff),tolerance=1e-11)
})

test_that("forwarded LSQ subsets retain their data mask and evaluate once", {
  set.seed(3620);n<-80L;d<-data.frame(x=runif(n),z=runif(n));d$y<-sin(4*d$x)+rnorm(n,sd=.2)
  b<-npregbw(y~x,data=d,bws=.2,bandwidth.compute=FALSE)
  wrapper<-function(...)nplsqreg(b,data=d,...,delta=.5,bandwidth.compute=FALSE)
  outer<-function(...)wrapper(...)
  for (fun in list(wrapper,outer)) {
    counter<-new.env();counter$n<-0L
    f<-fun(subset={counter$n<-counter$n+1L;x>.4},scale=rep(.3,n))
    keep<-d$x>.4
    oracle<-vapply(d$x[keep],function(e){w<-dnorm((d$x[keep]-e)/.2);sum(w*d$y[keep])/sum(w)},0.0)
    expect_identical(counter$n,1L);expect_equal(fitted(f),oracle,tolerance=1e-12)
  }
})

test_that("fit selectors use the first realization of forwarded controls", {
  set.seed(3621);n<-80L;x<-data.frame(x=runif(n));z<-data.frame(z=runif(n));y<-sin(4*x$x)+rnorm(n,sd=.2)
  for (tree in c(FALSE,TRUE)) {
    withr::local_options(np.tree=tree)
    counter<-new.env();counter$n<-0L
    f<-npreg(txdat=x,tydat=y,nmulti={counter$n<-counter$n+1L;1L},itmax=3L)
    expect_identical(counter$n,1L)
    counter$n<-0L
    f<-npreg(txdat=x,tydat=y,nomad={counter$n<-counter$n+1L;FALSE},nmulti=1L,itmax=3L)
    expect_identical(counter$n,1L)
  }
  counter<-new.env();counter$n<-0L
  f<-npscoef(txdat=x,tydat=y,tzdat=z,nmulti={counter$n<-counter$n+1L;1L})
  expect_identical(counter$n,1L)
})

test_that("native entry values and response labels survive dispatch", {
  withr::local_options(np.messages=FALSE)
  set.seed(3622);n<-80L;X<-data.frame(x=runif(n),z=rnorm(n));Y<-sin(4*X$x)+.2*X$z+rnorm(n,sd=.2)
  counter<-new.env();counter$n<-0L;s_m<-seq_len(40L);c_m<-'native call label'
  ref<-npregbw(xdat=X,ydat=Y,bws=c(.2,.3),bandwidth.compute=FALSE)
  b<-npregbw(xdat=X,ydat=Y,bws=c(.2,.3),bandwidth.compute=FALSE,subset={counter$n<-counter$n+1L;s_m})
  expect_identical(counter$n,1L);expect_identical(b$bw,ref$bw)
  # Native subset remains an accepted dot, not a formula row-selection request.
  set.seed(3623);a<-npregbw(xdat=X,ydat=Y,nmulti=1L,itmax=3L,subset=s_m)
  set.seed(3623);b<-npregbw(xdat=X,ydat=Y,nmulti=1L,itmax=3L)
  expect_identical(a$bw,b$bw)
  set.seed(3624);a<-npudensbw(dat=X,nmulti=1L,call=c_m)
  set.seed(3624);b<-npudensbw(dat=X,nmulti=1L)
  expect_identical(a$bw,b$bw)
  make<-function(){response<-Y;npindexbw(xdat=X,ydat=response,bws=c(1,.2,.3),bandwidth.compute=FALSE)}
  b<-make();expect_identical(b$ynames,'response')
  E<-data.frame(X,response=Y)
  for(obj in list(b,unserialize(serialize(b,NULL)))) {
    f<-npindex(obj,newdata=E,se=FALSE)
    expect_identical(f$diagnostics.sample,'evaluation')
    expect_equal(f$MSE,mean((Y-fitted(f))^2),tolerance=1e-13)
  }
})
