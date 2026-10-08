test_that("retained LSQ subsets require an explicit sample", {
  set.seed(4102); d <- data.frame(x=rnorm(80),z=runif(80)); d$y <- d$x + d$z + rnorm(80,.1)
  b <- npregbw(y~x+z,data=d,bws=c(.5,.5),bandwidth.compute=FALSE)
  expect_error(nplsqreg(b,subset=21:60,scale=rep(.3,40),bandwidth.compute=FALSE),"supply explicit 'data'")
  fit <- nplsqreg(b,data=d,subset=21:60,scale=rep(.3,40),bandwidth.compute=FALSE)
  oracle <- nplsqreg(b,txdat=d[21:60,c('x','z')],tydat=d$y[21:60],scale=rep(.3,40),bandwidth.compute=FALSE)
  expect_equal(fitted(fit),fitted(oracle),tolerance=0)
  expect_equal(fit$nobs,40)
  retained <- nplsqreg(b,scale=rep(.3,80),bandwidth.compute=FALSE)
  expect_equal(retained$nobs,80)
})

test_that("Ichimura scoring respects the factor training scale", {
  set.seed(5101); d <- data.frame(x=rnorm(90),z=rnorm(90))
  d$y <- factor(sample(0:3,90,TRUE)); e <- d[1:12,]
  b <- npindexbw(y~x+z,data=d,bws=c(1,.5,.5),bandwidth.compute=FALSE)
  f <- npindex(b,newdata=e,se=FALSE)
  expected <- mean((as.integer(e$y)-fitted(f))^2)
  expect_equal(f$MSE,expected,tolerance=1e-14)
  e$y <- as.numeric(as.character(e$y))
  unscored <- npindex(b,newdata=e,se=FALSE)
  expect_identical(unscored$diagnostics.sample,'unavailable')
  expect_true(is.na(unscored$MSE))
  expect_equal(fitted(unscored),fitted(f),tolerance=0)
  expect_error(npindex(b,newdata=e,y.eval=TRUE,se=FALSE),'must be a factor')
})

test_that("smooth coefficient predictions retain actual training MSE", {
  set.seed(5102); n <- 80L; x <- data.frame(x=rnorm(n)); z <- data.frame(z=runif(n))
  y <- (1+z$z)*x$x+rnorm(n,sd=.2)
  b <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=.5,bandwidth.compute=FALSE)
  f <- npscoef(b,txdat=x,tydat=y,tzdat=z,eydat=y+3)
  expected <- mean((y-fitted(f))^2)
  expect_equal(f$training.MSE,expected,tolerance=1e-14)
  expect_gt(abs(f$MSE-expected),1)
  pr <- predict(f,exdat=x[1:11,,drop=FALSE],ezdat=z[1:11,,drop=FALSE],se.fit=TRUE)
  expect_equal(pr$residual.scale,expected,tolerance=1e-14)
  ordinary <- npscoef(b,txdat=x,tydat=y,tzdat=z)
  expect_equal(fitted(f),fitted(ordinary),tolerance=0)
  expect_identical(ordinary$training.MSE,ordinary$MSE)
})

test_that("LSQ vector uncertainty reports a score for every tau", {
  set.seed(4103); n <- 70L; x <- data.frame(x=rnorm(n),z=runif(n)); y <- x$x+rnorm(n)
  b <- npregbw(xdat=x,ydat=y,bws=c(.5,.5),bandwidth.compute=FALSE)
  f <- nplsqreg(b,txdat=x,tydat=y,tau=c(.25,.5,.75),scale=rep(.3,n),delta=.5,bandwidth.compute=FALSE)
  pr <- predict(f,exdat=x[1:10,],se.fit=TRUE)
  expected <- vapply(f$tau.fits,function(one) one$fit$MSE,numeric(1))
  expect_equal(unname(pr$residual.scale),unname(expected),tolerance=0)
  expect_identical(names(pr$residual.scale),colnames(pr$fit))
})

test_that("plot-data surfaces and slices retain dimensions and names", {
  set.seed(321); n <- 70L
  d <- data.frame(x1=rnorm(n),x2=runif(n)); d$y <- sin(d$x1)+d$x2^2+rnorm(n,sd=.2)
  names(d)[1L] <- 'x one'
  for (family in c('npcdens','npcdist')) {
    b <- do.call(get(paste0(family,'bw')),list(xdat=d[1],ydat=d['y'],bws=c(.4,.5),bandwidth.compute=FALSE))
    fit <- get(family)(b)
    p <- plot(fit,plot.behavior='data',perspective=TRUE,neval=5)[[1L]]
    expect_identical(dim(p$xeval),c(25L,1L))
    expect_identical(dim(p$yeval),c(25L,1L))
    expect_identical(names(p$xeval),'x one')
    expect_identical(names(p$yeval),'y')
    expect_equal(p$nobs,25L)
    expect_length(fitted(p),25L)
    expect_error(capture.output(summary(p)),NA)
  }
  for (cols in list(1L,1:2)) {
    b <- npregbw(xdat=d[cols],ydat=d$y,bws=rep(.5,length(cols)),bandwidth.compute=FALSE)
    f <- npreg(b)
    for (perspective in if(length(cols)==1L) FALSE else c(FALSE,TRUE)) {
      p <- plot(f,plot.behavior='data',perspective=perspective,neval=5)[[1L]]
      expect_identical(names(p$eval),names(d)[cols])
      expect_identical(ncol(p$eval),length(cols))
      expect_equal(nrow(p$eval),length(fitted(p)))
      expect_error(capture.output(summary(p)),NA)
    }
  }
  b <- npscoefbw(xdat=d[1],ydat=d$y,zdat=d[2],bws=.4,bandwidth.compute=FALSE)
  p <- plot(npscoef(b),plot.behavior='data',perspective=FALSE,neval=5)[[1L]]
  fields <- c('R2','MSE','MAE','MAPE','CORR','SIGN')
  expect_true(all(fields %in% names(p)))
  expect_true(all(is.na(unlist(p[fields]))))
  expect_error(capture.output(summary(p)),NA)
})

test_that("LSQ dispatch uses formal promises once in their original scope", {
  set.seed(4101); n <- 70L; x <- data.frame(x=rnorm(n),z=runif(n)); y <- x$x+rnorm(n)
  b <- npregbw(xdat=x,ydat=y,bws=c(.5,.5),bandwidth.compute=FALSE)
  existed <- exists('.np_lsq_tau_test',.GlobalEnv,inherits=FALSE)
  if(existed) saved <- get('.np_lsq_tau_test',.GlobalEnv)
  on.exit(if(existed) assign('.np_lsq_tau_test',saved,.GlobalEnv) else rm('.np_lsq_tau_test',envir=.GlobalEnv),add=TRUE)
  assign('.np_lsq_tau_test',.9,.GlobalEnv)
  call <- function(.np_lsq_tau_test) nplsqreg(b,txdat=x,tydat=y,tau=.np_lsq_tau_test,
    scale=rep(.3,n),delta=.5,bandwidth.compute=FALSE)
  expect_identical(call(.25)$tau,.25)
  count <- 0L
  fit <- nplsqreg(b,txdat=x,tydat=y,scale={count<-count+1L;rep(.3,n)},delta=.5,bandwidth.compute=FALSE)
  expect_identical(count,1L)
  values <- list(bws=c(.3,.3),bandwidth.compute=FALSE)
  first <- nplsqregbw(b,xdat=x,ydat=y,pilot.args=values,delta=.5,bandwidth.compute=FALSE)
  alias <- values
  second <- nplsqregbw(b,xdat=x,ydat=y,pilot.args=alias,delta=.5,bandwidth.compute=FALSE)
  expect_identical(first$scale,second$scale)
})

test_that("native factor-index bandwidths retain data without caller frames", {
  set.seed(2521); n <- 80L
  d <- data.frame(x=rnorm(n),z=runif(n),g=factor(sample(c('A','B','C'),n,TRUE)))
  d$y <- sin(d$x)+d$z+rnorm(n,sd=.2)
  make <- function() {
    unused <- rnorm(2e5)
    npindexbw(xdat=d[1:3],ydat=d$y,bws=c(1,.5,.4,-.6,.3),bandwidth.compute=FALSE)
  }
  b <- make()
  expect_null(environment(b$call))
  expect_lt(length(serialize(b,NULL)),100000L)
  cold <- unserialize(serialize(b,NULL))
  expect_equal(fitted(npindex(cold,se=FALSE)),fitted(npindex(b,se=FALSE)),tolerance=0)
})

test_that("positional native refits replace training observations", {
  set.seed(260); n <- 80L
  x <- data.frame(x=rnorm(n)); y <- sin(x$x)+rnorm(n,sd=.2)
  X <- data.frame(x=rev(x$x)+.7); Y <- cos(X$x)+rnorm(n,sd=.2)
  b <- npregbw(xdat=x,ydat=y,bws=.5,bandwidth.compute=FALSE)
  positional <- npreg(b,X,Y)
  named <- npreg(b,txdat=X,tydat=Y)
  expect_equal(fitted(positional),fitted(named),tolerance=0)
  expect_true(positional$trainiseval)
  expect_gt(max(abs(fitted(positional)-fitted(npreg(b)))),.1)
  expect_equal(fitted(npreg(b,X,tydat=Y)),fitted(named),tolerance=0)
  E <- X[1:9,,drop=FALSE]
  expect_equal(predict(named,se.fit=TRUE,E),
               predict(named,se.fit=TRUE,exdat=E),tolerance=0)
  expect_equal(predict(named,FALSE,E,newdata=data.frame(wrong=1)),
               predict(named,exdat=E),tolerance=0)
  for (pair in list(list(npudensbw,npudens),list(npudistbw,npudist))) {
    bw <- pair[[1]](dat=x,bws=.5,bandwidth.compute=FALSE)
    fit <- pair[[2]](bw,X)
    expect_equal(fitted(fit),fitted(pair[[2]](bw,tdat=X)),tolerance=0)
    expect_equal(predict(fit,FALSE,E),predict(fit,edat=E),tolerance=0)
  }
})

test_that("smooth-coefficient positional refits agree with weighted least squares", {
  set.seed(261); n <- 80L
  x <- data.frame(x=rnorm(n)); z <- data.frame(z=runif(n)); y <- x$x+z$z+rnorm(n)
  X <- data.frame(x=rev(x$x)); Z <- data.frame(z=rev(z$z)); Y <- cos(X$x)-Z$z
  b <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=.5,bandwidth.compute=FALSE)
  fit <- npscoef(b,X,Y,Z)
  W <- cbind(1,X$x)
  expected <- vapply(seq_len(n),function(i) {
    w <- dnorm((Z$z-Z$z[i])/.5)
    sum(W[i,]*solve(crossprod(W,W*w),crossprod(W,w*Y)))
  },0.0)
  expect_equal(fitted(fit),expected,tolerance=1e-12)
})

test_that("argument matching preserves values and the native method order", {
  value <- quote(stop('must remain data'))
  args <- .np_match_native_args(list(1,2,value),npreghat.rbandwidth)
  expect_identical(names(args),c('txdat','exdat','y'))
  expect_identical(args$y,value)
  args <- .np_match_native_args(list(1,eydat=NULL),npreg.rbandwidth,c('txdat','tydat'))
  expect_identical(names(args),c('exdat','eydat'))
  expect_null(args$eydat)
  args <- .np_match_native_args(list(txd=1,2),npreg.rbandwidth)
  expect_identical(names(args),c('txdat','tydat'))
})
