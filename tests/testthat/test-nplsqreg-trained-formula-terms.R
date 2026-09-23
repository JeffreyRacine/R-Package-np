test_that("LSQ formula owners retain trained transformations at every tau", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  center <- function(x, location = mean(x))
    structure(x-location, location=location, class=c("numeric","r21lsqcenter"))
  make.center <- function(var, call) { call$location <- attr(var,"location"); call }
  previous <- getS3method("makepredictcall","r21lsqcenter",optional=TRUE)
  registerS3method("makepredictcall","r21lsqcenter",make.center,envir=asNamespace("stats"))
  on.exit({
    if(is.null(previous)) rm("makepredictcall.r21lsqcenter",
      envir=get(".__S3MethodsTable__.",asNamespace("stats"))) else
      registerS3method("makepredictcall","r21lsqcenter",previous,envir=asNamespace("stats"))
  },add=TRUE)
  d <- data.frame(x=seq(10,20,length.out=19))
  d$y <- sin(d$x)+.1*cos(seq_len(nrow(d)))
  d$y[4L] <- NA_real_
  e <- data.frame(x=c(11,12,18))
  formula <- y ~ center(x)
  tt <- attr(model.frame(formula,d,subset=x>11),"terms")
  ex <- model.frame(delete.response(tt),e)
  for(route in c("onecall","bandwidth")) for(tau in list(.3,c(.3,.7))) {
    aa <- list(bws=formula,data=d,scale=rep(1,nrow(d)),bw=.4,delta=.4,
      tau=tau,bandwidth.compute=FALSE,subset=quote(x>11))
    fit <- if(route=="onecall") do.call(nplsqreg,aa) else nplsqreg(do.call(nplsqregbw,aa))
    expect_identical(attr(fit$bws$terms,"predvars"),attr(tt,"predvars"))
    expect_identical(nrow(fit$bws$xdat),16L)
    p <- predict(fit,newdata=e)
    native <- predict(fit,exdat=ex)
    expect_equal(p,native,tolerance=0)
    expect_equal(predict(fit,newdata=e[1:2,,drop=FALSE]),
      if(is.matrix(p))p[1:2,,drop=FALSE]else p[1:2],tolerance=0)
    expect_equal(predict(unserialize(serialize(fit,NULL)),newdata=e),p,tolerance=0)
    expect_equal(predict(fit,newdata=data.frame(wrong=1),exdat=ex),native,tolerance=0)
    expect_error(predict(fit,newdata=data.frame(wrong=1)),"must contain")
    children <- if(length(tau)>1L)fit$tau.fits else list(fit)
    for(child in children) {
      expect_identical(attr(child$bws$terms,"predvars"),attr(tt,"predvars"))
      expect_equal(predict(child,newdata=e),predict(child,exdat=ex),tolerance=0)
    }
    direct <- do.call(nplsqreg,c(aa,list(newdata=e)))
    expect_equal(fitted(direct),native,tolerance=0)
    direct.native <- do.call(nplsqreg,c(aa,list(exdat=ex,newdata=data.frame(wrong=1))))
    expect_equal(fitted(direct.native),native,tolerance=0)
  }
})

test_that("LSQ trained terms preserve ordinary and omitted-row contracts", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  d<-data.frame(x=seq(.3,2,length.out=19));d$y<-sin(d$x)
  e<-data.frame(x=c(.4,NA_real_,1.7))
  for(ff in list(y~x,y~log(x))) {
    fit<-nplsqreg(bws=ff,data=d,scale=rep(1,nrow(d)),bw=.4,delta=.4,
      tau=.3,bandwidth.compute=FALSE)
    legacy<-fit;legacy$bws$terms<-NULL
    expect_equal(predict(fit,newdata=e),predict(legacy,newdata=e),tolerance=0)
    expect_true(is.na(predict(fit,newdata=e)[2L]))
    expect_identical(length(predict(fit,newdata=e)),3L)
  }
})
