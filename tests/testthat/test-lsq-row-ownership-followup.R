test_that("LSQ predictions preserve evaluation rows independently of training", {
  old <- options(np.messages=FALSE,na.action="na.exclude")
  on.exit(options(old),add=TRUE)
  set.seed(190926)
  d <- data.frame(x=runif(40),y=rnorm(40)); d$x[c(3,9)] <- NA
  e <- data.frame(x=c(.1,NA,.5,.8,NA))
  for (tau in list(.3,c(.7,.3))) {
    f <- nplsqreg(bws=y~x,data=d,scale=rep(1,40),bw=.5,delta=.4,
                  tau=tau,bandwidth.compute=FALSE,residuals=TRUE,se=TRUE)
    for (s in c(FALSE,TRUE)) {
      p <- predict(f,newdata=e,se.fit=s)
      control <- predict(f,newdata=e[c(1,3,4),,drop=FALSE],se.fit=s)
      vals <- if(s) list(p$fit,p$se.fit) else list(p)
      controls <- if(s) list(control$fit,control$se.fit) else list(control)
      for(j in seq_along(vals)) {
        expect_equal(NROW(vals[[j]]),5L)
        v <- as.matrix(vals[[j]]); ref <- as.matrix(controls[[j]])
        expect_true(all(is.na(v[c(2,5),,drop=FALSE])))
        expect_equal(unname(v[c(1,3,4),,drop=FALSE]),unname(ref),tolerance=1e-12)
      }
      native <- predict(f,exdat=e,newdata=data.frame(wrong=0),se.fit=s)
      expect_equal(native,p,tolerance=1e-12)
    }
    expect_equal(NROW(residuals(f)),40L)
    expect_true(all(is.na(as.matrix(residuals(f))[c(3,9),,drop=FALSE])))
    expect_equal(unname(as.matrix(residuals(f))[-c(3,9),,drop=FALSE]),
      unname(d$y[-c(3,9)]-as.matrix(fitted(f))[-c(3,9),,drop=FALSE]),tolerance=1e-12)
    expect_equal(predict(f),fitted(f),tolerance=1e-12)
    if(length(tau)>1L) {
      legacy <- f
      for(j in seq_along(legacy$tau.fits))
        legacy$tau.fits[[j]]$bws$omit <- NULL
      expect_equal(predict(legacy),fitted(f),tolerance=1e-12)
    }
    expect_equal(predict(unserialize(serialize(f,NULL)),newdata=e),
                 predict(f,newdata=e),tolerance=1e-12)
    options(na.action="na.pass")
    expect_equal(predict(f,newdata=e),predict(f,exdat=e),tolerance=1e-12)
    options(na.action="na.exclude")
    expect_error(residuals(nplsqreg(f$bws)), "not available")
  }
})
test_that("retained LSQ bandwidths preserve only their own training omissions", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(190927)
  d <- data.frame(x=runif(40),y=rnorm(40)); d$x[c(3,9)]<-NA
  for (tau in list(.3,c(.7,.3))) for (formula in c(FALSE,TRUE)) {
    args <- list(scale=rep(1,40),bw=.5,delta=.4,tau=tau,bandwidth.compute=FALSE)
    if(formula) args <- c(list(bws=y~x,data=d,na.action=na.exclude),args) else
      args <- c(list(xdat=d["x"],ydat=d$y),args)
    bw <- do.call(nplsqregbw,args)
    for (external in c(FALSE,TRUE)) {
      opts <- list(bws=bw,residuals=TRUE,gradients=TRUE,se=TRUE)
      if(external) opts$exdat <- data.frame(x=c(.1,NA,.8))
      f <- do.call(nplsqreg,opts)
      expect_equal(NROW(fitted(f)),if(external) 3L else 40L)
      expect_equal(NROW(se(f)),if(external) 3L else 40L)
      expect_equal(dim(gradients(f))[1L],if(external) 3L else 40L)
      expect_equal(NROW(residuals(f)),40L)
      expect_true(all(is.na(as.matrix(residuals(f))[c(3,9),,drop=FALSE])))
      if(length(tau)>1L) for(one in f$tau.fits)
        expect_equal(NROW(residuals(one)),40L)
      expect_equal(NROW(predict(f)),40L)
    }
    # Explicit compact data owns compact rows, not the constructor's original map.
    g <- nplsqreg(bw,txdat=bw$xdat,tydat=bw$ydat,residuals=TRUE)
    expect_equal(NROW(fitted(g)),38L)
    expect_equal(NROW(predict(g)),38L)
    expect_equal(NROW(fitted(nplsqreg(g$bws))),38L)
  }
})

test_that("LSQ omission maps survive subset and single-row evaluation", {
  old <- options(np.messages=FALSE,na.action="na.exclude")
  on.exit(options(old),add=TRUE)
  set.seed(190950)
  d <- data.frame(x=seq(.1,.9,length.out=40),y=rnorm(40))
  d$x[c(3,9)] <- NA
  for(tau in list(.3,c(.7,.3))) {
    f <- nplsqreg(bws=y~I(x^2),data=d,subset=seq_len(35),
      scale=rep(1,40),bw=.5,delta=.4,tau=tau,bandwidth.compute=FALSE,
      residuals=TRUE)
    expect_equal(NROW(fitted(f)),35L)
    expect_equal(NROW(residuals(f)),35L)
    e <- data.frame(x=c(NA,.5,NA))
    p <- predict(f,newdata=e,se.fit=TRUE)
    expect_equal(NROW(p$fit),3L)
    expect_equal(NROW(p$se.fit),3L)
    expect_true(all(is.na(as.matrix(p$fit)[c(1,3),,drop=FALSE])))
    options(na.action="na.fail")
    expect_error(predict(f,newdata=e),"missing values")
    options(na.action="na.exclude")
    g <- nplsqreg(bws=y~x,data=d,scale=rep(1,40),bw=.4,delta=.4,
      tau=tau,ckertype="epanechnikov",regtype="ll",bandwidth.compute=FALSE)
    warnings <- character()
    p <- withCallingHandlers(
      predict(g,newdata=data.frame(x=c(.5,NA,9)),se.fit=TRUE),
      warning=function(w) {
        warnings <<- c(warnings,conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_equal(NROW(p$fit),3L)
    expect_true(all(is.na(as.matrix(p$fit)[2:3,,drop=FALSE])))
    expect_length(warnings,1L)
    expect_match(warnings[[1]],"row\\(s\\) \\(3\\)")
  }
})
