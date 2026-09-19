test_that("stored uncertainty availability is independent of numeric NA values", {
  specs <- list(npdensity="derr", npdistribution="derr", condensity="conderr",
    condistribution="conderr", qregression="quanterr", smoothcoefficient="merr",
    lsqregression="quanterr")
  for (cls in names(specs)) {
    for (value in list(0, .2, NA_real_, NaN, c(NA_real_, .2), matrix(NA_real_,1,1))) {
      obj <- structure(setNames(list(value, TRUE), c(specs[[cls]], "se")), class=cls)
      expect_identical(se(obj), value)
      obj$se <- FALSE
      expect_error(se(obj), "standard errors were not computed")
    }
    for (value in list(NULL, numeric(0), NA)) {
      obj <- structure(setNames(list(value, TRUE), c(specs[[cls]], "se")), class=cls)
      expect_error(se(obj), "standard errors were not computed")
    }
  }
})

test_that("single unsupported smooth-coefficient predictions return computed NA errors", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(190984)
  d <- data.frame(x=runif(32),z=runif(32),y=rnorm(32))
  for (rt in c("lc","ll","lp")) {
    a <- list(formula=y~x|z,data=d,bws=.4,regtype=rt,
              ckertype="epanechnikov",bandwidth.compute=FALSE)
    if(rt=="lp") a$degree <- 2
    b <- do.call(npscoefbw,a)
    f <- npscoef(b,iterate=FALSE)
    for(native in c(FALSE,TRUE)) {
      args <- if(native) list(exdat=data.frame(x=9),ezdat=data.frame(z=9)) else
        list(newdata=data.frame(x=9,z=9))
      ev <- suppressWarnings(do.call(npscoef,c(list(bws=b,se=TRUE,iterate=FALSE),args)))
      expect_identical(ev$se,TRUE)
      expect_identical(se(ev),ev$merr)
      p <- suppressWarnings(do.call(predict,c(list(object=f,se.fit=TRUE),args)))
      expect_identical(p$fit,NA_real_)
      expect_identical(p$se.fit,NA_real_)
    }
    batch <- suppressWarnings(predict(f,newdata=data.frame(x=c(.4,9),z=c(.4,9)),se.fit=TRUE))
    expect_true(is.finite(batch$fit[1]))
    expect_true(is.finite(batch$se.fit[1]))
    expect_identical(batch$se.fit[2],NA_real_)
  }
})

test_that("LSQ scalar and vector tau extract computed unsupported-row uncertainty", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(190985)
  d <- data.frame(x=runif(32),y=rnorm(32))
  for(tau in list(.5,c(.3,.7))) {
    b <- nplsqregbw(bws=y~x,data=d,bw=.4,delta=.3,scale=rep(1,32),
      regtype="ll",ckertype="epanechnikov",bandwidth.compute=FALSE,tau=tau)
    f <- suppressWarnings(nplsqreg(b,exdat=data.frame(x=9),se=TRUE,gradients=TRUE))
    expect_identical(se(f),f$quanterr)
    expect_true(all(is.na(se(f))))
    expect_identical(se(unserialize(serialize(f,NULL))),se(f))
    if(length(tau)==1L) {
      rebuilt <- lsqregression(bws=f$bws,fit=f$fit,xeval=f$xeval,
        tau=f$tau,delta=f$delta,quantile=f$quantile,quanterr=f$quanterr,
        ntrain=f$ntrain)
      expect_identical(rebuilt$se,TRUE)
      expect_identical(se(rebuilt),f$quanterr)
    }
    expect_true(all(is.na(gradients(f,se=TRUE))))
    off <- suppressWarnings(nplsqreg(b,exdat=data.frame(x=9),se=FALSE))
    expect_error(se(off),"standard errors were not computed")
  }
})
