test_that("original higher derivatives preserve the admitted literal laws", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(61027);n <- 120L
  x <- data.frame(x=runif(n,-1,1),o=ordered(rep(1:2,n/2)))
  y <- x$x^3+.2*as.integer(x$o)+rnorm(n,sd=.3)
  b <- npregbw(xdat=x,ydat=y,bws=c(.6,.2),regtype="lp",degree=3L,
    bernstein.basis=TRUE,bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=x,tydat=y);saved <- fit
  state <- np:::.np_bias_regression_prepare(b,x,y,gradient.order=2L)
  for(law in c("inid","rademacher","fixed","geom")) {
    one <- np:::.np_bias_inference(state,x[1,,drop=FALSE],"d1",.05,63L,law)
    expect_identical(one$law,law)
    expect_identical(one$bands$pointwise,one$bands$simultaneous)
    expect_identical(one$bands$pointwise,one$bands$bonferroni)
    rng <- .Random.seed
    args <- list(x=fit,errors="bootstrap",center="bias-corrected",band="all",
      gradients=TRUE,gradient.order=2L,neval=2L,B=63L,output="data",
      bootstrap=if(law=="rademacher")"wild" else law)
    if(law=="rademacher")args$boot.control <- np_boot_control(wild="rademacher")
    out <- capture_plot_provider(do.call(plot,args))$provider
    for(panel in out) {
      z <- panel$inference
      expect_identical(z$law,law)
      expect_identical(z$settings$iterations,16L)
      expect_true(all(is.finite(z$bands$simultaneous)))
      expect_true(all(z$bands$pointwise[,2,]<=z$bands$simultaneous[,2,]))
      expect_true(all(z$bands$simultaneous[,2,]<=z$bands$bonferroni[,2,]))
    }
    expect_identical(rng,.Random.seed)
    expect_identical(fit,saved)
  }
})

test_that("pooled laws preserve scores and their admitted calibration", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(61028);n <- 120L
  x <- data.frame(x=runif(n,-1,1),u=factor(rep(1:3,each=n/3)))
  y <- x$x^3+as.integer(x$u)+rnorm(n,sd=.3)
  b <- npregbw(xdat=x,ydat=y,bws=c(.6,.2),regtype="lp",degree=3L,
    bernstein.basis=TRUE,bandwidth.compute=FALSE)
  state <- np:::.np_bias_regression_prepare(b,x,y,gradient.order=2L)
  expect_true(state$pooled)
  for(target in c("mean","d1","c2")) {
    z <- np:::.np_bias_inference(state,x[41,,drop=FALSE],target,.05,63L,"inid")
    expect_identical(z$law,"inid")
    expect_identical(z$covariance$calibration,"residual score process")
    expect_identical(z$bands$pointwise,z$bands$simultaneous)
    expect_identical(z$bands$pointwise,z$bands$bonferroni)
    expect_true(all(is.finite(z$bands$simultaneous)))
  }
  for(law in c("rademacher","fixed","geom")) {
    z <- np:::.np_bias_inference(state,x[41,,drop=FALSE],"mean",.05,63L,law)
    expect_identical(z$law,law)
    expect_identical(z$version,if(law=="rademacher")2L else 3L)
    expect_identical(z$settings$iterations,16L)
    expect_identical(z$bands$pointwise,z$bands$simultaneous)
    expect_identical(z$bands$pointwise,z$bands$bonferroni)
    expect_true(all(is.finite(z$bands$simultaneous)))
  }
})

test_that("categorical-only additional laws record their actual calibration", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(61029);x <- data.frame(u=factor(rep(1:3,each=30L)));y <- as.integer(x$u)+rnorm(90L)
  b <- npregbw(xdat=x,ydat=y,bws=.2,regtype="lc",bandwidth.compute=FALSE)
  state <- np:::.np_bias_regression_prepare(b,x,y)
  for(law in c("inid","rademacher","fixed","geom"))for(target in c("mean","c1")) {
    z <- np:::.np_bias_inference(state,x[31,,drop=FALSE],target,.05,63L,law)
    expect_identical(z$law,law)
    expect_identical(z$bands$pointwise,z$bands$simultaneous)
    expect_identical(z$bands$pointwise,z$bands$bonferroni)
    expect_identical(z$covariance$variance,if(law %in% c("fixed","geom"))
      "empirical residual-score process variance" else "within-profile residual variance")
    expect_identical(z$covariance$studentization,if(law=="rademacher")
      "within-profile residual SD" else "empirical residual-score process SD")
    expect_true(all(is.finite(z$bands$simultaneous)))
  }
})
