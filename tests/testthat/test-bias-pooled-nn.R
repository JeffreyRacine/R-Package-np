test_that("pooled NN rows retain the fitted owner and polynomial functionals", {

  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(61003);n <- 120L
  x <- data.frame(x=runif(n,-1,1),u=factor(rep(1:3,each=40)),o=ordered(rep(1:2,60)))
  signal <- x$x^3+.2*as.integer(x$u);y <- signal+rnorm(n,sd=.3)
  for(type in c("adaptive_nn","generalized_nn")) {
    b <- npregbw(xdat=x,ydat=y,bws=c(60,.2,.2),bwtype=type,regtype="lp",degree=3,
      bernstein.basis=TRUE,bandwidth.compute=FALSE)
    a <- np:::.np_bias_regression_prepare(b,x,y,gradient.order=2L)
    g <- x[1:2,,drop=FALSE];g$x <- c(-.3,.3);g$u <- factor(2,levels=levels(x$u))
    v <- np:::.np_bias_regression_rows(a,g,expand.grid(row=1:2,target=c("mean","d1","c2"),stringsAsFactors=FALSE))
    fit <- npreg(b,txdat=x,tydat=y,exdat=g,gradients=TRUE,gradient.order=2L)
    expect_equal(v$original_center[1:4],c(fitted(fit),gradients(fit)[,1L]),tolerance=1e-9)
    expect_equal(drop(v$W %*% signal),c(g$x^3+.4,6*g$x,rep(.2,2)),tolerance=1e-8)
    expect_identical(a$settings$iterations,16L)
    expect_identical(a$settings$main$type,type)
  }
})
