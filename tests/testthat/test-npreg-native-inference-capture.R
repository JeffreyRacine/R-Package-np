test_that("captured global inference matches independent polynomial OLS", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=TRUE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=64L))
  y <- x$x+x$x^2+cos(seq_len(64L))*.2
  e <- data.frame(x=c(-.7,-.3,0,.3,.7))
  bw <- npregbw(xdat=x,ydat=y,bws=1000,bandwidth.compute=FALSE,
    regtype="lp",degree=2,bernstein.basis=TRUE)
  fit <- .np_regression_direct(bw,x,y,e,se=TRUE,gradients=TRUE,fit.inference=TRUE)
  ordinary <- .np_regression_direct(bw,x,y,e,se=TRUE,gradients=TRUE)
  expect_identical(fit[names(ordinary)],ordinary)
  Q <- cbind(1,x$x,x$x^2); q <- qr(Q); U <- qr.Q(q)
  A <- backsolve(qr.R(q),t(U))
  z <- (y-drop(Q%*%qr.coef(q,y)))/sqrt(1-rowSums(U^2))
  for(derivative in c(FALSE,TRUE)) {
    V <- if(derivative)cbind(0,1,2*e$x)else cbind(1,e$x,e$x^2)
    H <- V%*%A; C <- tcrossprod(sweep(H,2L,z,"*"))
    got <- .np_plot_regression_joint_covariance(bw,x,y,e,fit=fit,
      deriv=if(derivative)1L else NULL,block.bytes=64*64*2)
    expect_equal(got$covariance,C,tolerance=1e-12)
    expect_equal(sqrt(diag(C)),if(derivative)fit$gerr[,1L]else fit$merr,
                 tolerance=1e-12)
    expect_gt(max(abs(C-diag(diag(C)))),1e-6)
  }
  expect_error(.np_plot_regression_joint_covariance(bw,x,y+1,e,fit=fit),"fit request")
  expect_error(.np_plot_regression_joint_covariance(bw,x,y,e,fit=fit,deriv=2L),
               "derivative order")
})

test_that("captured partial-global rows preserve request scope and permutations", {
  old <- options(np.messages=FALSE,np.largeh=TRUE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=48L),z=sin(seq_len(48L)))
  y <- x$x*x$z+cos(seq_len(48L))*.2
  e <- x[c(2,11,24,32,47),,drop=FALSE]
  for(tree in list(FALSE,TRUE,"auto")) {
    options(np.tree=tree)
    bw <- npregbw(xdat=x,ydat=y,bws=c(.45,1000),bandwidth.compute=FALSE,
      regtype="lp",degree=c(1,2),bernstein.basis=TRUE)
    f <- .np_regression_direct(bw,x,y,e,se=TRUE,gradients=TRUE,fit.inference=TRUE)
    a <- .np_plot_regression_joint_covariance(bw,x,y,e,fit=f,block.bytes=64*48)
    b <- .np_plot_regression_joint_covariance(bw,x,y,e,fit=f)
    expect_equal(a$covariance,b$covariance,tolerance=1e-13)
    expect_equal(sqrt(diag(a$covariance)),f$merr,tolerance=1e-12)
    ix <- rev(seq_len(nrow(x))); ie <- c(5L,2L,1L,4L,3L)
    fp <- .np_regression_direct(bw,x[ix,,drop=FALSE],y[ix],e[ie,,drop=FALSE],
      se=TRUE,gradients=TRUE,fit.inference=TRUE)
    cp <- .np_plot_regression_joint_covariance(bw,x[ix,,drop=FALSE],y[ix],
      e[ie,,drop=FALSE],fit=fp)
    expect_equal(cp$covariance,a$covariance[ie,ie],tolerance=1e-12)
  }
})

test_that("large-bandwidth public bands retain marginal intervals and nesting", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-2,2,length.out=60L)); y <- x$x+cos(seq_len(60L))
  bw <- npregbw(xdat=x,ydat=y,bws=700,regtype="ll",bandwidth.compute=FALSE)
  for(grad in c(FALSE,TRUE)) {
    a <- plot(bw,errors="asymptotic",band="all",gradients=grad,
      output="data",neval=5L)[[1L]]
    p <- plot(bw,errors="asymptotic",band="pointwise",gradients=grad,
      output="data",neval=5L)[[1L]]
    expect_identical(a$mean,p$mean); expect_identical(a$merr,p$merr)
    bands <- a$band.errors
    expect_true(all(is.finite(bands$simultaneous)))
    expect_true(all(bands$pointwise<=bands$simultaneous))
    expect_true(all(bands$simultaneous<=bands$bonferroni))
  }
})

test_that("uncaptured contrasts cannot inherit captured residual state", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=48L),f=factor(rep(letters[1:3],16)))
  y <- x$x+as.numeric(x$f)*.2+cos(seq_len(48L))*.1
  bw <- npregbw(xdat=x,ydat=y,bws=c(.4,.2),regtype="ll",bandwidth.compute=FALSE)
  original <- .np_plot_regression_joint_covariance
  contrast.residuals <- list()
  testthat::local_mocked_bindings(.np_plot_regression_joint_covariance=function(...) {
    args <- list(...)
    if (!is.null(args$categorical.index))
      contrast.residuals[[length(contrast.residuals)+1L]] <<- list(value=args$residual)
    original(...)
  }, .package="np")
  plot(bw,errors="asymptotic",band="all",gradients=TRUE,
       output="data",neval=5L)
  expect_length(contrast.residuals,1L)
  expect_null(contrast.residuals[[1L]]$value)
})
