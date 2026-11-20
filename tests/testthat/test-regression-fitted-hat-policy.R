test_that("hat matrices and apply preserve the fitted accepted ridge policy", {
  if(!spawn_mpi_slaves())skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 120L; i <- seq_len(n)
  x <- data.frame(x=seq(-1,1,length.out=n),z=sin(i*sqrt(2)))
  y <- cos(i*.37)+.1*sin(i*sqrt(2))
  h <- exp(seq(log(.02),log(.4),length.out=40L))[6L]
  bw <- npregbw(xdat=x,ydat=y,bws=rep(h,2L),regtype="lp",
    degree=c(2L,2L),bernstein.basis=TRUE,bandwidth.compute=FALSE)
  fit <- .np_regression_direct(bw,x,y,x,se=TRUE,gradients=TRUE,fit.inference=TRUE)
  for(j in 0:2) {
    s <- as.integer(seq_len(2L)==j)
    H <- npreghat(bw,txdat=x,exdat=x,s=s)
    target <- if(j)fit$grad[,j]else fit$mean
    expect_identical(attr(H,"ridge.used"),fit$fit.inference[[3L]])
    expect_equal(drop(H%*%y),target,tolerance=2e-11)
    yy <- cbind(response=y,other=rev(y))
    app <- npreghat(bw,txdat=x,exdat=x,s=s,y=yy,output="apply")
    expect_identical(colnames(app),colnames(yy))
    expect_equal(unname(app),unname(H%*%yy),tolerance=2e-11)
    if(j) {
      factor <- .npRmpi_with_local_regression(.np_plot_lsq_derivative_factor(bw,x,y,x,coordinate=j,
        estimate=target,se=fit$gerr[,j]))
      expect_equal(factor$estimate,target,tolerance=2e-11)
      expect_equal(factor$A,sweep(H,2L,factor$residual,"*"),
        tolerance=2e-11,ignore_attr=TRUE)
    }
  }
})

test_that("lean row capture leaves the existing request unchanged", {
  if(!spawn_mpi_slaves())skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=40L),z=sin(seq_len(40L)))
  y <- sin(x$x)+x$z
  bw <- npregbw(xdat=x,ydat=y,bws=c(.2,.2),regtype="lp",degree=c(2L,2L),
    bernstein.basis=TRUE,bandwidth.compute=FALSE)
  capture <- function(...) .npRmpi_with_local_regression(.np_plot_lsq_gradient_projection(...))
  legacy <- capture(bw,x,y,x)
  lean <- capture(bw,x,y,x,capture.gradient=TRUE)
  expect_identical(lean[names(legacy)],legacy)
  full <- .np_regression_direct(bw,x,y,x,se=TRUE,gradients=TRUE,fit.inference=TRUE)
  expect_identical(lean$ridge,full$fit.inference[[3L]])
  expect_identical(lean$gate,full$fit.inference[[4L]])
  expect_identical(lean$gradient.projection,full$fit.inference[[1L]])
  mean <- capture(bw,x,y,x,capture.gradient=TRUE,capture.mean.only=TRUE)
  expect_null(mean$gradient.projection)
  expect_identical(mean$mean,lean$mean)
  expect_identical(mean$ridge,lean$ridge)
})

test_that("fitted hat capture preserves constant and empty external rows", {
  if(!spawn_mpi_slaves())skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=64L));y<-sin(x$x)
  for(h in c(.01,1e6)) {
    bw<-npregbw(xdat=x,ydat=y,bws=h,regtype="lp",degree=2L,
      bernstein.basis=TRUE,bandwidth.compute=FALSE)
    ex<-if(h<1)data.frame(x=c(-1,0,1,1e6))else x[c(1L,32L,64L),,drop=FALSE]
    for(s in 0:2) {
      H <- suppressWarnings(npreghat(bw,txdat=x,exdat=ex,s=s))
      app <- suppressWarnings(npreghat(bw,txdat=x,exdat=ex,s=s,y=cbind(y,y),output="apply"))
      expect_equal(unname(app),unname(H%*%cbind(y,y)),tolerance=2e-11)
      if(h<1) {
        expect_true(all(is.na(H[4L,])))
        expect_true(all(is.na(app[4L,])))
      } else {
        fit<-.np_regression_direct(bw,x,y,ex,se=FALSE,gradients=s>0,
          gradient.order=max(1L,s))
        expect_equal(drop(H%*%y),if(s)fit$grad[,1L]else fit$mean,tolerance=2e-11)
        expect_true(all(attr(H,"ridge.used")>=0))
      }
    }
  }
})

test_that("higher derivatives and LSQ factors use the same fitted map", {
  if(!spawn_mpi_slaves())skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=64L),z=sin(seq_len(64L)))
  y <- cos(seq_len(64L)*.37)
  ex <- x[c(1L,5L,40L,64L),,drop=FALSE]
  for(h in c(.05,1e6)) {
    bw <- npregbw(xdat=x,ydat=y,bws=rep(h,2L),regtype="lp",
      degree=c(3L,3L),bernstein.basis=TRUE,bandwidth.compute=FALSE)
    fit <- .np_regression_direct(bw,x,y,ex,se=TRUE,gradients=TRUE,
      gradient.order=3L)
    H <- npreghat(bw,txdat=x,exdat=ex,s=c(0L,3L))
    expect_equal(drop(H%*%y),fit$grad[,2L],tolerance=2e-11)
    factor <- .npRmpi_with_local_regression(.np_plot_lsq_derivative_factor(bw,x,y,ex,coordinate=2L,
      derivative.order=3L,estimate=fit$grad[,2L],se=fit$gerr[,2L]))
    expect_equal(factor$estimate,fit$grad[,2L],tolerance=2e-11)
    expect_equal(factor$A,sweep(H,2L,factor$residual,"*"),
      tolerance=2e-11,ignore_attr=TRUE)
    expect_error(npreghat(bw,txdat=x,exdat=ex,s=c(0L,4L)),
      "exceeds local polynomial degree",fixed=TRUE)
  }
})


test_that("public mean capture preserves its owner and all-large ridge", {
  old <- options(np.messages=FALSE,np.tree=TRUE)
  on.exit(options(old),add=TRUE)
  if(!spawn_mpi_slaves())skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  n <- 64L; i <- seq_len(n)
  x <- data.frame(x=seq(-1,1,length.out=n))
  x$z <- x$x+1e-9*sin(i)
  y <- cos(i*.37)+.1*sin(i*sqrt(2))
  expect_warning(bw <- npregbw(xdat=x,ydat=y,bws=rep(1e6,2),regtype="lp",degree=c(1L,1L),
    bernstein.basis=TRUE,bandwidth.compute=FALSE), "regression design matrix is ill-conditioned", fixed=TRUE)
  original <- getFromNamespace(".np_plot_lsq_gradient_projection","npRmpi")
  mean.requests <- 0L
  local_mocked_bindings(.np_plot_lsq_gradient_projection=function(...,capture.mean.only=FALSE) {
    expect_true(capture.mean.only)
    mean.requests <<- mean.requests+1L
    original(...,capture.mean.only=capture.mean.only)
  },.package="npRmpi")
  expect_warning(H <- npreghat(bw,txdat=x,exdat=x), "regression design matrix is ill-conditioned", fixed=TRUE)
  expect_identical(mean.requests,1L)
  expect_warning(fit <- npreg(bw,txdat=x,tydat=y,exdat=x,se=FALSE), "regression design matrix is ill-conditioned", fixed=TRUE)
  expect_true(any(attr(H,"ridge.used")>0))
  expect_equal(drop(H%*%y),fit$mean,tolerance=2e-11)
  for(tree in list(FALSE,"auto")) {
    options(np.tree=tree)
    expect_warning(E <- npreghat(bw,txdat=x,exdat=x), "regression design matrix is ill-conditioned", fixed=TRUE)
    expect_equal(E,H,tolerance=0)
  }
})
