test_that("ill-conditioned accepted rows reproduce affine derivatives and their covariance", {
  old <- options(np.messages=FALSE,np.tree=FALSE,matprod="default")
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=64L),z=sin(seq_len(64L)))
  y <- x$z
  e <- x[c(5L,20L,40L,60L),,drop=FALSE]
  bw <- npregbw(xdat=x,ydat=y,bws=c(.03,.03),regtype="lp",
    degree=c(1L,1L),bernstein.basis=TRUE,bandwidth.compute=FALSE)
  fit <- .np_regression_direct(bw,x,y,e,se=TRUE,gradients=TRUE,fit.inference=TRUE)
  plain <- .np_regression_direct(bw,x,y,e,se=FALSE,gradients=TRUE)
  expect_identical(fit$mean,plain$mean)
  expect_identical(fit$grad,plain$grad)
  unregularized <- fit$fit.inference[[3L]]==0
  expect_true(any(unregularized))
  expect_equal(fit$mean[unregularized],e$z[unregularized],tolerance=1e-10)
  expect_lt(max(abs(fit$grad[unregularized,1L])),1e-10)
  expect_lt(max(abs(fit$grad[unregularized,2L]-1)),1e-10)
  policy <- attr(fit$fit.projection[[2L]],".np.fit.precision",exact=TRUE)
  expect_true(any(policy[,2L]>0))
  for(j in 0:2) {
    der <- if(j==0L)NULL else as.integer(seq_len(2L)==j)
    op <- .np_plot_regression_capture_operator(fit,bw,x,y,e,der,NULL,64*1024^2)
    H <- op$operator(e,der,seq_len(nrow(e)))
    tiled <- do.call(rbind,lapply(seq_len(nrow(e)),function(i)
      op$operator(e[i,,drop=FALSE],der,i)))
    expect_identical(H,tiled)
    target <- if(j==0L)fit$mean else fit$grad[,j]
    expect_equal(drop(H%*%y),target,tolerance=1e-10)
    C <- tcrossprod(sweep(H,2L,op$residual,"*"))
    got <- .np_plot_regression_joint_covariance(bw,x,y,e,fit=fit,deriv=der)
    expect_equal(got$covariance,C,tolerance=1e-12)
    expect_equal(sqrt(diag(C)),if(j==0L)fit$merr else fit$gerr[,j],tolerance=1e-12)
    for(mode in c("default","internal")) {
      options(matprod=mode)
      hat <- npreghat(bw,txdat=x,exdat=e,deriv=if(j==0L)integer(2L)else der)
      expect_equal(drop(hat%*%y),target,tolerance=1e-9)
    }
  }
})

test_that("precision replay policy follows reordered and singleton target slices", {
  x <- matrix(seq_len(15L)*.1,5L,3L)
  attr(x,".np.fit.precision") <- cbind(c(0,0,.1,0,.2),seq_len(5L))
  ix <- c(5L,2L,4L)
  y <- .np_plot_projection_slice(x,ix,c(3L,1L))
  expect_identical(unname(y),unname(structure(x[ix,c(3L,1L),drop=FALSE],
    .np.fit.precision=attr(x,".np.fit.precision")[ix,,drop=FALSE])))
  expect_identical(dim(attr(.np_plot_projection_slice(x,2L),".np.fit.precision")),c(1L,2L))
  attr(x,".np.fit.precision") <- matrix(0,4L,2L)
  expect_error(.np_plot_projection_slice(x,1L),"invalid accepted precision")
})

test_that("bounded Gaussian precision preserves accepted ridge scale and capture", {
  old<-options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(-1,1,length.out=64L),z=sin(seq_len(64L)),w=cos(seq_len(64L)*.71))
  y<-x$w+.02*sin(seq_len(64L)*.5);e<-x[c(1L,5L,40L,64L),,drop=FALSE]
  for(bound in c("fixed","range","one-sided","partial")) {
    args<-list(xdat=x,ydat=y,bws=rep(.05,3L),regtype="lp",degree=rep(2L,3L),
      bernstein.basis=TRUE,ckerbound=if(bound=="range")"range"else"fixed",bandwidth.compute=FALSE)
    if(bound!="range") {
      lower<-vapply(x,min,0);upper<-vapply(x,max,0)
      if(bound=="one-sided")upper[]<-Inf
      if(bound=="partial"){lower[2:3]<--Inf;upper[2:3]<-Inf}
      args<-c(args,list(ckerlb=lower,ckerub=upper))
    }
    bw<-do.call(npregbw,args)
    fit<-.np_regression_direct(bw,x,y,e,se=TRUE,gradients=TRUE,fit.inference=TRUE)
    plain<-.np_regression_direct(bw,x,y,e,se=FALSE,gradients=TRUE)
    expect_identical(fit$mean,plain$mean);expect_identical(fit$grad,plain$grad)
    policy<-attr(fit$fit.projection[[2L]],".np.fit.precision",exact=TRUE)
    expect_true(any(policy[,1L]>0 & policy[,2L]>0))
    expect_true(.np_regression_precision_scope(bw))
    expect_identical(.np_plot_regression_capture_supported(bw),bound!="range")
    for(j in 0:3) {
      deriv<-if(j)as.integer(seq_len(3L)==j)else integer(3L)
      hat<-npreghat(bw,txdat=x,exdat=e,deriv=deriv)
      expect_equal(drop(hat%*%y),if(j)fit$grad[,j]else fit$mean,tolerance=1e-9)
    }
    if(bound!="range")for(j in 0:3) {
      deriv<-if(j)as.integer(seq_len(3L)==j)else NULL
      op<-.np_plot_regression_capture_operator(fit,bw,x,y,e,deriv,NULL,64*1024^2)
      H<-op$operator(e,deriv,seq_len(nrow(e)))
      expect_equal(drop(H%*%y),if(j)fit$grad[,j]else fit$mean,tolerance=1e-10)
      expect_equal(sqrt(rowSums(sweep(H,2L,op$residual,"*")^2)),
        if(j)fit$gerr[,j]else fit$merr,tolerance=1e-10)
    }
    expect_identical(.np_bias_regression_kernel_class(bw,NULL),
      if(bound!="range")"incumbent"else"unsupported")
  }
})


test_that("bounded Gaussian correction keeps its continuous fixed-bandwidth scope", {
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(-1,1,length.out=64L));y<-sin(x$x)
  for(type in c("generalized_nn","adaptive_nn")) {
    bw<-npregbw(xdat=x,ydat=y,bws=10,regtype="lp",degree=2L,
      bwtype=type,ckerbound="fixed",ckerlb=-1,ckerub=1,bandwidth.compute=FALSE)
    expect_false(.np_regression_precision_scope(bw))
    expect_false(.np_plot_regression_capture_supported(bw))
    expect_identical(.np_bias_regression_kernel_class(bw,NULL),"unsupported")
    expect_error(.np_regression_direct(bw,x,y,x[1:2,,drop=FALSE],
      se=TRUE,fit.inference=TRUE),"compact regression inference requires")
  }
  for(ordered in c(FALSE,TRUE)) {
  x$g<-factor(rep(c("a","b"),32L),ordered=ordered)
  bw<-npregbw(xdat=x,ydat=y,bws=c(.3,.2),regtype="lp",degree=2L,
    ckerbound="fixed",ckerlb=-1,ckerub=1,bandwidth.compute=FALSE)
  expect_false(.np_regression_precision_scope(bw))
  expect_false(.np_plot_regression_capture_supported(bw))
  expect_identical(.np_bias_regression_kernel_class(bw,NULL),"unsupported")
  expect_error(.np_bias_regression_prepare(bw,x,y),"continuous predictors only")
  expect_error(.np_regression_direct(bw,x,y,x[1:2,,drop=FALSE],
    se=TRUE,fit.inference=TRUE),"bounded kernels require continuous predictors only")
  }
})


test_that("bounded Gaussian preparation preserves support and active adaptation", {
  old<-options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(-1,1,length.out=64L));y<-sin(3*x$x)+.1*cos(seq_len(64L))
  make<-function(h)npregbw(xdat=x,ydat=y,bws=h,regtype="ll",
    ckerbound="fixed",ckerlb=-1,ckerub=1,bandwidth.compute=FALSE)
  bw<-make(.2);held<-bw
  st<-.np_bias_regression_prepare(bw,x,y,adapt=FALSE,gradient.order=1L,
    law="inid",calibration="basic-local-envelope")
  expect_identical(st$bws[["ckerbound"]],"fixed")
  expect_equal(st$bws[["ckerlb"]],(-1-st$offset)/st$scale,ignore_attr=TRUE)
  expect_equal(st$bws[["ckerub"]],(1-st$offset)/st$scale,ignore_attr=TRUE)
  grid<-x[c(8L,32L,56L),,drop=FALSE]
  fit<-.np_regression_direct(bw,x,y,grid,se=TRUE,gradients=TRUE,fit.inference=TRUE)
  capture<-.np_bias_capture_prepare(st,fit,grid,c("mean","d1"))
  expect_s3_class(capture,"np_bias_capture_prototype")
  expect_identical(capture$state,st)
  expect_identical(bw,held)
  active<-.np_bias_regression_prepare(make(3),x,y,adapt=TRUE,
    gradient.order=1L,law="inid",calibration="basic-local-envelope")
  expect_true(isTRUE(active$settings$adaptation$active))
  expect_true(is.list(active$adaptation))
})

test_that("bounded Gaussian hat precision reproduces ill-conditioned affine derivatives", {
  old<-options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(-1,1,length.out=64L),z=sin(seq_len(64L)))
  y<-x$z;e<-x[c(5L,20L,40L,60L),,drop=FALSE]
  bw<-npregbw(xdat=x,ydat=y,bws=c(.03,.03),regtype="lp",
    degree=c(1L,1L),bernstein.basis=TRUE,ckerbound="fixed",
    ckerlb=vapply(x,min,0),ckerub=vapply(x,max,0),bandwidth.compute=FALSE)
  fit<-.np_regression_direct(bw,x,y,e,se=TRUE,gradients=TRUE,fit.inference=TRUE)
  unregularized<-fit$fit.inference[[3L]]==0
  expect_true(any(unregularized))
  expect_equal(fit$mean[unregularized],e$z[unregularized],tolerance=1e-10)
  expect_lt(max(abs(fit$grad[unregularized,1L])),1e-10)
  expect_lt(max(abs(fit$grad[unregularized,2L]-1)),1e-10)
  for(j in 0:2) {
    deriv<-if(j)as.integer(seq_len(2L)==j)else integer(2L)
    H<-npreghat(bw,txdat=x,exdat=e,deriv=deriv)
    expect_equal(drop(H%*%y),if(j)fit$grad[,j]else fit$mean,tolerance=1e-9)
  }
})
