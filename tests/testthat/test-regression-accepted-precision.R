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
