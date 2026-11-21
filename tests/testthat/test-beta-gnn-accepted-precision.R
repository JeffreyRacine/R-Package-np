test_that("beta GNN fitted rows preserve output requests and fitted-policy hats", {
  if(!spawn_mpi_slaves())skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  i<-seq_len(64L);x<-data.frame(x=(i-.5)/64,z=(sin(i*sqrt(2))+1)/2)
  x[5L,]<-x[4L,];y<-cos(i*.37)+.1*x$z
  bw<-npregbw(xdat=x,ydat=y,bws=c(6,6),bwtype="generalized_nn",ckertype="beta",
    ckerbound="fixed",ckerlb=c(0,0),ckerub=c(1,1),regtype="lp",degree=c(2L,2L),
    bernstein.basis=TRUE,bandwidth.compute=FALSE)
  for(e in list(NULL,x,x[c(5L,20L,40L,60L),,drop=FALSE])) {
    evaluation<-if(is.null(e))list()else list(exdat=e)
    fitrun<-function(...)do.call(npreg,c(list(bws=bw,txdat=x,tydat=y),evaluation,list(...)))
    hatrun<-function(...)do.call(npreghat,c(list(bws=bw,txdat=x),evaluation,list(...)))
    fit<-fitrun(gradients=TRUE,se=FALSE)
    sefit<-fitrun(gradients=TRUE,se=TRUE)
    meanfit<-fitrun(gradients=FALSE,se=FALSE)
    expect_identical(fit$mean,sefit$mean);expect_identical(fit$grad,sefit$grad)
    expect_identical(fit$mean,meanfit$mean)
    for(j in 0:2) {
      s<-as.integer(seq_len(2L)==j)
      H<-hatrun(s=s)
      expect_equal(drop(H%*%y),if(j)fit$grad[,j]else fit$mean,tolerance=2e-11)
      yy<-cbind(y,rev(y));app<-hatrun(s=s,y=yy,output="apply")
      expect_equal(unname(app),unname(H%*%yy),tolerance=2e-11)
    }
  }
})

test_that("beta GNN certificate preserves slices and rejects invalid fitted policy", {
  if(!spawn_mpi_slaves())skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  x<-data.frame(x=(seq_len(64L)-.5)/64);y<-sin(x$x*3)
  bw<-npregbw(xdat=x,ydat=y,bws=6,bwtype="generalized_nn",ckertype="beta",
    ckerbound="fixed",ckerlb=0,ckerub=1,regtype="lp",degree=1L,
    bernstein.basis=TRUE,bandwidth.compute=FALSE)
  e<-x[c(5L,25L,55L),,drop=FALSE]
  localrun<-function(expr).npRmpi_with_local_regression(expr)
  cap<-localrun(.np_plot_lsq_gradient_projection(bw,x,y,e,capture.gradient=TRUE))
  P<-cap$fit.projection[[2]];B<-cap$fit.projection[[1]]
  for(order in 0:1) {
    full<-localrun(.np_plot_lsq_beta_moment_weights(bw,x,y,B,e,accepted.projection=P,derivative=c(1L,order)))
    part<-localrun(.np_plot_lsq_beta_moment_weights(bw,x,y,B,e,index=c(3L,1L),accepted.projection=P,derivative=c(1L,order)))
    expect_identical(part$derivative.weights,full$derivative.weights[,c(3L,1L),drop=FALSE])
    expect_identical(part$operator.error,full$operator.error[,c(3L,1L),drop=FALSE])
    expect_true(all(full$certified));expect_true(all(full$operator.relative.error<=2e-11))
    expect_equal(drop(crossprod(full$derivative.weights,y)),if(order)cap$gradient[,1L]else cap$mean,tolerance=2e-11)
  }
  bad<-P;attr(bad,".np.fit.precision")<-matrix(0,2,2)
  expect_error(localrun(.np_plot_lsq_beta_moment_weights(bw,x,y,B,e,accepted.projection=bad,derivative=c(1L,0L))),"invalid accepted precision")
  bad<-P;bad[1L,1L]<-bad[1L,1L]+.01
  expect_error(localrun(.np_plot_lsq_beta_moment_weights(bw,x,y,B,e,accepted.projection=bad,derivative=c(1L,0L))),"level adjoint differs")
})

test_that("beta GNN derivative factors preserve the numerical certificate branch", {
  if(!spawn_mpi_slaves())skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.plot.progress=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=(seq_len(64L)-.5)/64)
  y <- sin(x$x*3)
  bw <- npregbw(xdat=x,ydat=y,bws=6,bwtype="generalized_nn",ckertype="beta",
    ckerbound="fixed",ckerlb=0,ckerub=1,regtype="lp",degree=1L,
    bernstein.basis=TRUE,bandwidth.compute=FALSE)
  e <- x[c(5L,25L,55L),,drop=FALSE]
  fit <- npreg(bws=bw,txdat=x,tydat=y,exdat=e,gradients=TRUE,se=TRUE)
  original <- .np_plot_lsq_beta_moment_weights
  seen <- list()
  observe <- function(...) {
    args <- list(...)
    ans <- original(...)
    direction <- args[["derivative",exact=TRUE]]
    if(!is.null(direction) && direction[2L]>0L)
      seen[[length(seen)+1L]] <<- list(error=ans[["operator.error",exact=TRUE]],
        certified=ans[["certified",exact=TRUE]])
    ans
  }
  local_mocked_bindings(.np_plot_lsq_beta_moment_weights=observe,.package="npRmpi")
  localrun <- function(expr).npRmpi_with_local_regression(expr)
  # Observe the actual row certificate, independently of whether the unchanged
  # subsequent factor certificate can establish sufficient precision today.
  value <- tryCatch(localrun(.np_plot_lsq_derivative_factor(bw,x,y,e,coordinate=1L,
    estimate=fit$grad[,1L],se=fit$gerr[,1L])),np_plot_joint_unavailable=identity)
  expect_length(seen,1L)
  expect_identical(dim(seen[[1L]]$error),c(nrow(x),nrow(e)))
  expect_true(all(is.finite(seen[[1L]]$error)))
  expect_true(all(seen[[1L]]$error>=0))
  expect_true(all(seen[[1L]]$certified))
})
