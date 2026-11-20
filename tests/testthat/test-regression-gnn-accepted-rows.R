test_that("GNN fitted precision follows the saved native policy and original design", {
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  z<-readRDS(test_path('fixtures','gnn-accepted-policy.rds'))
  bw<-npregbw(xdat=z$x,ydat=z$y,bws=4,regtype='lp',degree=2L,
    bernstein.basis=TRUE,bwtype='generalized_nn',ckertype='epanechnikov',bandwidth.compute=FALSE)
  capture<-(.np_plot_lsq_gradient_projection(bw,z$x,z$y,z$ex,capture.gradient=TRUE))
  expect_identical(capture$ridge,z$native.ridge)
  expect_identical(z$native.ridge,0)
  B<-capture$fit.projection[[1]]
  K<-(.np_plot_lsq_gate_weights(bw,.np_plot_lsq_gate_geometry(bw,z$x,z$ex)))
  w<-drop(K);A<-B*sqrt(w/max(w));Y<-z$y*sqrt(w/max(w))
  scales<-sqrt(colSums(A^2));coef<-qr.coef(qr(sweep(A,2,scales,'/'),LAPACK=TRUE),Y)/scales
  fit<-npreg(bw,txdat=z$x,tydat=z$y,exdat=z$ex,gradients=TRUE,se=FALSE)
  for(order in 0:1){V<-W.lp(z$x,exdat=z$ex,degree=2L,gradient.vec=if(order)1L else NULL,basis='glp',bernstein.basis=TRUE)
    target<-if(order)fit$grad[,1]else fit$mean
    expect_equal(target,drop(V%*%coef),tolerance=2e-11)
    H<-npreghat(bw,txdat=z$x,exdat=z$ex,s=order)
    expect_identical(attr(H,'ridge.used'),z$native.ridge)
    expect_equal(drop(H%*%z$y),target,tolerance=2e-11)
    spec<-list(bws=bw,x=z$x,xcon=z$x,numeric.cols=1L,d=1L,gradient.order=1L)
    corrected.main<-.np_bias_main_rows(spec,z$ex,order)
    expect_identical(attr(corrected.main,"ridge.used"),z$native.ridge)
    expect_equal(drop(corrected.main%*%z$y),target,tolerance=2e-11)
  }
  with.se<-npreg(bw,txdat=z$x,tydat=z$y,exdat=z$ex,gradients=TRUE,se=TRUE)
  expect_identical(with.se$mean,fit$mean);expect_identical(with.se$grad,fit$grad)
})

test_that("GNN capture distinguishes training occurrences under tree and dense owners", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  n<-64L;i<-seq_len(n);x<-data.frame(x=seq(-1,1,length.out=n),z=sin(i*sqrt(2)))
  x[5L,]<-x[4L,];x[41L,]<-x[40L,];y<-cos(i*.37)+.1*sin(i*sqrt(2))
  for(kernel in c('gaussian','epanechnikov'))for(tree in c(FALSE,TRUE)){
    options(np.tree=tree)
    bw<-npregbw(xdat=x,ydat=y,bws=c(6,6),bwtype='generalized_nn',ckertype=kernel,
      regtype='lp',degree=c(2L,2L),bernstein.basis=TRUE,bandwidth.compute=FALSE)
    training<-npreg(bw,txdat=x,tydat=y,gradients=TRUE,se=FALSE)
    explicit<-npreg(bw,txdat=x,tydat=y,exdat=x,gradients=TRUE,se=FALSE)
    expect_gt(max(abs(training$mean-explicit$mean)),1e-4)
    for(s in list(c(0L,0L),c(1L,0L))){H<-npreghat(bw,txdat=x,s=s);E<-npreghat(bw,txdat=x,exdat=x,s=s)
      expect_equal(drop(H%*%y),if(any(s))training$grad[,1]else training$mean,tolerance=2e-11)
      expect_equal(drop(E%*%y),if(any(s))explicit$grad[,1]else explicit$mean,tolerance=2e-11)
    }
  }
})


test_that("mixed GNN operators preserve their separate predecessor policy", {
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  local_mocked_bindings(.np_plot_lsq_gradient_projection=function(...)
    stop("mixed operators must retain their existing owner"),.package="np")
  i<-seq_len(64L);x<-data.frame(x=seq(-1,1,length.out=64),z=sin(i*sqrt(2)))
  y<-cos(i*.37);ex<-x[c(4L,17L,40L,61L),]
  for(kernel in c('gaussian','epanechnikov')){
    bw<-npregbw(xdat=x,ydat=y,bws=c(15,15),bwtype='generalized_nn',ckertype=kernel,
      regtype='lp',degree=c(2L,2L),bernstein.basis=TRUE,bandwidth.compute=FALSE)
    expect_no_condition(H<-npreghat(bw,txdat=x,exdat=ex,s=c(1L,1L)))
    expect_true(all(is.finite(H)))
  }
})
