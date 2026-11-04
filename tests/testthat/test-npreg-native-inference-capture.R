local({
  direct <- getFromNamespace('.np_regression_direct','npRmpi')
  covariance <- getFromNamespace('.np_plot_regression_joint_covariance','npRmpi')
  .np_regression_direct <- function(...) .npRmpi_with_local_regression(direct(...))
  .np_plot_regression_joint_covariance <- function(...) .npRmpi_with_local_regression(covariance(...))
test_that("captured global inference matches independent polynomial OLS", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
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
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
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

test_that("accepted covariance replay retains cancellation terms", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=64L),z=sin(seq_len(64L)))
  y <- sin(x$x)+cos(seq_len(64L))*.2
  e <- x[c(5L,20L,40L,60L),,drop=FALSE]
  bw <- npregbw(xdat=x,ydat=y,bws=c(.03,.03),regtype="lp",
    degree=c(1L,1L),bernstein.basis=TRUE,bandwidth.compute=FALSE)
  fit <- .np_regression_direct(bw,x,y,e,se=TRUE,gradients=TRUE,fit.inference=TRUE)
  ordinary <- .np_regression_direct(bw,x,y,e,se=TRUE,gradients=TRUE)
  expect_identical(fit[names(ordinary)],ordinary)
  expect_true(any(fit$fit.inference[[3L]]>0))
  for(j in 0:2) {
    derivative <- if(j==0L)NULL else as.integer(seq_len(2L)==j)
    got <- .np_plot_regression_joint_covariance(bw,x,y,e,fit=fit,deriv=derivative)
    target <- if(j==0L)fit$merr else fit$gerr[,j]
    expect_lt(max(abs(sqrt(diag(got$covariance))-target)),1e-12)
  }
  # This tests covariance ownership. Ill-conditioned forward/adjoint point
  # disagreement is separate; the public plot's existing guard is retained.
})

})

test_that("distributed inference capture retains the public fitted owner", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(n=3L))
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=TRUE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=64L))
  y <- x$x+x$x^2+cos(seq_len(64L))*.2
  e <- data.frame(x=c(-.6,.4)) # Fewer evaluation rows than ranks.
  bw <- npregbw(xdat=x,ydat=y,bws=1000,bandwidth.compute=FALSE,
    regtype="lp",degree=2,bernstein.basis=TRUE)
  fit <- npreg(bw,txdat=x,tydat=y,exdat=e,se=TRUE,gradients=TRUE,
    .np.fit.inference=TRUE,warn.glp.gradient=FALSE)
  ordinary <- npreg(bw,txdat=x,tydat=y,exdat=e,se=TRUE,gradients=TRUE,
    warn.glp.gradient=FALSE)
  for(nm in c("mean","merr","grad","gerr"))
    expect_identical(fit[[nm]],ordinary[[nm]])
  Q <- cbind(1,x$x,x$x^2); q <- qr(Q); U <- qr.Q(q)
  A <- backsolve(qr.R(q),t(U))
  z <- (y-drop(Q%*%qr.coef(q,y)))/sqrt(1-rowSums(U^2))
  for(derivative in c(FALSE,TRUE)) {
    V <- if(derivative)cbind(0,1,2*e$x)else cbind(1,e$x,e$x^2)
    C <- tcrossprod(sweep(V%*%A,2L,z,"*"))
    got <- .np_plot_regression_joint_covariance(bw,x,y,e,fit=fit,
      deriv=if(derivative)1L else NULL)
    expect_equal(got$covariance,C,tolerance=1e-12)
  }
  expect_error(npreg(bw,txdat=x,tydat=y,se=TRUE,.np.fit.inference=TRUE),
    "external fixed")
  recovered <- npreg(bw,txdat=x,tydat=y,exdat=e,se=TRUE)
  expect_identical(recovered$mean,ordinary$mean)
  expect_null(recovered$fit.inference)
})

test_that("distributed large-bandwidth public bands preserve marginal widths", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-2,2,length.out=60L)); y <- x$x+cos(seq_len(60L))
  bw <- npregbw(xdat=x,ydat=y,bws=700,regtype="ll",bandwidth.compute=FALSE)
  for(grad in c(FALSE,TRUE)) {
    a <- plot(bw,errors="asymptotic",band="all",gradients=grad,
      output="data",neval=5L)
    p <- plot(bw,errors="asymptotic",band="pointwise",gradients=grad,
      output="data",neval=5L)
    expect_identical(plot_target(a, band="pointwise")$estimate, plot_target(p)$estimate)
    pointwise <- plot_band(a, "pointwise")
    simultaneous <- plot_band(a, "simultaneous")
    bonferroni <- plot_band(a, "bonferroni")
    expect_identical(pointwise, plot_band(p, "pointwise"))
    expect_true(all(is.finite(simultaneous)))
    expect_true(all(pointwise[,1L]>=simultaneous[,1L] & pointwise[,2L]<=simultaneous[,2L]))
    expect_true(all(simultaneous[,1L]>=bonferroni[,1L] & simultaneous[,2L]<=bonferroni[,2L]))
  }
})

test_that("formula capture provenance survives MPI call materialization", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(42)
  x <- runif(60,-2,2); y <- x+rnorm(60)
  model <- npreg(y~x,regtype="ll")
  bands <- plot(model,errors="asymptotic",band="all",output="data",neval=5L)
  expect_true(all(is.finite(plot_band(bands, "simultaneous"))))
})
