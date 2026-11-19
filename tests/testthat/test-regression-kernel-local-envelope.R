test_that("kernel envelope admission is explicit, fixed and second order", {
  b <- list(type="fixed",ckertype="epanechnikov",ckerorder=2L,ckerbound="none")
  expect_identical(.np_bias_regression_kernel_class(b,NULL),"unsupported")
  expect_identical(.np_bias_regression_kernel_class(b,"basic-local-envelope"),"local-envelope")
  for(type in c("generalized_nn","adaptive_nn")) {
    q<-b;q[["type"]]<-type
    expect_identical(.np_bias_regression_kernel_class(q,"basic-local-envelope"),"unsupported")
  }
  for(order in c(4L,6L,8L)) {
    q<-b;q[["ckerorder"]]<-order
    expect_identical(.np_bias_regression_kernel_class(q,"basic-local-envelope"),"unsupported")
  }
  for(kernel in c("gaussian","uniform")) {
    q<-b;q[["ckertype"]]<-kernel;q[["ckerbound"]]<-"fixed"
    expect_identical(.np_bias_regression_kernel_class(q,"basic-local-envelope"),"unsupported")
  }
  q<-b;q[["ckerbound"]]<-"range"
  expect_identical(.np_bias_regression_kernel_class(q,"basic-local-envelope"),"unsupported")
  q<-b;q[["ckertype"]]<-"gaussian"
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    q[["type"]]<-type
    expect_identical(.np_bias_regression_kernel_class(q,NULL),"incumbent")
  }
})

test_that("kernel support transport preserves tight bounds and continuous positions", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable in this context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  opts<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(opts),add=TRUE)
  x<-data.frame(u=factor(rep(c("a","b"),40)),x=seq(-3,2,length.out=80))
  y<-sin(x$x)+.15*cos(37*x$x)+as.numeric(x$u)
  bw<-npregbw(xdat=x,ydat=y,bws=c(0,.8),regtype="ll",ckertype="epanechnikov",
    ckerbound="fixed",ckerlb=-3,ckerub=2,bandwidth.compute=FALSE)
  held<-bw;set.seed(104011);seed<-.Random.seed
  state<-.np_bias_regression_prepare(bw,x,y,np_bias_control(),adapt=TRUE,law="mammen",calibration="basic-local-envelope")
  expect_identical(bw,held);expect_identical(.Random.seed,seed)
  expect_equal(state$bws[["ckerlb"]][2],min(state$x$x),tolerance=0)
  expect_equal(state$bws[["ckerub"]][2],max(state$x$x),tolerance=0)
  expect_identical(state$bws[["ckerlb"]][1],bw[["ckerlb"]][1])
  expect_identical(state$pilot.bws[["ckerbound"]],"none")
  expect_equal(unname(state$bws[["bandwidth"]][["x"]][state$numeric.cols]),
    unname(bw[["bandwidth"]][["x"]][state$numeric.cols]/state$scale),tolerance=0)
  check.grid<-x[c(1L,40L,80L),,drop=FALSE]
  original<-.np_bias_regression_rows(state,check.grid,data.frame(row=1:3,target="mean"))$original_center
  expected<-fitted(npreg(bws=bw,txdat=x,tydat=y,exdat=check.grid,se=FALSE))
  expect_equal(original,expected,tolerance=1e-8)
  expect_identical(state$pilot.bws[["ckertype"]],"gaussian")
  expect_identical(state$settings$kernel.admission$lower,bw[["ckerlb"]])
  expect_identical(state$settings$kernel.admission$upper,bw[["ckerub"]])
  expect_identical(state$settings$kernel.admission,list(
    version="regression-kernel-envelope-v1",calibration="basic-local-envelope",
    kernel="epanechnikov",order=2,type="fixed",bound="fixed",
    lower=bw[["ckerlb"]],upper=bw[["ckerub"]]))
  expect_false(isTRUE(state$settings$adaptation$active))
  g<-x[c(1,40,80),,drop=FALSE]
  expect_error(.np_bias_inference(state,g,calibration="absolute",B=19L),"admitted basic-local-envelope")
  expect_error(.np_bias_inference(state,g,calibration="tail-envelope",law="inid",B=19L),"admitted basic-local-envelope")
  expect_error(.np_bias_inference(state,g,calibration="basic-envelope",law="inid",B=19L),"admitted basic-local-envelope")
  expect_error(.np_bias_inference(state,g,calibration="basic-envelope",law="fixed",
    blocklen=5L,refinement="pilot-block",B=19L),"admitted basic-local-envelope")
  expect_error(.np_bias_inference(state,g,calibration="separate-tails",law="inid",B=19L),"admitted basic-local-envelope")
  expect_identical(.Random.seed,seed)
  expect_error(.np_bias_regression_prepare(bw,x,y),"unbounded second-order Gaussian")
})

test_that("public Epanechnikov local envelopes preserve state and target grids", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable in this context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  opts<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(opts),add=TRUE)
  x<-data.frame(x=seq(-3,2,length.out=80));y<-sin(x$x)+.15*cos(37*x$x)
  bw<-npregbw(xdat=x,ydat=y,bws=.8,regtype="ll",ckertype="epanechnikov",bandwidth.compute=FALSE)
  fit<-npreg(bws=bw,txdat=x,tydat=y);held<-fit;set.seed(104012);seed<-.Random.seed
  run<-function(obj,...)suppressWarnings(plot(obj,errors="bootstrap",center="bias-corrected",
    bootstrap="inid",B=99L,neval=5L,band="all",gradients=TRUE,output="data",
    random.seed=104013L,boot.control=np_boot_control(calibration="basic-local-envelope"),...))
  a<-run(fit);b<-run(bw,xdat=x,ydat=y)
  for(band in c("pointwise","simultaneous","bonferroni")) {
    z<-as.data.frame(a,band=band)
    expect_equal(nrow(z),5L);expect_true(all(is.finite(z$lower)));expect_true(all(z$lower<=z$upper))
    expect_identical(z,as.data.frame(b,band=band))
  }
  expect_identical(fit,held);expect_identical(.Random.seed,seed)
  expect_error(plot(fit,errors="bootstrap",center="bias-corrected",bootstrap="inid",
    B=19L,neval=5L,band="all",output="data",boot.control=np_boot_control(calibration="basic-envelope")),"unbounded second-order Gaussian")
  expect_identical(.Random.seed,seed)
})


test_that("pooled kernel states retain their complete admission certificate", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable in this context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  opts<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(opts),add=TRUE)
  x<-data.frame(u=factor(rep(c("a","b","c"),40)),x=seq(-3,2,length.out=120))
  y<-sin(x$x)+.15*cos(37*x$x)+as.numeric(x$u)
  bw<-npregbw(xdat=x,ydat=y,bws=c(.2,.8),regtype="ll",ckertype="epanechnikov",
    ckerbound="fixed",ckerlb=-3,ckerub=2,bandwidth.compute=FALSE)
  held<-bw;set.seed(104014);seed<-.Random.seed
  state<-.np_bias_regression_prepare(bw,x,y,np_bias_control(),adapt=TRUE,law="inid",
    calibration="basic-local-envelope")
  expect_true(isTRUE(state$pooled))
  expect_identical(state$settings$kernel.admission,list(
    version="regression-kernel-envelope-v1",calibration="basic-local-envelope",
    kernel="epanechnikov",order=2,type="fixed",bound="fixed",
    lower=bw[["ckerlb"]],upper=bw[["ckerub"]]))
  expect_identical(state$bw[["ckerlb"]],bw[["ckerlb"]])
  expect_identical(state$bw[["ckerub"]],bw[["ckerub"]])
  before<-state
  expect_error(.np_bias_inference(state,x[c(1,60,120),,drop=FALSE],
    calibration="basic-envelope",law="inid",B=19L),"admitted basic-local-envelope")
  expect_identical(state,before);expect_identical(bw,held);expect_identical(.Random.seed,seed)
})


test_that("beta envelope admission is restricted to fixed second-order known support", {
  b<-list(type="fixed",ckertype="beta",ckerorder=2L,ckerbound="fixed")
  for(cal in list(NULL,"absolute","symmetric-sd","asymmetric-sd","basic-envelope","tail-envelope"))
    expect_identical(.np_bias_regression_kernel_class(b,cal),"unsupported")
  expect_identical(.np_bias_regression_kernel_class(b,"basic-local-envelope"),"local-envelope")
  for(field in c("type","ckerorder","ckerbound")) {
    values<-switch(field,type=c("generalized_nn","adaptive_nn"),ckerorder=c(4L,6L,8L),ckerbound=c("none","range"))
    for(value in values) {q<-b;q[[field]]<-value
      expect_identical(.np_bias_regression_kernel_class(q,"basic-local-envelope"),"unsupported")}
  }
})

test_that("beta support, derivative plots and state ownership survive admission", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable in this context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  opts<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(opts),add=TRUE)
  x<-data.frame(x=seq(-3,2,length.out=80));y<-sin(x$x)+.15*cos(37*x$x)
  args<-list(xdat=x,ydat=y,bws=.3,regtype="lp",degree=2L,ckertype="beta",
    ckerbound="fixed",ckerlb=-3,ckerub=2,bandwidth.compute=FALSE)
  bad<-args;bad$ckerlb<--Inf
  expect_error(do.call(npregbw,bad),"finite")
  bw<-do.call(npregbw,args);held<-bw;set.seed(104012);seed<-.Random.seed
  state<-.np_bias_regression_prepare(bw,x,y,np_bias_control(),adapt=TRUE,
    gradient.order=1L,law="inid",calibration="basic-local-envelope")
  expect_identical(state$bws[["ckertype"]],"beta")
  expect_equal(state$bws[["ckerlb"]],min(state$x$x),tolerance=0)
  expect_equal(state$bws[["ckerub"]],max(state$x$x),tolerance=0)
  expect_equal(unname(state$bws[["bandwidth"]][["x"]][state$numeric.cols]),
    unname(bw[["bandwidth"]][["x"]][state$numeric.cols]/state$scale),tolerance=0)
  check.grid<-x[c(1L,40L,80L),,drop=FALSE]
  original<-.np_bias_regression_rows(state,check.grid,data.frame(row=1:3,target="mean"))$original_center
  expected<-fitted(npreg(bws=bw,txdat=x,tydat=y,exdat=check.grid,se=FALSE))
  expect_equal(original,expected,tolerance=1e-8)
  expect_identical(state$pilot.bws[["ckertype"]],"gaussian")
  expect_identical(state$pilot.bws[["ckerbound"]],"none")
  expect_identical(state$settings$kernel.admission$kernel,"beta")
  expect_false(isTRUE(state$settings$adaptation$active))
  fit<-npreg(bws=bw,txdat=x,tydat=y);heldfit<-fit
  run<-function(obj,output,...)suppressWarnings(plot(obj,errors="bootstrap",center="bias-corrected",
    bootstrap="inid",B=99L,neval=5L,band="all",gradients=TRUE,output=output,
    random.seed=104013L,boot.control=np_boot_control(calibration="basic-local-envelope"),...))
  a<-run(fit,"data");b<-run(bw,"data",xdat=x,ydat=y)
  target<-file.path(tempdir(),"beta-envelope.pdf");grDevices::pdf(target)
  c<-tryCatch(run(fit,"plot-data"),finally=grDevices::dev.off())
  for(band in c("pointwise","simultaneous","bonferroni")) {
    z<-as.data.frame(a,band=band)
    expect_equal(nrow(z),5L);expect_true(all(is.finite(z$lower)))
    expect_true(all(z$lower<=z$upper));expect_identical(z,as.data.frame(b,band=band))
    expect_identical(z,as.data.frame(c,band=band))
  }
  expect_identical(fit,heldfit);expect_identical(bw,held);expect_identical(.Random.seed,seed)
  expect_error(.np_bias_regression_prepare(bw,x,y),"unbounded second-order Gaussian")
  expect_error(.np_bias_inference(state,x[c(1,40,80),,drop=FALSE],calibration="absolute",B=19L),"admitted basic-local-envelope")
  expect_identical(.Random.seed,seed)
})
