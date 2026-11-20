test_that("bounded Gaussian degree-zero NN keeps Gaussian defaults and fitted state", {
  skip_on_cran()
  if(!spawn_mpi_slaves(n=1L))skip("MPI workers unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opt),add=TRUE)
  set.seed(1126801);seed <- .Random.seed
  for(type in c("generalized_nn","adaptive_nn"))for(kind in c("lc","lp0","pooled")) {
    x <- data.frame(x=seq(-.97,.96,length.out=90L));y <- sin(3*x$x)+.1*cos(seq_len(90L))
    if(kind=="pooled"){x$z <- factor(rep(letters[1:3],length.out=90L));y <- y+.2*as.integer(x$z)}
    args <- list(xdat=x,ydat=y,bws=c(25,if(kind=="pooled").15),regtype=if(kind=="lp0")"lp"else"lc",bwtype=type,ckertype="gaussian",ckerbound="fixed",ckerlb=-1,ckerub=1,bandwidth.compute=FALSE)
    if(kind=="lp0")args$degree <- 0L
    bw <- do.call(npregbw,args);fit <- npreg(bw,txdat=x,tydat=y);held <- fit
    expect_true(.np_bias_regression_nn_zero(bw,"gaussian"))
    expect_false(.np_bias_regression_nn_zero(bw))
    for(law in c("gaussian","mammen","rademacher","inid","fixed","geom")) {
      st <- .np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,law=law,calibration="absolute")
      expect_null(st$settings$kernel.admission)
      expect_identical(st$settings$main$kernel,"gaussian")
      expect_false(isTRUE(st$settings$adaptation$active))
      expect_false(isTRUE(st$settings$mixed$active))
    }
    run <- function(law="mammen",cal=NULL,errors="bootstrap",band="all") {
      a <- list(x=fit,errors=errors,center="bias-corrected",neval=3L,
        gradients=TRUE,alpha=.2,random.seed=1126802L,output="data",band=band)
      if(errors=="bootstrap"){a$bootstrap <- if(law%in%c("mammen","rademacher"))"wild"else law;a$B <- 63L}else a$gaussian.draws <- 63L
      if(errors=="bootstrap"&&(!is.null(cal)||law=="rademacher"))a$boot.control <- np_boot_control(wild=if(law=="rademacher")law else"mammen",calibration=if(is.null(cal))"default"else cal)
      suppressWarnings(capture_plot_provider(do.call(plot,a)))
    }
    for(law in c("mammen","rademacher","inid")) {
      default <- run(law);explicit <- run(law,cal="absolute")
      a <- default$provider[[1L]]$inference;e <- explicit$provider[[1L]]$inference
      expect_identical(a$calibration.method,"absolute");expect_identical(a$calibration.requested,"default")
      a$calibration.requested <- e$calibration.requested
      expect_identical(a,e);expect_identical(default$result$panels,explicit$result$panels)
    }
    expect_s3_class(run(errors="asymptotic")$result,"np_plot_data")
    expect_s3_class(run(band="pmzsd")$result,"np_plot_data")
    expect_error(run(cal="basic-envelope"),"basic-envelope requires inid, fixed or geom")
    higher.args <- args;higher.args$regtype <- "lp";higher.args$degree <- 1L
    higher <- do.call(npregbw,higher.args)
    expect_identical(.np_bias_regression_kernel_class(higher,"absolute"),"unsupported")
    expect_error(.np_bias_regression_prepare(higher,x,y,law="mammen",calibration="absolute"),"requires a supported second-order kernel")
    higher$ckerbound <- "none"
    expect_identical(.np_bias_regression_kernel_class(higher,"absolute"),"incumbent")
    expect_identical(.Random.seed,seed);expect_identical(fit,held)
  }
})
