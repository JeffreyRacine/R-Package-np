test_that("bounded degree-zero NN inference preserves declared support and defaults", {
  skip_on_cran()
  if(!spawn_mpi_slaves(n=1L)) skip("MPI workers unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opt),add=TRUE)
  set.seed(1126201);seed <- .Random.seed
  for(type in c("generalized_nn","adaptive_nn")) for(kind in c("lc","lp0","pooled")) {
    x <- data.frame(x=seq(-.97,.96,length.out=90L));y <- sin(3*x$x)+.1*cos(seq_len(90L))
    if(kind=="pooled"){x$z <- factor(rep(letters[1:3],length.out=90L));y <- y+.2*as.integer(x$z)}
    args <- list(xdat=x,ydat=y,bws=c(25,if(kind=="pooled").15),regtype=if(kind=="lp0")"lp"else"lc",bwtype=type,ckertype="epanechnikov",ckerbound="fixed",ckerlb=-1,ckerub=1,bandwidth.compute=FALSE)
    if(kind=="lp0")args$degree <- 0L
    bw <- do.call(npregbw,args);fit <- npreg(bw,txdat=x,tydat=y);held <- fit
    for(law in c("gaussian","mammen","rademacher","inid","fixed","geom")) {
      st <- .np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,law=law,calibration="absolute")
      expect_identical(st$settings$kernel.admission$type,type)
      expect_identical(st$settings$kernel.admission$bound,"fixed")
      expect_identical(st$settings$kernel.admission$lower,bw$ckerlb)
      expect_identical(st$settings$kernel.admission$upper,bw$ckerub)
      expect_false(isTRUE(st$settings$adaptation$active))
      expect_false(isTRUE(st$settings$mixed$active))
    }
    for(law in list(NULL,NA_character_,c("inid","geom"),"unknown"))
      expect_error(.np_bias_regression_prepare(bw,x,y,law=law,calibration="basic-local-envelope"),"degree-zero kernel NN inference requires")
    run <- function(law="wild",wild="mammen",cal=NULL,errors="bootstrap",band="all") {
      a <- list(x=fit,errors=errors,center="bias-corrected",neval=3L,
        gradients=TRUE,alpha=.2,random.seed=1126202L,output="data",band=band)
      if(errors=="bootstrap"){a$bootstrap <- law;a$B <- 63L}else a$gaussian.draws <- 63L
      if(errors=="bootstrap"&&(!is.null(cal)||wild!="mammen"))a$boot.control <- np_boot_control(wild=wild,calibration=if(is.null(cal))"default"else cal)
      suppressWarnings(capture_plot_provider(do.call(plot,a)))
    }
    for(wild in c("mammen","rademacher")) {
      effective <- if(wild=="mammen")"basic-local-envelope"else"tail-envelope"
      default <- run(wild=wild);explicit <- run(wild=wild,cal=effective)
      a <- default$provider[[1L]]$inference;e <- explicit$provider[[1L]]$inference
      expect_identical(a$calibration.method,effective);expect_identical(a$calibration.requested,"default")
      a$calibration.requested <- e$calibration.requested
      expect_identical(a,e);expect_identical(default$result$panels,explicit$result$panels)
      value <- run(wild=wild,cal="symmetric-sd")
      expect_s3_class(value$result,"np_plot_data")
      expect_identical(value$provider[[1L]]$inference$calibration.method,"absolute")
    }
    expect_s3_class(run(errors="asymptotic")$result,"np_plot_data")
    expect_s3_class(run(band="pmzsd")$result,"np_plot_data")
    expect_s3_class(run(law="inid",cal="absolute")$result,"np_plot_data")
    expect_error(run(cal="basic-envelope"),"basic-envelope requires inid, fixed or geom")
    higher.args <- args;higher.args$regtype <- "lp";higher.args$degree <- 2L
    higher <- do.call(npregbw,higher.args)
    expect_identical(.np_bias_regression_kernel_class(higher,"absolute"),"unsupported")
    expect_error(.np_bias_regression_prepare(higher,x,y,law="mammen",calibration="basic-local-envelope"),"requires a supported second-order kernel")
    expect_identical(.Random.seed,seed);expect_identical(fit,held)
  }
})
