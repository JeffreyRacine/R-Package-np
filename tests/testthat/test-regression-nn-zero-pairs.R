test_that("degree-zero Epanechnikov NN exposes pairs without changing its rows", {
  skip_on_cran()
  if(!spawn_mpi_slaves(n=1L)) skip("MPI workers unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  opts <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opts),add=TRUE)
  x <- data.frame(x=seq(-.97,.96,length.out=120L))
  y <- sin(3*x$x)+.1*cos(seq_len(nrow(x)))
  set.seed(1125601); seed <- .Random.seed
  for(type in c("generalized_nn","adaptive_nn")) for(reg in c("lc","lp")) {
    args <- list(xdat=x,ydat=y,bws=35,regtype=reg,ckertype="epanechnikov",
                 bwtype=type,bandwidth.compute=FALSE)
    if(reg=="lp") args$degree <- 0L
    bw <- do.call(npregbw,args);fit <- npreg(bw,txdat=x,tydat=y);held <- fit
    for(i in 1:3) {
      cal <- c("basic-envelope","separate-tails","tail-envelope")[i]
      law <- c("inid","fixed","geom")[i]
      local <- .np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,
        law=law,calibration="basic-local-envelope")
      st <- .np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,
        law=law,calibration=cal)
      st$settings$kernel.admission$calibration <- "basic-local-envelope"
      expect_identical(st,local)
      value <- suppressWarnings(plot(fit,errors="bootstrap",center="bias-corrected",
        bootstrap=law,B=99L,neval=3L,band="all",gradients=TRUE,alpha=.1,
        random.seed=1125602L,output="data",boot.control=np_boot_control(
          calibration=cal,blocklen=if(law=="inid")NULL else 5L)))
      expect_s3_class(value,"np_plot_data")
      expect_true(all(vapply(value$inference,function(z)
        identical(z$calibration,cal)&&identical(z$resampling_law,law),TRUE)))
      expect_identical(fit,held);expect_identical(.Random.seed,seed)
      if(cal=="basic-envelope")
        expect_error(plot(fit,errors="bootstrap",center="bias-corrected",
          bootstrap="wild",band="all",output="data",boot.control=np_boot_control(calibration=cal)),
          "basic-envelope requires inid, fixed or geom")
    }
    run <- function(control=NULL) {
      a <- list(x=fit,errors="bootstrap",center="bias-corrected",bootstrap="inid",
        B=99L,neval=3L,band="all",gradients=TRUE,alpha=.1,random.seed=1125602L,output="data")
      if(!is.null(control)) a$boot.control <- control
      suppressWarnings(capture_plot_provider(do.call(plot,a)))
    }
    automatic <- run();explicit <- run(np_boot_control(calibration="basic-local-envelope"))
    a <- automatic$provider[[1L]]$inference;e <- explicit$provider[[1L]]$inference
    expect_identical(a$calibration.requested,"default")
    expect_identical(a$calibration.method,"basic-local-envelope")
    a$calibration.requested <- e$calibration.requested;expect_identical(a,e)
    expect_identical(automatic$result$panels,explicit$result$panels)
    expect_identical(.Random.seed,seed);expect_identical(fit,held)
    expect_identical(.np_bias_regression_kernel_class(bw,"absolute"),"local-envelope")
    bounded <- bw;bounded$ckerbound <- "fixed"
    expect_identical(.np_bias_regression_kernel_class(bounded,"tail-envelope"),"local-envelope")
    for(d in 2:3) {
      higher <- bw;higher$regtype <- "lp";higher$basis <- "glp";higher$degree <- d
      expect_identical(.np_bias_regression_kernel_class(higher,"tail-envelope"),"unsupported")
      expect_identical(.np_bias_regression_calibration(higher,"default","bootstrap",
        "simultaneous","inid","mammen"),"absolute")
    }
  }
})
