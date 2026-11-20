test_that("regression defaults resolve only qualified fixed kernel laws", {
  b <- list(ckertype="epanechnikov",ckerorder=2L,ckerbound="none",
            type="fixed",ncon=1L,nuno=0L,nord=0L)
  resolve <- function(bws=b, calibration="default", method="bootstrap",
                      band="simultaneous", bootstrap="wild", wild="mammen")
    .np_bias_regression_calibration(bws,calibration,method,band,bootstrap,wild)
  expect_identical(resolve(),"basic-local-envelope")
  for(law in c("inid","fixed","geom"))
    expect_identical(resolve(bootstrap=law),"basic-local-envelope")
  for(kernel in "gaussian") {
    z <- b;z$ckertype <- kernel
    expect_identical(resolve(z),"absolute")
  }
  z <- b;z$ckertype <- "beta";z$ckerbound <- "fixed"
  expect_identical(resolve(z),"basic-local-envelope")
  z$ckerbound <- "none";expect_identical(resolve(z),"absolute")
  for(type in c("generalized_nn","adaptive_nn")) {
    z <- b;z$type <- type;expect_identical(resolve(z),"absolute")
  }
  for(order in c(4L,6L,8L)) {
    z <- b;z$ckerorder <- order;expect_identical(resolve(z),"absolute")
  }
  z <- b;z$ncon <- 0L;expect_identical(resolve(z),"absolute")
  z <- b;z$ckerbound <- "range";expect_identical(resolve(z),"absolute")
  expect_identical(resolve(method="asymptotic"),"absolute")
  expect_identical(resolve(band="pmzsd"),"absolute")
  expect_identical(resolve(wild="rademacher"),"tail-envelope")
  for(cal in c("absolute","separate-tails","tail-envelope","basic-envelope",
               "basic-local-envelope","smoothed-percentile"))
    expect_identical(resolve(calibration=cal),cal)
  expect_identical(resolve(calibration="symmetric-sd"),"absolute")
  expect_identical(resolve(calibration="asymmetric-sd"),"separate-tails")
})

test_that("ordinary fixed Epanechnikov default uses the explicit complete family", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.plot.progress=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=40L));y <- sin(3*x$x)+.2*cos(seq_len(40L))
  b <- npregbw(xdat=x,ydat=y,bws=.4,bandwidth.compute=FALSE,
               regtype="lc",ckertype="epanechnikov")
  set.seed(112L);seed <- .Random.seed;held <- b
  run <- function(control=NULL) {
    args <- list(x=b,xdat=x,ydat=y,errors="bootstrap",center="bias-corrected",
      band="simultaneous",B=63L,neval=3L,output="data")
    if(!is.null(control))args$boot.control <- control
    capture_plot_provider(do.call(plot,args))
  }
  automatic <- run();explicit <- run(np_boot_control(calibration="basic-local-envelope"))
  a <- automatic$provider[[1L]]$inference;e <- explicit$provider[[1L]]$inference
  expect_identical(a$calibration.requested,"default")
  expect_identical(a$calibration.method,"basic-local-envelope")
  a$calibration.requested <- e$calibration.requested
  expect_identical(a,e)
  expect_identical(automatic$result$panels,explicit$result$panels)
  expect_identical(automatic$result$inference[[1L]]$calibration,"basic-local-envelope")
  expect_identical(automatic$result$inference[[1L]]$calibration_details$requested,"default")
  for(control in list(np_boot_control(),np_boot_control(calibration="default")))
    expect_identical(run(control)$provider,automatic$provider)
  expect_identical(.Random.seed,seed);expect_identical(b,held)
  expect_error(run(np_boot_control(calibration="basic-envelope")),"requires inid")
  symmetric <- run(np_boot_control(calibration="symmetric-sd"))
  absolute <- run(np_boot_control(calibration="absolute"))
  s <- symmetric$provider[[1L]]$inference
  a <- absolute$provider[[1L]]$inference
  expect_identical(s$calibration.method,"absolute")
  s$calibration.requested <- a$calibration.requested
  expect_identical(s,a)
  expect_identical(symmetric$result$panels,absolute$result$panels)
})
