test_that("plot metadata distinguishes protective defaults from explicit calibrations", {
  request<-list(errors="bootstrap",bootstrap="wild",center="bias-corrected",alpha=.05,gradients=FALSE)
  for(effective in c("absolute","basic-local-envelope","tail-envelope"))for(requested in c("default",effective)) {
    raw<-list(inference=list(calibration.requested=requested,calibration.method=effective,law="rademacher"))
    value<-.np_plot_data_metadata(raw,request)
    expect_identical(value$calibration,effective)
    expected<-if(requested=="default"&&effective%in%c("basic-local-envelope","tail-envelope"))"default"else NULL
    expect_identical(value$calibration_details$requested,expected)
  }
})

test_that("public Rademacher defaults retain requested metadata and identical intervals", {
  skip_on_cran()

  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(-.9,.9,length.out=60L));y<-sin(3*x$x)+.15*cos(seq_len(60L))
  for(kernel in c("epanechnikov","beta","uniform","gaussian")) {
    args<-list(xdat=x,ydat=y,bws=if(kernel=="beta").3 else .5,regtype=if(kernel=="uniform")"lc"else"ll",ckertype=kernel,bandwidth.compute=FALSE)
    if(kernel=="beta")args<-c(args,list(ckerbound="fixed",ckerlb=-1,ckerub=1));bw<-do.call(npregbw,args);held.bw<-bw
    set.seed(112871L);held<-.Random.seed
    run<-function(control,interval=NULL)capture_plot_provider(plot(bw,xdat=x,ydat=y,errors="bootstrap",center="bias-corrected",band="all",B=199L,neval=3L,alpha=.1,output="data",boot.control=control,interval=interval))
    automatic<-run(np_boot_control(wild="rademacher"));effective<-if(kernel=="gaussian")"absolute"else"tail-envelope";explicit<-run(np_boot_control(wild="rademacher"),interval=regression_test_interval(effective))
    a<-automatic$provider[[1L]]$inference;e<-explicit$provider[[1L]]$inference
    expect_identical(a$calibration.requested,"default");expect_identical(a$calibration.method,effective)
    a$calibration.requested<-e$calibration.requested;expect_identical(a,e)
    expect_identical(automatic$result$panels,explicit$result$panels)
    for(meta in automatic$result$inference){expect_identical(meta$components$simultaneous$calibration,effective);expect_identical(meta$resampling_law,"rademacher");expect_identical(meta$calibration_details$requested,if(kernel=="gaussian")NULL else"default")}
    for(meta in explicit$result$inference)expect_null(meta$calibration_details$requested)
    expect_identical(.Random.seed,held);expect_identical(bw,held.bw)
  }
})
