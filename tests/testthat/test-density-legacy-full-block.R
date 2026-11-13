test_that("corrected density legacy rejects a deterministic full circular block", {
  withr::local_options(np.messages=FALSE)
  ns<-asNamespace("np");check<-get(".np_udens_fixed_block_check",environment(get(".np_udens_b2_boot_request",ns)))
  n<-24L
  # Independent law oracle: every possible circular starting point visits each row once.
  for(start in seq_len(n))expect_identical(tabulate(((start-1L+0:(n-1L))%%n)+1L,nbins=n),rep(1L,n))
  expect_error(check("fixed",n,n),"has no resampling variation",fixed=TRUE)
  expect_null(check("fixed",n-1L,n));expect_null(check("geom",n,n));expect_null(check("inid",n,n))
  x<-data.frame(x=seq(-1.3,1.7,length.out=n));z<-data.frame(x=x$x,z=cos(seq_len(n)))
  mixed<-data.frame(x=x$x,g=factor(rep(c("a","b"),n/2L)))
  settings<-list(fixed=list(x=x,b=.4,type="fixed"),gnn=list(x=x,b=8,type="generalized_nn"),ann=list(x=x,b=8,type="adaptive_nn"),q2=list(x=z,b=c(.4,.5),type="fixed"),mixed=list(x=mixed,b=c(.4,.1),type="fixed"))
  for(k in names(settings)) {
    a<-settings[[k]];bw<-npudensbw(dat=a$x,bws=a$b,bwtype=a$type,bandwidth.compute=FALSE);fit<-npudens(bw,tdat=a$x)
    before<-serialize(fit,NULL)
    for(band in c("pointwise","bonferroni","simultaneous","all","pmzsd"))for(selector in 1:2) {
      ctrl<-if(selector==1L)list(boot.control=np_boot_control(calibration="legacy",blocklen=n)) else list(bias.control=np_udens_bias_control("legacy"),boot.control=np_boot_control(blocklen=n))
      set.seed(7301);seed<-.Random.seed
      args<-c(list(x=fit,errors="bootstrap",bootstrap="fixed",center="bias-corrected",B=79L,neval=5L,perspective=FALSE,band=band,output="data"),ctrl)
      expect_error(do.call(plot,args),"fixed block length equal to sample size has no resampling variation",fixed=TRUE,info=paste(k,band,selector))
      expect_identical(.Random.seed,seed);expect_identical(serialize(fit,NULL),before)
    }
  }
  # The incumbent corrected request uses exactly the same condition and message.
  request<-get(".np_udens_b2_boot_request",ns)
  expect_error(request("fixed",n,n),"fixed block length equal to sample size has no resampling variation",fixed=TRUE)
  expect_identical(request("fixed",n-1L,n),list(law="fixed",blocklen=n-1L))
})
