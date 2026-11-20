test_that("basic kernel admission is fixed and keeps explicit method certificates", {
  for(kernel in c("epanechnikov","beta")) {
    b <- list(type="fixed",ckertype=kernel,ckerorder=2L,
      ckerbound=if(kernel=="beta")"fixed" else "none",ncon=1L)
    expect_identical(.np_bias_regression_kernel_class(b,"basic-envelope"),"local-envelope")
    expect_identical(.np_bias_regression_calibration(b,"default","bootstrap",
      "simultaneous","inid","mammen"),"basic-local-envelope")
    for(type in c("generalized_nn","adaptive_nn")) {
      z<-b;z$type<-type
      expect_identical(.np_bias_regression_kernel_class(z,"basic-envelope"),"unsupported")
    }
    for(order in c(4L,6L,8L)) {
      z<-b;z$ckerorder<-order
      expect_identical(.np_bias_regression_kernel_class(z,"basic-envelope"),"unsupported")
    }
  }
})

test_that("fixed non-Gaussian basic envelopes match literal pairs roots", {
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(-.95,.95,length.out=60L));y<-sin(3*x$x)+.15*cos(seq_len(60L))
  grid<-x[c(7L,30L,54L),,drop=FALSE];B<-99L;alpha<-.2
  # Literal circular-block covariance, independent of the scale provider.
  variance<-function(F,L,law) {
    F<-F-rowMeans(F);n<-ncol(F)
    if(L==1L)return(rowSums(F^2))
    if(law=="fixed") {
      part<-function(w)if(w==0L)rep(0,nrow(F))else {
        Z<-sapply(seq_len(n),function(i)rowSums(F[,((i-1L+seq_len(w)-1L)%%n)+1L,drop=FALSE]))
        rowMeans(Z^2)
      }
      return((n%/%L)*part(L)+part(n%%L))
    }
    v<-rowSums(F^2)
    for(k in seq_len(n-1L))v<-v+2*(n-k)/n*(1-1/L)^k*
      rowSums(F*F[,((seq_len(n)-1L+k)%%n)+1L,drop=FALSE])
    v
  }
  set.seed(112311L);held<-.Random.seed
  for(kernel in c("epanechnikov","beta")) {
    args<-list(xdat=x,ydat=y,bws=if(kernel=="beta").3 else .5,
      regtype="ll",ckertype=kernel,bandwidth.compute=FALSE)
    if(kernel=="beta")args<-c(args,list(ckerbound="fixed",ckerlb=-1,ckerub=1))
    bw<-do.call(npregbw,args)
    for(law in c("inid","fixed","geom")) {
      st<-.np_bias_regression_prepare(bw,x,y,adapt=TRUE,law=law,calibration="basic-envelope")
      local<-.np_bias_regression_prepare(bw,x,y,adapt=TRUE,law=law,calibration="basic-local-envelope")
      expect_identical(st$settings$kernel.admission$calibration,"basic-envelope")
      same<-st;same$settings$kernel.admission$calibration<-"basic-local-envelope"
      expect_identical(same,local)
      expect_error(.np_bias_inference(st,grid,B=B,law=law,calibration="basic-local-envelope"),
        "admitted basic-envelope calibration")
      L<-if(law=="inid")1L else 5L
      z<-.np_with_seed(112312L,.np_bias_inference(st,grid,B=B,law=law,alpha=alpha,
        blocklen=if(L==1L)NULL else L,calibration="basic-envelope",
        refinement=if(L==1L)"none"else"pilot-block",resolution.warn=FALSE))
      rows<-.np_bias_regression_rows(st,grid,data.frame(row=1:3,target="mean"))
      F<-sweep(rows$W,2L,st$residual,"*");F<-F-rowMeans(F)
      counts<-.np_bias_stream(z$rng$start,if(law=="inid")
        rmultinom(B,nrow(x),rep(1/nrow(x),nrow(x)))-1 else
        .np_block_counts_drawer(nrow(x),B,L,law)(1L,B)-1)
      raw<-F%*%counts$value
      expected.se<-rows$se*sqrt(variance(F,L,law)/variance(F,1L,law))
      root<-expected.se*(raw/apply(raw,1L,sd))
      expect_equal(z$se,expected.se,tolerance=1e-10)
      expect_identical(counts$end,z$rng$end)
      sorted<-t(apply(root,1L,sort));absolute<-t(apply(abs(root),1L,sort))
      k<-1L+floor(B*alpha/2);kb<-1L+floor(B*alpha/(2*nrow(grid)))
      sim<-c(sort(apply(root,2L,min))[k],sort(apply(root,2L,max))[B+1L-k])
      abs.sim<-sort(apply(abs(root),2L,max))[ceiling(B*(1-alpha))]
      for(band in names(z$bands)) {
        limits<-switch(band,pointwise=sorted[,c(k,B+1L-k)],
          simultaneous=matrix(rep(sim,each=nrow(grid)),ncol=2L),
          bonferroni=cbind(pmin(sorted[,kb],sim[1L]),pmax(sorted[,B+1L-kb],sim[2L])))
        q<-switch(band,pointwise=absolute[,ceiling(B*(1-alpha))],
          simultaneous=rep(abs.sim,nrow(grid)),
          bonferroni=pmax(absolute[,ceiling(B*(1-alpha/nrow(grid)))],abs.sim))
        reference<-cbind(rows$center-pmax(limits[,2L],q),rows$center-pmin(limits[,1L],-q))
        expect_equal(unname(z$bands[[band]][,,1L]),unname(reference),tolerance=1e-10)
      }
      expect_error(.np_bias_inference(st,grid,law="mammen",calibration="basic-envelope"),
        "basic-envelope requires inid, fixed or geom")
      expect_identical(.Random.seed,held)
    }
  }
})
