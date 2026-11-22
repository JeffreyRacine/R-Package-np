test_that("basic derivative and contrast roots retain units and per-target families", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(100021);x<-data.frame(x=runif(120),u=factor(rep(1:3,40)))
  y<-sin(3*x$x)+as.numeric(x$u)/5+rnorm(120,sd=.3)
  bw<-npregbw(xdat=x,ydat=y,bws=c(.3,.2),regtype="lp",degree=3L,bandwidth.compute=FALSE)
  grid<-x[seq_len(26L),,drop=FALSE];B<-199L;alpha<-.05
  for(order in 1:2) {
    state<-.np_bias_regression_prepare(bw,x,y,np_bias_control(),gradient.order=order,law="inid")
    targets<-c("mean","d1","c2")
    z<-.np_with_seed(1901,.np_bias_inference(state,grid,targets=targets,B=B,law="inid",
      calibration="basic-envelope",resolution.warn=FALSE))
    expect_match(z$root.units,"endpoint units",fixed=TRUE)
    expect_identical(z$calibration.method,"basic-envelope")
    counts<-.np_bias_stream(z$rng$start,rmultinom(B,nrow(x),rep(1/nrow(x),nrow(x)))-1)
    expect_identical(counts$end,z$rng$end)
    for(target in targets) {
      C<-do.call(rbind,lapply(split(seq_len(nrow(grid)),ceiling(seq_len(nrow(grid))/25)),function(ii)
        .np_bias_regression_rows(state,grid[ii,,drop=FALSE],data.frame(row=seq_along(ii),target=target))$W))
      ii<-which(z$map$target==target);center<-z$center[ii];scale<-z$se[ii]
      expect_equal(rowSums(C),rep(if(target=="mean")1 else 0,nrow(grid)),tolerance=1e-8)
      if(target=="d1") {
        expect_equal(drop(C%*%x$x),rep(if(order==1L)1 else 0,nrow(grid)),tolerance=1e-8)
        expect_equal(drop(C%*%(x$x^2)),if(order==1L)2*grid$x else rep(2,nrow(grid)),tolerance=1e-8)
      }
      F<-sweep(C,2L,state$residual,"*");F<-F-rowMeans(F);R<-F%*%counts$value
      denom<-apply(R,1L,sd);expect_true(all(scale[denom==0]==0));R<-scale*(R/ifelse(denom==0,1,denom))
      rs<-t(apply(R,1L,sort));aa<-t(apply(abs(R),1L,sort));m<-nrow(R)
      k<-1L+floor(B*alpha/2);kb<-1L+floor(B*alpha/(2*m));ks<-ceiling(B*(1-alpha))
      lim<-c(sort(apply(R,2L,min))[k],sort(apply(R,2L,max))[B+1L-k]);aq<-sort(apply(abs(R),2L,max))[ks]
      for(band in c("pointwise","simultaneous","bonferroni")) {
        q<-switch(band,pointwise=rs[,c(k,B+1L-k)],simultaneous=cbind(rep(lim[1],m),rep(lim[2],m)),
          bonferroni=cbind(pmin(rs[,kb],lim[1]),pmax(rs[,B+1L-kb],lim[2])))
        ar<-switch(band,pointwise=aa[,ks],simultaneous=rep(aq,m),bonferroni=pmax(aa[,ceiling(B*(1-alpha/m))],aq))
        expected<-cbind(center-pmax(q[,2],ar),center-pmin(q[,1],-ar));expected[scale==0,]<-center[scale==0]
        expect_equal(unname(z$bands[[band]][ii,,1L]),unname(expected),tolerance=1e-10)
      }
      separate<-.np_with_seed(1901,.np_bias_inference(state,grid,targets=target,B=B,law="inid",calibration="basic-envelope",resolution.warn=FALSE))
      for(band in names(z$bands))expect_equal(unname(z$bands[[band]][ii,,1L]),unname(separate$bands[[band]][,,1L]),tolerance=1e-8)
      if(target=="c2")expect_true(any(scale==0))
    }
  }
})

test_that("public basic derivative plots preserve fitted objects and existing order restrictions", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  set.seed(100022);x<-data.frame(x=runif(100),u=factor(rep(1:2,50)));y<-x$x^3+as.numeric(x$u)+rnorm(100,sd=.3)
  bw<-npregbw(xdat=x,ydat=y,bws=c(.35,.2),regtype="lp",degree=3L,bandwidth.compute=FALSE)
  fit<-npreg(bw,txdat=x,tydat=y);before<-fit;rng<-.Random.seed
  for(law in c("inid","fixed","geom"))for(order in 1:2) {
    out<-withCallingHandlers(plot(fit,errors="bootstrap",center="bias-corrected",bootstrap=law,
      gradients=TRUE,gradient.order=order,band="all",neval=5L,B=199L,output="data",
      interval="basic-envelope", boot.control=np_boot_control(blocklen=if(law=="inid")NULL else 3L)),
      warning=function(w) {
        if(startsWith(conditionMessage(w),"Empirical band quantiles have limited draw resolution"))
          invokeRestart("muffleWarning")
      })
    expect_identical(out$targets$quantity,c("derivative","contrast"))
    expect_identical(out$targets$variable,c("x","u"))
    expect_identical(out$targets$order,c(as.integer(order),NA_integer_))
    for(panel in 1:2) {
      d<-as.data.frame(out,panel=panel,band="simultaneous")
      expect_true(all(is.finite(d$lower)));expect_true(all(d$lower<=d$upper))
      meta<-out$inference[[out$targets$inference_id[panel]]]
      expect_identical(meta$components$simultaneous$calibration,"basic-envelope")
      expect_identical(meta$resampling_law,law)
      expect_identical(meta$scope,"finite grid within each scalar-target family")
      if(panel==2L) {
        coordinate<-which(out$panels[[panel]]$variables$name=="u")
        base<-as.integer(out$panels[[panel]]$grid[[coordinate]])==1L
        expect_true(any(base))
        expect_true(all(d$lower[base]==0 & d$upper[base]==0))
        expect_true(all(out$panels[[panel]]$values$bias_corrected[base]==0))
        expect_true(all(out$panels[[panel]]$values$std_error[base]==0))
      }
    }
    expect_identical(fit,before);expect_identical(.Random.seed,rng)
  }
  expect_error(plot(fit,errors="bootstrap",center="bias-corrected",bootstrap="inid",gradients=TRUE,
    gradient.order=3L,band="simultaneous",neval=5L,B=199L,output="data",
    interval="basic-envelope"),"degree - 1",fixed=TRUE)
  expect_identical(.Random.seed,rng)
})


test_that("basic envelopes retain nearest-neighbour refusals for means and derivatives", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  set.seed(100023);x<-data.frame(x=runif(80));y<-sin(3*x$x)+rnorm(80,sd=.3)
  for(type in c("generalized_nn","adaptive_nn")) {
    bw<-npregbw(xdat=x,ydat=y,bws=25,regtype="ll",bwtype=type,bandwidth.compute=FALSE)
    fit<-npreg(bws=bw,txdat=x,tydat=y);before<-fit;rng<-.Random.seed
    for(gradient in c(FALSE,TRUE)) {
      expect_error(plot(fit,errors="bootstrap",center="bias-corrected",bootstrap="inid",
        gradients=gradient,band="simultaneous",neval=5L,B=199L,output="data",
        interval="basic-envelope"),
        "basic-envelope requires inid, fixed or geom and fixed-bandwidth regression",fixed=TRUE)
      expect_identical(fit,before);expect_identical(.Random.seed,rng)
    }
  }
})
