test_that("local envelope extrema match an independent oracle and preserve identities", {
check<-function(...)for(value in list(...))expect_true(value)
set.seed(102004)
B<-199L;alpha<-c(.2,.05,.1);m<-31L
map<-data.frame(row=seq_len(m),target='mean')
scale<-c(0,exp(seq(-3,3,length.out=m-1L)))
Z<-matrix(sample(c(-4,-1,0,0,2,7),m*B,replace=TRUE),m);Z[1,]<-0
R<-Z*scale
oracle<-function(R,scale,alpha) {
 Z<-R/ifelse(scale==0,1,scale);B<-ncol(Z)
 lo<-sort(apply(Z,2,min));hi<-sort(apply(Z,2,max));ab<-sort(apply(abs(Z),2,max))
 sapply(alpha,function(a) {k<-1+floor(B*a/2);ka<-ceiling(B*(1-a));c(min(lo[k],-ab[ka]),max(hi[B+1-k],ab[ka]))})
}
p<-.np_band_local_prepare(B,map,alpha)
for(ii in list(15:31,1:8,9:14))p<-.np_band_local_update(p,R[ii,,drop=FALSE],scale[ii],map[ii,,drop=FALSE],seq_len(B))
f<-.np_band_local_finish(p)
check(identical(unname(f$cutoffs),unname(oracle(R,scale,alpha))))
q<-.np_band_local_prepare(B,map,alpha)
q<-.np_band_local_update(q,R*7,scale*7,map,seq_len(B))
check(isTRUE(all.equal(.np_band_local_finish(q)$cutoffs,f$cutoffs,tolerance=1e-14)))
err<-function(expr,pattern){a<-tryCatch({force(expr);NULL},error=identity);check(inherits(a,'error'),grepl(pattern,conditionMessage(a)))}
err(.np_band_local_finish(.np_band_local_prepare(B,map,alpha)),'missing')
err(.np_band_local_update(p,R,scale,map,seq_len(B)),'duplicate')
err(.np_band_local_update(q,R,scale,map,rev(seq_len(B))),'identities')
bad<-R;bad[1,1]<-1
err(.np_band_local_update(.np_band_local_prepare(B,map,alpha),bad,scale,map,seq_len(B)),'zero target')
err(.np_band_local_prepare(B,map,alpha,byte.limit=1),'budget')
for(sc in c(0,2)) {
 mm<-map[1,,drop=FALSE];z<-.np_band_local_prepare(B,mm,alpha)
 rr<-matrix(if(sc==0)0 else seq(-2,4,length.out=B),1,B)
 z<-.np_band_local_update(z,rr,sc,mm,seq_len(B))
 check(isTRUE(all.equal(unname(.np_band_local_finish(z)$cutoffs),unname(oracle(rr,sc,alpha)),tolerance=0)))
}
# Nine families, inverted endpoints, pointwise identity, alpha order and BF nesting.
fullmap<-expand.grid(row=1:3,target=paste0('d',1:9),stringsAsFactors=FALSE)
ss<-rep(c(0,.1,7),9);local<-vector('list',9)
for(g in 1:9) {ii<-which(fullmap$target==paste0('d',g));m1<-fullmap[ii,,drop=FALSE]
 st<-.np_band_local_prepare(B,m1,alpha);rr<-R[1:3,]*ss[ii]
 st<-.np_band_local_update(st,rr,ss[ii],m1,seq_len(B));local[[g]]<-.np_band_local_finish(st)
}
limits<-setNames(lapply(1:3,function(j){a<-array(0,c(27,2,3));a[,1,]<- -j;a[,2,]<-j;a}),c('pointwise','simultaneous','bonferroni'))
protected<-.np_band_local_protect(limits,local,fullmap,ss)
check(identical(protected$pointwise,limits$pointwise))
for(b in c('simultaneous','bonferroni'))check(all(protected[[b]][,1,]<=limits[[b]][,1,]),all(protected[[b]][,2,]>=limits[[b]][,2,]))
check(all(protected$bonferroni[,1,]<=protected$simultaneous[,1,]),all(protected$bonferroni[,2,]>=protected$simultaneous[,2,]))


})

test_that("local envelope controls preserve defaults and have an explicit name", {
  expect_identical(np_boot_control()$calibration,"default")
  expect_identical(np_boot_control(calibration="basic-local-envelope")$calibration,"basic-local-envelope")
  expect_error(np_boot_control(calibration="basic-local"),"calibration")
})

test_that("public local envelopes enclose both constituents without changing pointwise bounds", {
  old.options<-options(np.messages=FALSE,np.plot.progress=FALSE)
  on.exit(options(old.options),add=TRUE)
  set.seed(102051);x<-data.frame(x=runif(80));y<-sin(3*x$x)+rnorm(80,sd=.4)
  bw<-npregbw(xdat=x,ydat=y,bws=.3,regtype="ll",bandwidth.compute=FALSE)
  fit<-npreg(bws=bw,txdat=x,tydat=y);held<-fit;seed<-.Random.seed
  run<-function(cal,law)withCallingHandlers(plot(fit,errors="bootstrap",center="bias-corrected",
    bootstrap=law,B=199L,neval=6L,band="all",output="data",random.seed=102052L,
    boot.control=np_boot_control(calibration=cal,blocklen=if(law=="inid")NULL else 5L)),
    warning=function(w)if(grepl("resolution",conditionMessage(w),ignore.case=TRUE))invokeRestart("muffleWarning"))
  for(law in c("inid","fixed","geom")) {
    raw<-run("basic-envelope",law);local<-run("basic-local-envelope",law)
    standardized<-run("tail-envelope",law)
    expect_identical(as.data.frame(raw,band="pointwise"),as.data.frame(local,band="pointwise"))
    for(band in c("simultaneous","bonferroni")) {
      a<-as.data.frame(raw,band=band);b<-as.data.frame(local,band=band)
      # The local simultaneous constituent also protects Bonferroni.
      c<-as.data.frame(standardized,band="simultaneous")
      expect_equal(b$lower,pmin(a$lower,c$lower),tolerance=1e-12)
      expect_equal(b$upper,pmax(a$upper,c$upper),tolerance=1e-12)
      expect_true(all(b$lower<=a$lower));expect_true(all(b$upper>=a$upper))
    }
    for(panel in seq_along(local$panels)) {
      i<-local$targets$inference_id[panel]
      expect_identical(local$inference[[i]]$calibration,"basic-local-envelope")
      expect_identical(local$inference[[i]]$resampling_law,law)
    }
    expect_identical(fit,held);expect_identical(.Random.seed,seed)
  }
  expect_error(run("basic-envelope","wild"),"requires inid")
  expect_identical(.Random.seed,seed)
})


test_that("only explicit regression normalization admits the local envelope", {
  dots<-list(errors="bootstrap",center="bias-corrected",bootstrap="inid",
    boot.control=np_boot_control(calibration="basic-local-envelope"))
  for(context in c("plot.lsqregression","plot.npdensity","plot.condensity",
    "plot.condistribution","plot.qregression","plot.conmode",
    "plot.smoothcoefficient","plot.plregression","plot"))
    expect_error(.np_plot_normalize_public_dots(dots,context=context),"requires bias-corrected regression")
  z<-.np_plot_normalize_public_dots(dots,allow_basic_local_envelope=TRUE)
  expect_identical(z$plot.errors.boot.calibration,"basic-local-envelope")
})


test_that("nearest-neighbour local envelopes retain fixed selected indices", {
  opts<-options(np.messages=FALSE,np.plot.progress=FALSE)
  on.exit(options(opts),add=TRUE)
  set.seed(103071);x<-data.frame(x=seq(-1,1,length.out=80));y<-sin(3*x$x)+rnorm(80,sd=.3)
  for(type in c("generalized_nn","adaptive_nn")) {
    bw<-npregbw(xdat=x,ydat=y,bws=45,regtype="ll",bwtype=type,bandwidth.compute=FALSE)
    fit<-npreg(bws=bw,txdat=x,tydat=y);held<-fit;seed<-.Random.seed
    run<-function(cal,object=fit,...)suppressWarnings(plot(object,errors="bootstrap",
      center="bias-corrected",bootstrap="inid",B=199L,neval=6L,band="all",output="data",
      random.seed=103072L,boot.control=np_boot_control(calibration=cal),...))
    z<-run("basic-local-envelope");zz<-run("basic-local-envelope",bw,xdat=x,ydat=y)
    for(b in c("pointwise","simultaneous","bonferroni"))
      expect_identical(as.data.frame(z,band=b),as.data.frame(zz,band=b))
    expect_identical(fit,held);expect_identical(.Random.seed,seed)
    expect_error(run("basic-envelope"),"fixed-bandwidth")
    for(wild in c("mammen","rademacher"))
      expect_error(plot(fit,errors="bootstrap",center="bias-corrected",
        bootstrap="wild",B=199L,neval=6L,band="all",output="data",
        boot.control=np_boot_control(calibration="basic-local-envelope",wild=wild)),
        "fixed-bandwidth")
    expect_identical(.Random.seed,seed)
  }
})
