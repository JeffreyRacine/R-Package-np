test_that("asymmetric Q6 has an explicit selected-coordinate contract", {
  a <- list(errors="bootstrap",bootstrap="smoothed-q6",center="bias-corrected",
    band="simultaneous",B=199L,boot.control=np_boot_control(calibration="asymmetric-sd"))
  expect_true(.np_udens_q6_options(a,TRUE))
  expect_error(.np_udens_q6_options(a),"only.*density")
  for (band in c("pointwise","bonferroni","all","pmzsd")) {
    b<-a;b$band<-band;expect_error(.np_udens_q6_options(b,TRUE),"requires explicit")
  }
  for (design in list(NULL,list(vary=1.1,xtrim=0),list(vary=NA,xtrim=0),
      list(varying=1,xtrim=0),list(vary=1,xtrim=.1),list(vary=1,xtrim=-.1),
      list(vary=1,xtrim=list(c(.1,.9)))))
    expect_error(.np_udens_q6_asymmetric_design(design),"requires")
  expect_null(.np_udens_q6_asymmetric_design(list(vary=1,xtrim=0)))
  expect_null(.np_udens_q6_asymmetric_design(list(vary=1,xtrim=list(c(0,1)))))
  expect_error(.np_udens_q6_asymmetric_budget(1000,10000,999),"workspace bound")
  expect_gt(.np_udens_q6_asymmetric_budget(20,3,17),.np_udens_q6_budget(20,3,17))
})

test_that("asymmetric Q6 signed reduction retains negative critical values", {
  r<-list(index=1:2,center=c(1,2),se=c(.1,.2))
  z<-.np_udens_q6_asymmetric_intervals(r,seq(-3,-1,length.out=99),seq(3,5,length.out=99),.1)
  expect_lt(z$critical[["lower"]],0)
  expect_gt(z$lower[1],r$center[1])
  expect_identical(z$lower,r$center-unname(quantile(seq(-3,-1,length.out=99),.95))*r$se)
  mirror<-r;mirror$center<--r$center
  q<-.np_udens_q6_asymmetric_intervals(mirror,seq(3,5,length.out=99),seq(-3,-1,length.out=99),.1)
  expect_identical(q$lower,-z$upper);expect_identical(q$upper,-z$lower)
  expect_error(.np_udens_q6_asymmetric_intervals(r,NA_real_,1,.05),"calibration")
})

test_that("selected Q6 roots match literal complete-sample roots", {
  x<-qnorm((seq_len(24)-.5)/24);h<-.6;grid<-seq(min(x),max(x),length.out=5)
  design<-list(vary=1L,xtrim=0)
  for(kernel in c("uniform","epanechnikov")) {
    b<-structure(list(bw=h,scaling=FALSE,type="fixed",ckerbound="none",
      ckertype=kernel,ckerorder=2L,ndim=1L,ncon=1L),class="bandwidth")
    state<-.np_udens_q6_asymmetric_prepare(b,data.frame(x=x),42L,design)
    old<-.np_udens_q6_prepare(b,data.frame(x=x),42L)
    expected<-old;expected$plot.calibration<-"asymmetric-sd";expected$construction<-state$construction
    expect_identical(state,expected);expect_false(.np_udens_q6_active(state))
    expect_true(.np_udens_q6_asymmetric_active(state))
    family<-.np_udens_q6_asymmetric_family(state,data.frame(x=grid),design)
    expect_error(.np_udens_q6_asymmetric_family(state,data.frame(x=grid+.01),design),"does not match")
    r<-.np_udens_q6_reference(x,h,grid,kernel);r$fractions<-family$fractions
    tasks<-.np_with_seed(42L,.np_udens_q6_plan(.Random.seed,length(x),17L))
    Y<-do.call(cbind,lapply(tasks,.np_udens_q6_samples,x=x,h=h))
    literal<-matrix(0,ncol(Y),length(grid))
    for(i in seq_len(ncol(Y)))for(j in seq_along(grid)) {
      target<-min(Y[,i])+r$fractions[j]*diff(range(Y[,i]))
      density<-mean(dt((target-x)/(h*sqrt(2/3)),6))/(h*sqrt(2/3))
      literal[i,j]<-(.np_udens_q6_curve(Y[,i,drop=FALSE],target,h,kernel,64L)-density)/r$se[j]
    }
    z<-.np_udens_q6_asymmetric_maximum(Y,r)
    expect_equal(z$plus,apply(literal,1,max),tolerance=1e-11)
    expect_equal(z$minus,apply(-literal,1,max),tolerance=1e-11)
    fixed<-.np_udens_q6_maximum(Y,r)$maximum
    expect_false(isTRUE(all.equal(z$plus,fixed,tolerance=1e-8)))
    split<-do.call(rbind,lapply(tasks,.np_udens_q6_asymmetric_chunk,reference=r))
    expect_equal(split[,1],z$plus,tolerance=1e-12);expect_equal(split[,2],z$minus,tolerance=1e-12)
    perm<-c(5,1,5,3,2,4)
    expect_identical(.np_udens_q6_asymmetric_family(state,data.frame(x=grid[perm]),design),family)
    # Target construction is affine-equivariant, without changing the fitted law.
    rr<-r;rr$x<-3+2*x;rr$h<-2*h;rr$grid<-3+2*grid;rr$se<-r$se/2
    zz<-.np_udens_q6_asymmetric_maximum(3+2*Y,rr)
    expect_equal(zz$plus,z$plus,tolerance=1e-10);expect_equal(zz$minus,z$minus,tolerance=1e-10)
  }
})
