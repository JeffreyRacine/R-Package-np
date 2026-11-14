test_that("selected asymmetric companions admit only their explicit public choices", {
  a<-list(errors="bootstrap",bootstrap="smoothed-q6",center="bias-corrected",band="all",boot.control=np_boot_control(calibration="asymmetric-sd"))
  for(band in c("pointwise","simultaneous","bonferroni","all")){a$band<-band;expect_true(.np_udens_q6_options(a,TRUE))}
  for(band in c("pmzsd","none",NA_character_)){a$band<-band;expect_error(.np_udens_q6_options(a,TRUE),"requires explicit")}
  for(value in list(c(named="asymmetric-sd"),c(named="symmetric-sd"),factor("asymmetric-sd"),c("asymmetric-sd","asymmetric-sd"),NA_character_,list("asymmetric-sd"))){
    z<-a;z$band<-"all";z$boot.control<-NULL;z$plot.errors.boot.calibration<-value
    expect_error(.np_udens_q6_options(z,TRUE),"requires explicit")
  }
  expect_error(.np_udens_q6_asymmetric_pointwise_block(1000,10000,999),"workspace bound")
  expect_gte(.np_udens_q6_asymmetric_pointwise_block(30,5,39),1)
  for(B in c(1L,19L,39L,40L,999L)){
    z<-.np_udens_resolution(B,5,.05,c("pointwise","simultaneous","bonferroni"),"selected-asymmetric")
    expect_equal(z$nominal_tail_draws,B*c(.025,.025,.005))
    expect_identical(z$resolution_limited,z$nominal_tail_draws<1)
  }
})

test_that("selected marginal roots preserve the signed density target and all components", {
  manual<-function(x,p){x<-sort(x);h<-1+(length(x)-1)*p;i<-floor(h);if(i==length(x))x[i]else x[i]+(h-i)*(x[i+1]-x[i])}
  x<-qnorm((seq_len(24)-.5)/24);h<-.6;g<-seq(min(x),max(x),length.out=5)
  for(kernel in c("uniform","epanechnikov")){
    r<-.np_udens_q6_reference(x,h,g,kernel);r$fractions<-seq(0,1,length.out=5)
    tasks<-.np_with_seed(42L,.np_udens_q6_plan(.Random.seed,length(x),39L))
    Y<-do.call(cbind,lapply(tasks,.np_udens_q6_samples,x=x,h=h));literal<-matrix(0,39L,5L)
    for(i in seq_len(ncol(Y)))for(j in seq_along(g)){
      t<-min(Y[,i])+r$fractions[j]*diff(range(Y[,i]));density<-mean(dt((t-x)/(h*sqrt(2/3)),6))/(h*sqrt(2/3))
      literal[i,j]<-.np_udens_q6_curve(Y[,i,drop=FALSE],t,h,kernel,64L)-density
    }
    split<-do.call(rbind,lapply(tasks,.np_udens_q6_asymmetric_pointwise_chunk,reference=r))
    expect_equal(unname(split[,1:5]),literal,tolerance=1e-11)
    # Pure reducer proof uses literal task execution; actual MPI ownership has its separate installed gate.
    f<-.np_udens_q6_asymmetric_pointwise_reduce;local<-new.env(parent=environment(f))
    local$.np_udens_q6_asymmetric_pointwise_dispatch<-function(tasks,reference)do.call(rbind,lapply(tasks,.np_udens_q6_asymmetric_pointwise_chunk,reference=reference))
    environment(f)<-local
    z<-f(r,tasks,.05,TRUE);p<-f(r,tasks,.05,FALSE)
    expect_identical(z$pointwise,p$pointwise)
    expect_equal(z$pointwise$lower,r$center-apply(literal,2,manual,p=.975),tolerance=1e-11)
    expect_equal(z$pointwise$upper,r$center-apply(literal,2,manual,p=.025),tolerance=1e-11)
    expect_equal(z$raw.bonferroni$lower,r$center-apply(literal,2,manual,p=.995),tolerance=1e-11)
    expect_equal(z$raw.bonferroni$upper,r$center-apply(literal,2,manual,p=.005),tolerance=1e-11)
    old<-.np_udens_q6_asymmetric_maximum(Y,r)
    expect_equal(z$maximum[,1],old$plus,tolerance=1e-12);expect_equal(z$maximum[,2],old$minus,tolerance=1e-12)
    sim<-.np_udens_q6_asymmetric_intervals(r,old$plus,old$minus,.05)
    expect_equal(z$simultaneous$lower,sim$lower,tolerance=1e-12);expect_equal(z$simultaneous$upper,sim$upper,tolerance=1e-12)
    expect_equal(z$bonferroni$lower,pmin(z$raw.bonferroni$lower,sim$lower),tolerance=1e-12)
    expect_equal(z$bonferroni$upper,pmax(z$raw.bonferroni$upper,sim$upper),tolerance=1e-12)
    # Force two-target blocks without changing scientific draws; same complete maximum.
    local$.np_udens_q6_asymmetric_pointwise_block<-function(n,m,B)2L
    block<-f(r,tasks,.05,TRUE);expect_identical(block,z)
    expect_error(f(r,tasks,NA_real_,TRUE),"alpha")
  }
})

test_that("signed marginal reduction retains negative roots and rejects incomplete outputs", {
  f<-.np_udens_q6_asymmetric_pointwise_reduce;local<-new.env(parent=environment(f));environment(f)<-local
  r<-list(x=1:24,grid=1:2,fractions=c(0,1),center=c(1,2),se=c(.1,.2),index=1:2)
  tasks<-list(list(bsz=19L))
  roots<-cbind(seq(-3,-1,length.out=19),seq(-5,-2,length.out=19))
  local$.np_udens_q6_asymmetric_pointwise_dispatch<-function(tasks,reference)cbind(roots,0)
  z<-f(r,tasks,.1,TRUE);expect_true(all(z$pointwise$lower>r$center));expect_true(all(z$simultaneous$lower>r$center))
  local$.np_udens_q6_asymmetric_pointwise_dispatch<-function(tasks,reference)matrix(NA_real_,19,3)
  expect_error(f(r,tasks,.1,TRUE),"incomplete")
  local$.np_udens_q6_asymmetric_pointwise_dispatch<-function(tasks,reference)matrix(0,18,3)
  expect_error(f(r,tasks,.1,TRUE),"incomplete")
})
