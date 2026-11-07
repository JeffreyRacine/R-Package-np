# Small deterministic reducer tests; public/coverage/timing matrices live outside
# the default lane. Existing beta tests protect the original SD reducers.
test_that('bounded SD companions preserve components and invert signed roots', {
  roots<-rbind(seq(-4,2,length.out=99),seq(-1,5,length.out=99),rep(c(-2,0,3),33))
  store<-.np_udens_beta_store();on.exit(.np_udens_beta_store_close(store),add=TRUE)
  root<-.np_udens_beta_root_store(store,3L,99L)
  root<-.np_udens_beta_root_write(root,roots,1L)
  for(law in c('symmetric-sd','asymmetric-sd')) {
    region<-list(lower=c(-.2,.1,.3),upper=c(.4,.5,.8))
    center<-c(.2,.7,.4);rhat<-c(-.1,.3,.2);scale<-c(.5,2,.8);alpha<-.1
    got<-.np_udens_bounded_sd(root,rhat,list(scale),region,center,alpha,'all',law)
    pw<-.np_udens_beta_stored_reduce(root,rhat,scale,region,center,alpha,'pointwise',law)
    si<-.np_udens_beta_stored_reduce(root,rhat,scale,region,center,alpha,'simultaneous',law)
    expect_identical(got$bands$pointwise,list(lower=pw$lower,upper=pw$upper))
    expect_identical(got$bands$simultaneous,list(lower=si$lower,upper=si$upper))
    if(law=='asymmetric-sd') {
      k<-floor(ncol(roots)*(alpha/3)/4)+1L
      low<-apply(roots,1,sort)[k,];high<-apply(roots,1,sort)[100L-k,]
      lower<-rhat+region$lower-scale*high;upper<-rhat+region$upper-scale*low
    } else {
      k<-ceiling(ncol(roots)*(1-(alpha/3)/2))
      z<-apply(abs(roots),1,sort)[k,]
      raw.lower<-rhat+region$lower-scale*z;raw.upper<-rhat+region$upper+scale*z
      h<-.np_udens_beta_hull(center,raw.lower,raw.upper);lower<-h[,1];upper<-h[,2]
    }
    expect_identical(got$bands$bonferroni$lower,pmin(lower,si$lower))
    expect_identical(got$bands$bonferroni$upper,pmax(upper,si$upper))
    expect_true(all(got$bands$bonferroni$lower<=si$lower & got$bands$bonferroni$upper>=si$upper))
    direct<-.np_udens_bounded_sd(root,rhat,list(scale),region,center,alpha,'bonferroni',law)
    expect_identical(direct$lower,got$bands$bonferroni$lower)
    expect_identical(direct$upper,got$bands$bonferroni$upper)
    expect_true(direct$sd_companion$resolution.limited)
    expect_equal(direct$sd_companion$expected.tail.draws,99*alpha/(if(law=='asymmetric-sd')12 else 6))
  }
})

test_that('bounded SD singleton and constituent envelope preserve the complete intervals', {
  store<-.np_udens_beta_store();on.exit(.np_udens_beta_store_close(store),add=TRUE)
  root<-.np_udens_beta_root_store(store,1L,19L)
  root<-.np_udens_beta_root_write(root,matrix(seq(1,4,length.out=19),1),1L)
  for(law in c('symmetric-sd','asymmetric-sd')) {
    region<-list(lower=2,upper=3)
    got<-.np_udens_bounded_sd(root,1,list(.5,2),region,3,.05,'all',law)
    expect_identical(got$bands$pointwise,got$bands$simultaneous)
    expect_identical(got$bands$bonferroni,got$bands$simultaneous)
    pieces<-lapply(c(.5,2),function(s).np_udens_beta_stored_reduce(root,1,s,region,3,.05,'simultaneous',law))
    expect_identical(got$lower,pmin(pieces[[1]]$lower,pieces[[2]]$lower))
    expect_identical(got$upper,pmax(pieces[[1]]$upper,pieces[[2]]$upper))
  }
  expect_error(.np_udens_bounded_sd(root,1,list(1),list(lower=0,upper=1),0,.05,'pmzsd','symmetric-sd'),'invalid bounded')
})
