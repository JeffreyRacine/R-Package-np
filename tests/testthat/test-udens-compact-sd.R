test_that('compact SD admission requires the explicit commissioned route', {
  b<-structure(list(ckertype='uniform',ndim=1L,ncon=1L,type='fixed',ckerorder=2L,ckerbound='fixed'),class='bandwidth')
  for(cal in c('symmetric-sd','asymmetric-sd')){
    expect_true(.np_udens_compact_request(b,NULL,'bias-corrected','bootstrap',cal,TRUE))
    expect_false(.np_udens_compact_request(b,NULL,'bias-corrected','bootstrap',cal,FALSE))
    expect_false(.np_udens_compact_request(b,NULL,'estimate','bootstrap',cal,TRUE))
    expect_false(.np_udens_compact_request(b,list(version=3L),'bias-corrected','bootstrap',cal,TRUE))
    for(key in c('ndim','ncon')){z<-b;z[[key]]<-2L;expect_false(.np_udens_compact_request(z,NULL,'bias-corrected','bootstrap',cal,TRUE))}
    for(change in list(list(type='generalized_nn'),list(ckerorder=4L),list(ckerbound='none'),list(ckertype='gaussian'))){z<-b;z[names(change)]<-change;expect_false(.np_udens_compact_request(z,NULL,'bias-corrected','bootstrap',cal,TRUE))}
  }
  for(cal in list(NULL,'legacy','pmzsd','smoothed-percentile'))expect_false(.np_udens_compact_request(b,NULL,'bias-corrected','bootstrap',cal,TRUE))
  expect_true(.np_udens_compact_request(b,NULL,'bias-corrected','asymptotic'))
})

test_that('compact SD uses complete signed interval envelopes and shared all components', {
  x<-qbeta((1:20-.35)/21,2,3);grid<-c(.14,.43,.81);B<-7L
  for(kernel in c('uniform','epanechnikov'))for(cal in c('symmetric-sd','asymmetric-sd')){
    call<-function(band,g=grid)suppressWarnings(.np_udens_compact_hidden(x,.18,g,kernel,.1,band,B,19L,variant=cal))
    all<-call('all')
    expect_identical(all$family_size,3L)
    expect_identical(names(all$bands),c('pointwise','simultaneous','bonferroni'))
    for(b in names(all$bands)){
      z<-call(b);expect_identical(z$center,all$center)
      expect_identical(z$lower,all$bands[[b]]$lower);expect_identical(z$upper,all$bands[[b]]$upper)
    }
    dup<-call('all',grid[c(3,1,1,2)])
    expect_identical(dup$center,all$center[c(3,1,1,2)])
    expect_identical(dup$lower,all$lower[c(3,1,1,2)])
    expect_identical(dup$family_size,3L)
    expect_true(all(all$bands$bonferroni$lower<=all$bands$simultaneous$lower))
    expect_true(all(all$bands$bonferroni$upper>=all$bands$simultaneous$upper))
    expect_true(all(is.finite(c(all$lower,all$upper))))
    expect_true(all$sd_companion$resolution.limited)
  }
  store<-.np_udens_beta_store();on.exit(.np_udens_beta_store_close(store))
  z<-rbind(as.double(1:7),as.double(2:8));root<-.np_udens_beta_root_store(store,2L,7L);root<-.np_udens_beta_root_write(root,z,1L)
  scales<-list(c(1,2),c(4,3));region<-list(lower=c(-1,-1),upper=c(1,1))
  got<-.np_udens_bounded_sd(root,c(0,0),scales,region,c(0,0),.1,'all','asymmetric-sd')$bands$pointwise
  lo<-pmin(-1-scales[[1]]*c(7,8),-1-scales[[2]]*c(7,8))
  hi<-pmax(1-scales[[1]]*c(1,2),1-scales[[2]]*c(1,2))
  expect_identical(got$lower,lo);expect_identical(got$upper,hi)
  expect_false(identical(got$upper,1-pmax(scales[[1]],scales[[2]])*c(1,2)))
})

test_that('compact SD shared projection certifies tiles and fails closed on resources', {
  plan<-.np_udens_compact_operator_plan(c(.01,.35,.8),.13,'epanechnikov')
  theta<-rbind(c(0,0),c(1,2),c(-.3,.4))
  a<-.np_udens_compact_common(plan,theta,128L,fit.block=1L,target.block=1L)
  b<-.np_udens_compact_common(plan,theta,256L,fit.block=2L,target.block=2L)
  expect_true(max(abs(a$mu-b$mu))<=1e-12+1e-10*max(abs(b$mu)))
  expect_identical(a$diagnostics$base,apply(theta,2,min))
  expect_lte(a$diagnostics$bytes,256*1024^2)
  expect_error(.np_udens_compact_common(plan,theta,128L,held.bytes=1024^3),'workspace')
  expect_error(.np_udens_compact_common(plan,matrix(-1,2,2),128L),'projection')
  store<-.np_udens_beta_store();on.exit(.np_udens_beta_store_close(store))
  root<-.np_udens_beta_root_store(store,2L,7L)
  expect_error(.np_udens_beta_root_read(root,1L),'incomplete')
  x<-qbeta((1:20-.35)/21,2,3)
  expect_error(.np_udens_compact_hidden(x,.13,c(.2,.8),'uniform',.1,'all',1e9,19L,variant='asymmetric-sd'),'budget')
  expect_error(.np_udens_compact_hidden(x,.13,0,'uniform',.1,'all',7L,19L,variant='asymmetric-sd'),'request')
})
