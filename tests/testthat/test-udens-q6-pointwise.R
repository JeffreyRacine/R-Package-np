test_that('Q6 marginal roots equal direct conditional curves without a maximum', {
 .np_with_seed(519L, {
  x<-qnorm((1:20-.5)/20);h<-.6;B<-31L;start<-.Random.seed
  tasks<-.np_udens_q6_plan(start,20L,B,7L)
  Y<-do.call(cbind,lapply(tasks,.np_udens_q6_samples,x=x,h=h))
  for(kernel in c('uniform','epanechnikov')) {
   ref<-.np_udens_q6_reference(x,h,c(-.8,0,.6),kernel)
   roots<-sapply(ref$grid,function(t).np_udens_q6_curve(Y,t,h,kernel,64L))
   roots<-abs(sweep(sweep(roots,2,ref$reference.center),2,ref$se,'/'))
   chunks<-do.call(rbind,lapply(tasks,.np_udens_q6_pointwise_chunk,reference=ref))
   expect_equal(chunks[,1:3],roots,tolerance=1e-12)
   expect_equal(apply(roots,1,max),.np_udens_q6_maximum(Y,ref)$maximum,tolerance=1e-12)
   marginal<-apply(roots,2,quantile,.95,type=7)
   expect_true(all(marginal<=quantile(apply(roots,1,max),.95,type=7)))
   # A target permutation is purely a reindexing, including duplicate rows.
   rp<-ref;rp$grid<-ref$grid[c(3,1,3)];rp$reference.center<-ref$reference.center[c(3,1,3)];rp$se<-ref$se[c(3,1,3)]
   cp<-do.call(rbind,lapply(tasks,.np_udens_q6_pointwise_chunk,reference=rp))
   expect_identical(cp[,1:3],chunks[,c(3,1,3)])
  }
 })
})
test_that('pointwise admission leaves other unsupported bands rejected', {
 a<-list(errors='bootstrap',bootstrap='smoothed-q6',center='bias-corrected',band='pointwise',boot.control=np_boot_control(calibration='symmetric-sd'))
 expect_true(.np_udens_q6_options(a,TRUE))
 expect_error(.np_udens_q6_options(a,FALSE),'only.*density')
 for(b in c('pmzsd')){z<-a;z$band<-b;expect_error(.np_udens_q6_options(z,TRUE),'requires explicit')}
 z<-a;z$boot.control<-np_boot_control(calibration='asymmetric-sd');expect_error(.np_udens_q6_options(z,TRUE),'requires explicit')
 expect_error(.np_udens_q6_pointwise_budget(100,50,1e8),'workspace bound')
 expect_true(.np_udens_q6_pointwise_budget(100,50,9999)<192*1024^2)
})

test_that('Q6 pointwise blocks respect the same fixed memory envelope', {
 for(n in c(20L,200L,1000L)) for(B in c(31L,999L,9999L,120000L)) {
  m<-50L
  k<-.np_udens_q6_pointwise_block(n,m,B)
  expect_true(k>=min(8L,m) && k<=m)
  expect_true(.np_udens_q6_pointwise_budget(n,m,B,k)<=192*1024^2)
  if(k<m) expect_error(.np_udens_q6_pointwise_budget(n,m,B,k+1L),'workspace bound')
 }
 expect_equal(.np_udens_q6_pointwise_block(200L,50L,9999L),50L)
 expect_error(.np_udens_q6_pointwise_block(100L,50L,1e8),'workspace bound')
 expect_error(.np_udens_q6_pointwise_block(20L,5L,0L),'replicate count')
})
