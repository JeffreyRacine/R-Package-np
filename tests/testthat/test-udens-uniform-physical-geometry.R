test_that('uniform compact observed rows retain physical open support', {
  x<-21:65;grid<-26:60;h<-5;lo<-20;width<-46
  plan<-.np_udens_compact_operator_plan((grid-lo)/width,h/width,'uniform')
  physical<-list(x=matrix(x,ncol=1L),grid=matrix(grid,ncol=1L),h=h)
  z<-outer(grid,x,function(t,v)(v-t)/h)
  expected<-3/8*(3-5*z^2)/h
  expected[abs(outer(grid,x,'-'))>=h]<-0
  old<-.np_udens_compact_operator(plan,(x-lo)/width)/width
  got<-.np_udens_compact_operator(plan,(x-lo)/width,physical)/width
  expect_gt(max(abs(old-expected)),.01)
  expect_equal(got,expected,tolerance=1e-12)
  # Slicing must carry training and target identity independently.
  geo<-.np_udens_compact_geometry_slice(physical,rows=3:5,observations=4:9)
  pl<-.np_udens_compact_operator_plan((grid[3:5]-lo)/width,h/width,'uniform')
  expect_equal(.np_udens_compact_operator(pl,(x[4:9]-lo)/width,geo)/width,
    expected[3:5,4:9],tolerance=1e-12)
  expect_error(.np_udens_compact_operator(pl,(x-lo)/width,geo),'shape')
})

test_that('physical target identity and q2 observed-axis rows survive normalization', {
  raw<-matrix(c(.5,.5+.Machine$double.eps,.5),ncol=1L)
  z<-.np_udens_compact_physical_grid(rep(.5,3),list(x=matrix(1:10,ncol=1L),grid=raw,h=.1))
  expect_length(z$grid,2L);expect_identical(z$restore,c(1L,2L,1L))
  x<-as.matrix(expand.grid(a=(1:7)/8,b=(1:7)/8));h<-c(.125,.125)
  physical<-list(x=x,grid=x,h=h)
  da<-lapply(1:2,function(j)list(lv=log(x[,j]),l1v=log1p(-x[,j])))
  cache<-.np_udens_compact_normal_q2_axis('uniform',physical)(x,da,h)
  observed<-.np_udens_beta_joint_block_native(cache,seq_len(nrow(x)),seq_len(nrow(x)))
  parts<-lapply(1:2,function(j){
    pl<-.np_udens_compact_operator_plan(x[,j],h[j],'uniform')
    .np_udens_compact_operator(pl,exp(log(x[,j])),
      .np_udens_compact_geometry_slice(physical,columns=j))
  })
  expect_equal(observed,parts[[1]]*parts[[2]]-1,tolerance=1e-12)
  factory<-.np_udens_compact_observed_axis(x,h,x,'uniform',physical)
  subset<-.np_udens_beta_joint_block_native(factory(c(4L,1L)),seq_len(nrow(x)),seq_len(nrow(x)))
  expect_identical(subset,observed[c(4L,1L),,drop=FALSE])
})
