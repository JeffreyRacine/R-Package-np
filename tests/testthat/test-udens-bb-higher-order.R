test_that("beta-binomial kernel order is exact and propagated", {
  for(order in c(4L,6L,8L)) {
    plan <- .np_udens_bb_guide_plan(c(.01,.5,.99),.12,order)
    expect_identical(plan$kernel.order,order)
    expect_identical(dim(plan$operator),c(3L,order+6L))
    value <- .np_udens_bb_guide(plan,c(2,3),100,48L)$value
    expect_equal(value,dbeta(plan$grid,2,3),tolerance=1e-10)
  }
  for(order in list(NULL,NA_real_,Inf,2L,5L,6.5,c(6L,8L),"6",TRUE,matrix(6),factor(6)))
    expect_error(.np_udens_bb_guide_plan(.5,.1,order),"order must be")
  plan <- .np_udens_bb_guide_plan(.5,.1,6L)
  plan$kernel.order <- NULL
  expect_error(.np_udens_bb_evaluate(plan,log(.4),log(.6),1),"order must be")
})

test_that("beta-binomial public spec adds only orders six and eight", {
  b <- structure(list(type="fixed",ckerbound="fixed",ckertype="epanechnikov",
    ckerorder=4L,ndim=1L,ncon=1L,ckerlb=0,ckerub=1),class="bandwidth")
  for(order in c(4L,6L,8L)) {b$ckerorder<-order;expect_null(.np_udens_bb_spec(b))}
  for(order in list(NULL,NA_real_,Inf,"6",TRUE,c(6,8),6.5,matrix(6))) {
    b$ckerorder <- order;expect_error(.np_udens_bb_spec(b),"requires q=1")
  }
  b$ckerorder<-6L
  for(key in c("type","ckerbound","ckertype")) {
    bad<-b;bad[[key]]<-switch(key,type="adaptive_nn",ckerbound="range",ckertype="uniform")
    expect_error(.np_udens_bb_spec(bad),"requires q=1")
  }
})

test_that("degree-aware native guide rejects malformed scalar degrees", {
  p <- .np_udens_bb_guide_plan(.5,.1,8L)
  for(degree in list(NA_real_,Inf,9.5,c(7,9),"9",TRUE,NULL))
    expect_error(.Call(C_np_plot_udens_bb_guide_order,p$grid,p$h,log(.4),log(.6),1,p$operator,degree),"degree must be")
  expect_error(.Call(C_np_plot_udens_bb_guide_order,p$grid,p$h,log(.4),log(.6),1,p$operator,7L),"dimensions")
})
