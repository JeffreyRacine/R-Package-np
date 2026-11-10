test_that("compact q2 normal admission is exact and leaves neighboring owners intact", {
  b<-structure(list(ckertype="uniform",type="fixed",ckerbound="fixed",ckerorder=2L,ndim=2L,ncon=2L),class="bandwidth")
  request<-function(b,control=NULL,center="bias-corrected",method="asymptotic",band="pmzsd",calibration=NULL)
    .np_udens_compact_normal_q2_request(b,control,center,method,band,calibration)
  expect_true(request(b));expect_true(request(modifyList(b,list(ckertype="epanechnikov"))))
  for(change in list(list(ckertype="gaussian"),list(ckertype="beta"),list(ckerorder=4L),
    list(type="adaptive_nn"),list(type="generalized_nn"),list(ckerbound="range"),
    list(ckerbound="none"),list(ncon=1L),list(ndim=1L),list(ndim=3L)))
    expect_false(request(modifyList(b,change)))
  for(band in c("pointwise","simultaneous","bonferroni","all"))expect_false(request(b,band=band))
  expect_false(request(b,method="bootstrap"));expect_false(request(b,center="estimate"))
  expect_false(request(b,control=np_udens_bias_control()));expect_false(request(b,calibration="symmetric-sd"))
  expect_false(.np_udens_compact_normal_q2_active(NULL))
  expect_false(.np_udens_compact_normal_q2_active(list(construction="density-fixed-compact-cubic-v1")))
  # q1 admission remains the old explicit family; the new selector does not subsume it.
  b1<-modifyList(b,list(ndim=1L,ncon=1L));expect_true(.np_udens_compact_request(b1,NULL,"bias-corrected","asymptotic"))
})

test_that("compact q2 capacity refuses before oversized target tensors", {
  cap<-.np_udens_compact_normal_q2_capacity
  x<-cap(20,9,c(1600,1600));remaining<-192*1024^2-x
  expect_lte(cap(20,9,c(1600,1600),remaining),192*1024^2)
  expect_error(cap(20,9,c(1600,1600),remaining+1),"workspace")
  expect_error(cap(20,9,c(1601,1600)),"node cap")
  expect_error(cap(20,10),"geometry");expect_error(cap(1e9,9),"workspace")
})

# Saved-bank, default surfaces, literal influence, MPI and timings are external
# qualification. These small tests add unique admission and allocation contracts.
