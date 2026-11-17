test_that("M4 admission preserves q1 unbounded Epanechnikov normal and maximum owners", {
  b <- structure(list(type="fixed",ckertype="epanechnikov",ckerorder=2L,
    ckerbound="none",ndim=1L,ncon=1L),class="bandwidth")
  request <- function(z=b,control=NULL,center="bias-corrected",method="asymptotic",
                       band="simultaneous",calibration=NULL)
    .np_udens_m4_request(z,control,center,method,band,calibration)
  expect_true(request())
  for (band in c("pointwise","bonferroni","all")) expect_true(request(band=band))
  expect_false(request(band="pmzsd"))
  for (method in c("none","bootstrap")) expect_false(request(method=method))
  for (cal in c("legacy","symmetric-sd","asymmetric-sd","pmzsd"))
    expect_false(request(calibration=cal))
  expect_false(request(control=np_udens_bias_control()))
  expect_false(request(center="estimate"))
  for (change in list(list(ckertype="gaussian"),list(ckertype="uniform"),list(ckertype="beta"),
      list(type="generalized_nn"),list(type="adaptive_nn"),list(ckerbound="fixed"),
      list(ckerbound="range"),list(ckerorder=4L),list(ckerorder=6L),list(ckerorder=8L),
      list(ndim=2L,ncon=2L),list(ndim=3L,ncon=3L),list(ncon=0L)))
    expect_false(request(z=modifyList(b,change)))
  expect_false(.np_udens_m4_active(NULL))
  expect_false(.np_udens_m4_active(list(construction="density-d8-positive-pilot-v2")))
})

test_that("M4 compact correction has the stated independent polynomial moments", {
  span <- sqrt(5)
  # Independent monomial integration, with support [-sqrt(5), sqrt(5)].
  coef <- 3/(4*span)*c(2.734375, -3.28125-2.734375/5,
    .721875+3.28125/5, -.721875/5)
  moment <- function(k)sum(coef*2*span^(k+2*(0:3)+1)/(k+2*(0:3)+1))
  expect_equal(moment(0L),1,tolerance=1e-12)
  expect_equal(moment(2L),0,tolerance=1e-12)
  expect_equal(moment(4L),0,tolerance=1e-12)
  expect_gt(abs(moment(6L)),.1)
  expect_identical(.np_udens_m4_kernel(c(-2*span,-span,span,2*span)),rep(0,4))
  x <- seq(-span,span,length.out=9L)
  expect_identical(.np_udens_m4_kernel(x),.np_udens_m4_kernel(-x))
})

test_that("M4 input and workspace failures precede numerical allocation", {
  expect_error(.np_udens_m4_moments(1:9,.1,0),"geometry")
  expect_error(.np_udens_m4_moments(rep(1,10),.1,0),"normal guide")
  expect_error(.np_udens_m4_moments(1:10,0,0),"geometry")
  expect_error(.np_udens_m4_budget(100,100000),"workspace")
  expect_error(.np_udens_m4_budget(100,50,B=1e9),"workspace")
  expect_lt(.np_udens_m4_budget(1000,50,B=9999),192*1024^2)
})

# Saved-original coverage, complete covariance, installed MPI and lifecycle
# qualification are retained externally, not repeated in the default test lane.


test_that("compact normal companions use marginal and raw Bonferroni quantiles", {
  z <- list(center=c(.1,.3,.1),se=c(.2,.1,.2),bands=list(),
    settings=list(grid.size=2L,family.rows=c(1L,2L,1L)))
  p <- .np_udens_whole_normal_assemble(z,.05,"pointwise")
  b <- .np_udens_whole_normal_assemble(z,.05,"bonferroni")
  expect_equal(p$inference$bands$pointwise[,1L,1L],z$center-qnorm(.025,lower.tail=FALSE)*z$se)
  expect_equal(b$inference$bands$bonferroni[,2L,1L],z$center+qnorm(.05/4,lower.tail=FALSE)*z$se)
  expect_identical(b$inference$settings$family.rows,c(1L,2L,1L))
  z$bands$simultaneous <- array(c(z$center-2,z$center+2),c(3L,2L,1L))
  all <- .np_udens_whole_normal_assemble(z,.05,"all")
  expect_identical(all$inference$bands$pointwise,p$inference$bands$pointwise)
  expect_identical(all$inference$bands$bonferroni,b$inference$bands$bonferroni)
  expect_identical(all$inference$bands$simultaneous,z$bands$simultaneous)
  expect_lt(max(all$boot.all.err$bonferroni),min(all$boot.all.err$simultaneous))
})
