test_that("higher beta orders retain the existing q1 calibrated construction", {
  x <- data.frame(x = seq(.05, .95, length.out = 20L))
  base <- structure(list(type="fixed", ckertype="beta", ckerorder=2L,
    ckerbound="fixed", ndim=1L, ncon=1L, bw=.2, scaling=FALSE,
    ckerlb=0, ckerub=1), class="bandwidth")
  for (method in c("asymptotic", "bootstrap")) {
    for (calibration in if (method == "bootstrap")
         c("symmetric-sd", "asymmetric-sd") else NA_character_) {
      cal <- if (is.na(calibration)) NULL else calibration
      for (band in c("pointwise", "simultaneous", "bonferroni", "all")) {
        reference <- .np_udens_beta_public_prepare(base,x,method,"inid",band,cal,.05,42L)
        for (order in c(4L,6L,8L)) {
          b <- base; b$ckerorder <- order; snapshot <- serialize(b,NULL)
          got <- .np_udens_beta_public_prepare(b,x,method,"inid",band,cal,.05,42L)
          expect_identical(got,reference)
          expect_identical(serialize(b,NULL),snapshot)
        }
      }
    }
  }
})

test_that("higher beta admission retains normal q2 and unsupported-owner refusals", {
  base <- structure(list(type="fixed",ckertype="beta",ckerorder=4L,
    ckerbound="fixed",ndim=1L,ncon=1L,bw=.2,scaling=FALSE,
    ckerlb=0,ckerub=1),class="bandwidth")
  x <- data.frame(x=seq(.05,.95,length.out=20L))
  set.seed(145);before <- .Random.seed
  for (order in c(0L,1L,3L,5L,7L,9L,NA_integer_)) {
    b <- base;b$ckerorder <- order
    expect_error(.np_udens_beta_public_prepare(b,x,"asymptotic","inid",
      "all",NULL,.05,42L),"beta corrected bands require fixed bandwidths")
  }
  for (order in c(4L,6L,8L)) {
    b <- base;b$ckerorder <- order
    for (band in "none")
      expect_error(.np_udens_beta_public_prepare(b,x,"asymptotic","inid",
        band,NULL,.05,42L),"beta corrected bands require fixed bandwidths")
    for (change in list(list(ndim=2L,ncon=2L),list(ndim=3L,ncon=3L),
                       list(ncon=0L),list(type="adaptive_nn"),
                       list(type="generalized_nn"))) {
      bad <- utils::modifyList(b,change)
      expect_error(.np_udens_beta_public_prepare(bad,x,"asymptotic","inid",
        "all",NULL,.05,42L),"beta corrected bands require fixed bandwidths")
    }
    bad <- b;bad$ckerbound <- "range"
    expect_error(.np_udens_beta_public_prepare(bad,x,"asymptotic","inid",
      "all",NULL,.05,42L),"known ckerbound='fixed'",fixed=TRUE)
    for (law in c("fixed","geom","smoothed-q6"))
      expect_error(.np_udens_beta_public_prepare(b,x,"bootstrap",law,
        "all","symmetric-sd",.05,42L),"require bootstrap='inid'",fixed=TRUE)
  }
  expect_identical(.Random.seed,before)
  base$ckerorder <- 2L
  expect_true(.np_udens_beta_order_request(base,"pmzsd"))
  base$ndim <- 2L;expect_true(.np_udens_beta_order_request(base,"all"))
})


test_that("higher beta pmzsd uses the unchanged marginal normal state", {
  x <- data.frame(x=seq(.05,.95,length.out=20L))
  b <- structure(list(type="fixed",ckertype="beta",ckerorder=2L,
    ckerbound="fixed",ndim=1L,ncon=1L,bw=.2,scaling=FALSE,
    ckerlb=0,ckerub=1),class="bandwidth")
  prepare <- function(b,method="asymptotic",calibration=NULL,alpha=.05)
    .np_udens_beta_public_prepare(b,x,method,"inid","pmzsd",calibration,alpha,42L)
  ref <- prepare(b)
  expect_true(ref$normal.pmzsd)
  for (order in c(4L,6L,8L)) {
    b$ckerorder <- order
    expect_identical(prepare(b),ref)
    expect_error(prepare(b,"bootstrap","symmetric-sd"),"asymptotic")
    expect_error(prepare(b,calibration="symmetric-sd"),"calibration")
    for (alpha in c(0,.5,1)) expect_error(prepare(b,alpha=alpha),"alpha")
    b2 <- b;b2$ndim <- b2$ncon <- 2L
    expect_error(prepare(b2),"fixed bandwidths")
  }
})
