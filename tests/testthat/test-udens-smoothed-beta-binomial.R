test_that("beta-binomial guide preserves polynomial reproduction at boundaries", {
  grid <- c(.001,.5,.999)
  plan <- .np_udens_bb_guide_plan(grid,.08)
  for(ab in list(c(1,1),c(2,3),c(4,2))) {
    low <- .np_udens_bb_guide(plan,ab,80,32L)$value
    high <- .np_udens_bb_guide(plan,ab,100,48L)$value
    expect_equal(low,dbeta(grid,ab[1],ab[2]),tolerance=1e-10)
    expect_equal(high,low,tolerance=1e-10)
  }
  expect_error(.np_udens_bb_guide_plan(c(0,.5),.1),"interior")
  expect_error(.np_udens_bb_guide_plan(.5,0),"positive bandwidth")
  expect_error(.np_udens_bb_guide(plan,c(0,1),100,48L),"guide rule")
  expect_error(.np_udens_bb_guide(plan,c(1,1),100,64L),"guide rule")
})

test_that("beta-binomial complete factor refuses malformed mathematical inputs", {
  expect_error(.np_udens_bb_factor(c(.5,.5),matrix(-1,2,2),matrix(0,1,2),
                                  matrix(0,1,2),c(-1,-1),10),"dimensions")
  expect_error(.np_udens_bb_check(1,1.01,"oracle"),"refinement failed")
  expect_error(.np_udens_bb_check(c(1,2),1,"oracle"),"invalid")
})

test_that("beta-binomial dual-log draws preserve the R stream and reject before RNG", {
  result <- .np_with_seed(271828L, {
    entry <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
    a <- rep(c(.4,1,2,10),10);b <- rev(a)
    original <- rbeta(length(a),a,b);terminal <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
    assign(".Random.seed",entry,.GlobalEnv)
    paired <- .np_udens_bb_draw(a,b)
    expect_identical(paired[,1L],original)
    expect_identical(get(".Random.seed",.GlobalEnv,inherits=FALSE),terminal)
    expect_true(all(is.finite(paired[,2:3])))
    expect_lt(max(abs(exp(paired[,2])+exp(paired[,3])-1)),1e-14)
    expect_error(.np_udens_bb_draw(c(1,-1),c(1,1)),"positive")
    expect_identical(get(".Random.seed",.GlobalEnv,inherits=FALSE),terminal)
    TRUE
  })
  expect_true(result)
})

test_that("paired-log guide adapter preserves the ordinary complete fit", {
  x <- (1:20)/21
  expect_identical(.np_udens_bb_fit(x,log(x),log1p(-x)),.np_udens_beta_q1_fit(x))
  expect_error(.np_udens_bb_fit(x,log(x),log1p(-x)+.01),"paired log")
})

test_that("beta-binomial admission adds only its explicit statistical tuple", {
  dots <- list(errors="bootstrap",bootstrap="smoothed-beta-binomial",
    center="bias-corrected",band="all",B=39L,
    boot.control=np_boot_control(calibration="symmetric-sd"))
  expect_true(.np_udens_bb_options(dots,TRUE))
  expect_error(.np_udens_bb_options(dots,FALSE),"unconditional density")
  for(key in c("errors","center","band","boot.control")) {
    bad <- dots;bad[[key]] <- NULL
    expect_error(.np_udens_bb_options(bad,TRUE),"requires explicit")
  }
  for(cal in c("asymmetric-sd","default","absolute","smoothed-percentile")) {
    bad <- dots;bad$boot.control <- np_boot_control(calibration=cal)
    expect_error(.np_udens_bb_options(bad,TRUE),"symmetric-sd")
  }
  for(band in c("pointwise","simultaneous","bonferroni","all")) {
    good <- dots;good$band <- band;expect_true(.np_udens_bb_options(good,TRUE))
  }
  for(value in list(0,1.5,Inf,NA_real_)) {
    bad <- dots;bad$B <- value;expect_error(.np_udens_bb_options(bad,TRUE),"integer")
  }
  b <- structure(list(type="fixed",ckerbound="fixed",ckertype="epanechnikov",
    ckerorder=4L,ndim=1L,ncon=1L,ckerlb=0,ckerub=1),class="bandwidth")
  expect_null(.np_udens_bb_spec(b))
  changes <- list(type=c("adaptive_nn","generalized_nn"),ckerbound=c("none","range"),
    ckertype=c("uniform","gaussian","beta"),ckerorder=c(2L,6L,8L),ndim=c(2L,3L),ncon=0L)
  for(key in names(changes))for(value in changes[[key]]) {
    bad <- b;bad[[key]] <- value;expect_error(.np_udens_bb_spec(bad),"requires q=1")
  }
  for(value in list(Inf,NA_real_,-1)) {
    bad <- b;bad$ckerub <- value;expect_error(.np_udens_bb_spec(bad),"finite ordered")
  }
})

test_that("beta-binomial tail assembly uses complete absolute roots and explicit hull", {
  r <- list(grid=c(.1,.9));draw <- cbind(c(-9,-2,1,3,4),c(1,-3,2,7,5),0)
  z <- .np_udens_bb_reduce(r,draw,.1,"all")
  absolute <- abs(draw[,1:2])
  expect_identical(z$pointwise,unname(apply(absolute,2,quantile,probs=.9,type=7)))
  expect_identical(z$simultaneous,unname(quantile(apply(absolute,1,max),.9,type=7)))
  expect_identical(z$bonferroni.raw,unname(apply(absolute,2,quantile,probs=.95,type=7)))
  expect_identical(z$bonferroni,pmax(z$bonferroni.raw,z$simultaneous))
  expect_error(.np_udens_bb_reduce(r,draw,.05,"pmzsd"),"incomplete")
  expect_error(.np_udens_bb_budget(100,100000,999),"workspace")
})
