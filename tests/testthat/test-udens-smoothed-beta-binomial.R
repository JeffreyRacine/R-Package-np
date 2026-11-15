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
