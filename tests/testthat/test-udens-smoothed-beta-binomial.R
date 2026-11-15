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
