test_that("conditional inference is an explicit, validated opt-in", {
  expect_identical(np_conditional_bias_control(), np_conditional_bias_control("analytic"))
  expect_error(np_conditional_bias_control("unqualified"), "arg")
  f <- .np_conditional_control
  expect_false(f(NULL, NULL, TRUE, "bootstrap", "estimate", "default", "inid", FALSE, FALSE, FALSE))
  expect_error(f(NULL, NULL, TRUE, "bootstrap", "estimate", "default", "conditional-pilot", FALSE, FALSE, FALSE), "requires")
  bad <- np_conditional_bias_control(); bad$extra <- 1
  expect_error(f(bad, NULL, TRUE, "bootstrap", "bias-corrected", "symmetric-sd", "conditional-pilot", FALSE, FALSE, FALSE), "invalid")
})

test_that("conditional signed extremes are inverted in the correct direction", {
  # Hand-ranked roots: a large positive error must widen the lower endpoint.
  roots <- rbind(c(-2,-1,0,3,9), c(-4,-2,0,6,18))
  a <- .np_conditional_bootstrap_slice(c(10,20), roots, c(1,2), "asymmetric-sd", .4)
  expect_equal(unname(a$bands$pointwise), rbind(c(7,11),c(14,22)))
  expect_equal(unname(a$bands$simultaneous), rbind(c(7,11),c(14,22)))
  expect_equal(unname(a$bands$bonferroni), rbind(c(1,12),c(2,24)))
  # Reversing the roots reflects the interval; it must not keep the old tails.
  b <- .np_conditional_bootstrap_slice(c(10,20), -roots, c(1,2), "asymmetric-sd", .4)
  expect_equal(unname(b$bands$simultaneous), rbind(c(9,13),c(18,26)))
})

test_that("conditional workspace guards precede random generation", {
  x <- seq(-1,1,length.out=10); grid <- data.frame(x=0,y=0)
  expect_error(.np_conditional_workspace(x,x,grid,c(.3,.3),1.5), "integer")
  expect_error(.np_conditional_workspace(x,x,grid,c(.3,.3),1), "at least 2")
  expect_error(.np_conditional_workspace(x,x,data.frame(x=rep(0,257),y=1:257),c(.3,.3),20), "256 target")
})
