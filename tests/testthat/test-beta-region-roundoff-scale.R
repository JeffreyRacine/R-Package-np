test_that("sphere roundoff uses the objective arithmetic scale", {
 f<-.np_udens_beta_sphere_accept
 # Saved cancellation witness: the proposed step meets the original gtol.
 v<--.07019149167247;p<--2.358932035509313e-14
 expect_false(f(v,v+6.66e-16,p,7.3e-13,1.192083e-7,1e-7))
 expect_true(f(v,v+6.66e-16,p,7.3e-13,1.192083e-7,1e-7,45.55))
 expect_false(f(v,v+1e-7,p,7.3e-13,1.192083e-7,1e-7,45.55))
 expect_false(f(v,v+6.66e-16,p,1.1e-7,1.192083e-7,1e-7,45.55))
 expect_false(f(v,v+6.66e-16,-1e-5,7.3e-13,1.192083e-7,1e-7,45.55))
 # Ordinary sufficient decrease is independent of the roundoff scale.
 expect_true(f(1,.9,-.1,1,2,1e-7))
 expect_identical(f(1,.9,-.1,1,2,1e-7),f(1,.9,-.1,1,2,1e-7,100))
})
