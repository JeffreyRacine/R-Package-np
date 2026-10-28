test_that("GNN bootstrap collapsed margins use the exact CDF limit", {
  delta <- rbind(c(-2, 0, 3), c(-1, .2, 1))
  z <- .np_udist_gnn_zero_margin(delta, c(0, .5))
  expect_identical(z$H[1, ], c(0, .5, 1))
  expect_identical(z$D[1, ], rep(0, 3))
  expect_identical(z$H[2, ], pnorm(delta[2, ]/.5))
  expect_identical(z$D[2, ], (delta[2, ]/.5)*dnorm(delta[2, ]/.5))
  prepared <- list(n=3L, category=rep(1,3), corrected=TRUE,
    map=list(1L), cache=list(list(delta=matrix(c(0,-1,-2),1),
    order=matrix(1:3,1),k=2L)))
  expect_error(.np_udist_gnn_center(prepared,c(3,0,0)),"nonpositive")
  prepared$zero.limit <- TRUE
  expect_identical(.np_udist_gnn_center(prepared,c(3,0,0))$center,.5)
  expect_identical(.np_udist_gnn_center(prepared,c(3,0,0))$radii,matrix(0,1,1))
})
