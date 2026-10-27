test_that("local calibration neighborhoods preserve units, rows and constant variance", {
  owner <- np:::.np_bias_local_neighborhood
  x <- matrix(seq(-1,1,length.out=24),ncol=1)
  q <- seq(.25,1,length.out=24)
  a <- owner(x,q)
  b <- owner(7+3*x,q)
  ii <- c(1L,12L,24L)
  expect_equal(a$raw(ii),b$raw(ii),tolerance=1e-13)
  expect_equal(drop(a$raw(ii) %*% (4*q)),rep(4,3),tolerance=1e-13)
  permutation <- rev(seq_len(24))
  c <- owner(x[permutation,,drop=FALSE],q[permutation])
  expect_equal(c$raw(match(ii,permutation))[,match(seq_len(24),permutation)],
               a$raw(ii),tolerance=1e-13)
  expect_error(a$raw(seq_len(17)),"invalid local calibration neighborhood tile")
  expect_error(owner(matrix(1,24,1),q),"invalid local calibration neighborhood geometry")
})
