test_that("beta q2 prepared axes bind their complete family and measure", {
  theta <- c(0,0,0,0,0); h <- c(.2,.3)
  grid <- rbind(c(.25,.3),c(.75,.8))
  p <- .np_udens_beta_q2_axes(theta,h,grid,32L)
  expect_identical(.np_udens_beta_q2_axes_check(p,theta,h,grid[2:1,,drop=FALSE],32L),p$axes)
  expect_equal(vapply(p$axes,function(z)sum(z$prob),0.),c(1,1),tolerance=1e-14)
  expect_error(.np_udens_beta_q2_axes_check(p,theta,h*2,grid,32L),"incompatible")
  expect_error(.np_udens_beta_q2_axes_check(p,theta,h,grid,36L),"incompatible")
  expect_error(.np_udens_beta_q2_axes_check(p,theta,h,matrix(c(.1,.2),1L),32L),"invalid")
  bad <- p; bad$axes[[1L]]$prob[1L] <- -1
  expect_error(.np_udens_beta_q2_axes_check(bad,theta,h,grid,32L),"invalid")
  bad <- p; bad$axes[[1L]]$mean[1L] <- bad$axes[[1L]]$mean[1L]+.01
  expect_error(.np_udens_beta_q2_axes_check(bad,theta,h,grid,32L),"prepared axis mean")
  expect_error(.np_udens_beta_q2_axes(theta,h,matrix(c(0,.5),1L),32L),"invalid")
  expect_error(.np_udens_beta_q2_axes(theta,c(1e-8,.3),grid,32L),"node budget")
  # Existing native owners admit -Inf log weights for negligible tail mass.
  concentrated <- .np_udens_beta_q2_axes(c(99,99,0,0,0),h,grid,32L)
  expect_true(any(concentrated$axes[[1L]]$prob==0))
  expect_true(all(is.finite(concentrated$axes[[1L]]$T)))
  # Higher calibrated q2 admission is separate from the retained normal refusal.
  bws <- list(ckerorder=4L,ndim=2L,ncon=2L,ckertype="beta",type="fixed",ckerbound="fixed")
  expect_true(.np_udens_beta_order_request(bws,"all"))
  expect_false(.np_udens_beta_normal_q2_request(bws))
})

test_that("beta q2 preparation storage is exact and invocation owned", {
  store <- .np_udens_beta_store()
  on.exit(.np_udens_beta_store_close(store),add=TRUE)
  prep <- list(theta=c(0,0,0,0,0),h=c(.2,.3),grid=matrix(c(.2,.3),1L),
               inner=64L,outer=32L,ratio=matrix(pi,3L,4L))
  record <- .np_udens_beta_q2_preparation_write(prep,store,"preparation-test.rds")
  expect_identical(.np_udens_beta_q2_preparation_read(record),prep)
  bad <- record; bad$h <- bad$h*2
  expect_error(.np_udens_beta_q2_preparation_read(bad),"identity")
  bad <- record; bad$size <- bad$size+1
  expect_error(.np_udens_beta_q2_preparation_read(bad),"size")
  path <- store$path
  .np_udens_beta_store_close(store)
  expect_false(dir.exists(path))
})


test_that("beta q2 axes preserve coincident rounded nodes without changing weights", {
  theta <- rep(0,5); h <- c(.2,.15); grid <- rbind(c(.25,.3),c(.75,.8))
  p <- .np_udens_beta_q2_axes(theta,h,grid,32L)
  expect_true(any(diff(p$axes[[2L]]$lv-p$axes[[2L]]$l1v)==0))
  expect_identical(.np_udens_beta_q2_axes_check(p,theta,h,grid,32L),p$axes)
  bad <- p; z <- bad$axes[[2L]]; ix <- c(2L,1L,seq.int(3L,z$nodes))
  for(k in c("lv","l1v","prob","log.weights"))z[[k]] <- z[[k]][ix]
  z$T <- z$T[ix,,drop=FALSE]; bad$axes[[2L]] <- z
  expect_error(.np_udens_beta_q2_axes_check(bad,theta,h,grid,32L),"invalid")
})
