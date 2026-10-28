# Tiny default-lane contracts. Public/MPI matrices stay in the retained harness.
test_that("bootstrap calibration validates names and preserves old controls", {
  expect_identical(np_boot_control()$calibration,"default")
  old <- np_boot_control();old$calibration <- NULL
  expect_identical(.np_plot_boot_calibration(old),"default")
  for (cal in c("absolute","separate-tails","tail-envelope"))
    expect_identical(np_boot_control(calibration=cal)$calibration,cal)
  expect_error(np_boot_control(calibration="rank"),"calibration")
  expect_error(np_boot_control(calibration=NULL),"calibration")
  old$calibration <- c("default","absolute")
  expect_error(.np_plot_boot_calibration(old),"calibration")
})

test_that("signed tails invert orientation and preserve singleton nesting", {
  x <- matrix(c(-4,-1,0,2,7),1L)
  map <- data.frame(row=1L,target="mean")
  z <- .np_band_tail_prepare(5L,map,.05)
  z <- .np_band_tail_update(z,x,map,seq_len(5L))
  z <- .np_band_tail_finish(z)
  expect_identical(as.numeric(z$limits$pointwise),c(-4,7))
  expect_identical(z$limits$pointwise,z$limits$simultaneous)
  expect_identical(z$limits$pointwise,z$limits$bonferroni)
  y <- .np_band_rank_invert(10,2,z$limits)
  expect_identical(as.numeric(y$pointwise),c(-4,18))
  expect_error(.np_band_tail_prepare(0,map,.05),"draw count")
})

test_that("low resolution diagnoses all bands and warns once", {
  z <- .np_band_resolution(1L,c("mean","mean"),.05,c(TRUE,TRUE),
    c("pointwise","simultaneous","bonferroni"),TRUE,matrix(1L,1L))
  expect_true(all(z$computed & z$resolution_limited))
  expect_equal(nrow(z),3L)
  expect_warning(.np_band_resolution_warn(z),class="np_band_resolution_warning")
  expect_silent(.np_band_resolution_warn(transform(z,resolution_limited=FALSE)))
})

test_that("pilot block scale has the independent-block identity", {
  F <- rbind(c(-2,-1,1,2),rep(0,4))
  for (law in c("fixed","geom")) {
    z <- .np_bias_pilot_block_scale(F,c(3,0),1L,law)
    expect_equal(z$se,c(3,0),tolerance=1e-14)
    expect_equal(z$variance.ratio,c(1,0),tolerance=1e-14)
  }
  expect_error(.np_bias_pilot_block_scale(F,c(3,0),4L,"fixed"),"variance")
})
