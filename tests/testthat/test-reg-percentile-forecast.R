test_that("draw-aware capacity is checked before corrected-row capture", {
  expect_error(.np_reg_percentile_prepare(list(d=2L,y=rep(0,150)),NULL,
    data.frame(x=seq_len(2500)),B=999L),"memory forecast")
  expect_equal(.np_reg_percentile_forecast(100,100,999),
    8*(10*100*999+10*100*39+6*100*100+20*100*256))
  expect_error(.np_reg_percentile_forecast(100,Inf,999),"dimensions")
})


test_that("requested percentile families validate before preparation", {
  expect_error(.np_reg_percentile_prepare(list(d=2L),NULL,data.frame(x=1),
    targets="d3"),"target family")
  expect_error(.np_reg_percentile_prepare(list(d=2L),NULL,data.frame(x=1),
    targets=c("mean","mean")),"target family")
})
