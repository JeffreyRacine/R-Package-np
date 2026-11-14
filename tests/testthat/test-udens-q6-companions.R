test_that('Q6 symmetric companions retain explicit law and unsupported refusals', {
 a<-list(errors='bootstrap',bootstrap='smoothed-q6',center='bias-corrected',band='all',boot.control=np_boot_control(calibration='symmetric-sd'))
 for(band in c('pointwise','simultaneous','bonferroni','all')){z<-a;z$band<-band;expect_true(.np_udens_q6_options(z,TRUE))}
 for(method in c('asymptotic','none')){z<-a;z$errors<-method;expect_error(.np_udens_q6_options(z,TRUE),'requires explicit')}
 for(center in c('estimate',NA_character_)){z<-a;z$center<-center;expect_error(.np_udens_q6_options(z,TRUE),'requires explicit')}
 z<-a;z$boot.control<-np_boot_control(calibration='asymmetric-sd');expect_true(.np_udens_q6_options(z,TRUE))
 z<-a;z$band<-'pmzsd';expect_error(.np_udens_q6_options(z,TRUE),'requires explicit')
 z<-a;z$gaussian.draws<-9999L;expect_error(.np_udens_q6_options(z,TRUE),'does not accept')
 expect_error(.np_udens_q6_options(a,FALSE),'only.*density')
 for(B in c(1L,19L,999L)){
  resolution<-.np_band_resolution(B,rep('density',4L),.05,rep(TRUE,4L),c('pointwise','simultaneous','bonferroni'))
  expect_identical(resolution$B,rep(B,3L));expect_equal(resolution$nominal_tail_draws,B*c(.05,.05,.05/4))
  expect_identical(resolution$resolution_limited,resolution$nominal_tail_draws<1)
 }
})
