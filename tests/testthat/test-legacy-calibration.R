test_that("legacy calibration resolves only the retained family controls", {
  control <- np_boot_control(calibration="legacy")
  expect_identical(.np_plot_boot_calibration(control),"legacy")
  constructors <- list(regression=np_bias_control,density=np_udens_bias_control,
                       distribution=np_udist_bias_control)
  for(family in names(constructors)) {
    old <- constructors[[family]](construction="legacy")
    select <- function(bias=NULL,center="bias-corrected",method="bootstrap")
      .np_plot_legacy_bias_control("legacy",bias,family,center,method)
    expect_identical(select(),old)
    expect_identical(select(old),old)
    expect_error(select(constructors[[family]]()),"conflicts")
    expect_error(select(center="estimate"),"requires")
    expect_error(select(method="asymptotic"),"requires")
    expect_null(.np_plot_legacy_bias_control("default",NULL,family,"estimate","none"))
  }
})
