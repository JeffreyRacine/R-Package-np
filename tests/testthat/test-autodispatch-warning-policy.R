test_that("warning conversion participates in the existing option policy", {
 ns<-asNamespace("npRmpi");snap<-get(".npRmpi_autodispatch_option_snapshot",ns)
 sync<-get(".npRmpi_autodispatch_should_sync_options",ns)
 old<-options();on.exit(options(old),add=TRUE)
 options(npRmpi.autodispatch.option.sync="onchange",npRmpi.autodispatch.option.snapshot=NULL)
 for(w in c(-1,0,1,2,3)) {
  options(warn=w);s<-snap();expect_identical(s$warn,as.integer(w))
  expect_true(sync(s));expect_false(sync(s))
 }
 options(npRmpi.autodispatch.option.sync="never",warn=0)
 expect_false(sync(snap()))
 options(npRmpi.autodispatch.option.sync="always")
 expect_true(sync(snap()));expect_true(sync(snap()))
})
