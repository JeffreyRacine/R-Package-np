test_that("warning conversion participates in the existing option policy", {
 ns<-asNamespace("npRmpi");snap<-get(".npRmpi_autodispatch_option_snapshot",ns)
 sync<-get(".npRmpi_autodispatch_should_sync_options",ns)
 old<-options();on.exit(options(old),add=TRUE)
 state<-get(".npRmpi_option_agreement",ns);saved<-state$snapshot
 on.exit(state$snapshot<-saved,add=TRUE);state$snapshot<-NULL
 options(npRmpi.autodispatch.option.sync="onchange")
 for(w in c(-1,0,1,2,3)) {
  options(warn=w);s<-snap();expect_identical(s$warn,as.integer(w))
  expect_true(sync(s));expect_false(sync(s))
 }
 options(npRmpi.autodispatch.option.sync="never",warn=0)
 expect_false(sync(snap()))
 options(npRmpi.autodispatch.option.sync="always")
 expect_true(sync(snap()));expect_true(sync(snap()))
})
