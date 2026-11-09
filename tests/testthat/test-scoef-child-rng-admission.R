test_that("smooth-coefficient child admission precedes local work", {
  ns<-asNamespace("npRmpi");e<-new.env(parent=ns);hits<-0L
  refusal<-structure(list(message="admission sentinel",call=NULL),class=c("np_rng_scope_refusal","error","condition"))
  e$.npRmpi_rng_scope_collective<-function()TRUE
  e$.np_rng_scope_capture<-function(collective,deferred){expect_true(collective);expect_true(deferred);refusal}
  e$.npRmpi_with_local_regression<-function(...)stop("local work entered")
  f<-get(".npscoefbw_run_fixed_degree_bcast_payload",ns);environment(f)<-e
  expect_identical(f(stop("x forced"),stop("y forced"),NULL,NULL,NULL,NULL),refusal)
})
test_that("smooth-coefficient coordinator signals drained refusal", {
  ns<-asNamespace("npRmpi");e<-new.env(parent=ns)
  refusal<-structure(list(message="drained sentinel",call=NULL),class=c("np_rng_scope_refusal","error","condition"))
  e$.npRmpi_has_active_slave_pool<-function(...)TRUE
  e$.npRmpi_autodispatch_called_from_bcast<-function()FALSE
  e$.npRmpi_bcast_cmd_expr<-function(...)refusal
  f<-get(".npscoefbw_run_fixed_degree_collective",ns);environment(f)<-e
  old<-options(npRmpi.local.regression.mode=FALSE);on.exit(options(old))
  expect_error(f(1,1,1,1,list(),list()),"drained sentinel",class="np_rng_scope_refusal")
})
