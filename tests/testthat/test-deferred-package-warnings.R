test_that('only named warnings defer; conditions, visibility and error precedence survive', {
 # Fast local condition contract; MPI/public ownership is checked separately.
 h <- function(n) getFromNamespace(n, 'npRmpi')
 scope <- h('.npRmpi_with_deferred_package_warnings')
 activate <- h('.npRmpi_warning_scope_activate')
 defer <- h('.npRmpi_defer_named_warning')
 old <- options(warn=0);on.exit(options(old))
 w1 <- structure(list(message='first',call=quote(original(x)),detail=7L),class=c('named_warning','warning','condition'))
 w2 <- simpleWarning('second',quote(second(x)))
 err <- structure(list(message='primary',call=quote(failing(x)),detail=9L),class=c('original_error','error','condition'))
 stages <- character(); caught <- list()
 answer <- withCallingHandlers(scope({activate();defer(warning(w1));defer(warning(w2));stages<-c(stages,'done');invisible(42L)}), warning=function(w){expect_identical(stages,'done');caught[[length(caught)+1L]]<<-w;invokeRestart('muffleWarning')})
 expect_identical(answer,42L);expect_identical(caught,list(w1,w2))
 expect_false(withVisible(scope(invisible(42L)))$visible)
 got <- withCallingHandlers(tryCatch(scope({activate();defer(warning(w1));stop(err)}),error=identity),warning=function(w)stop('must not signal deferred warning'))
 expect_identical(class(got),class(err));expect_identical(conditionMessage(got),conditionMessage(err));expect_identical(conditionCall(got),conditionCall(err));expect_identical(got$detail,9L)
 expect_identical(attr(got,'.npRmpi.deferred.warnings'),list(w1))
 attr(got,'.npRmpi.deferred.warnings')<-NULL;expect_identical(got,err)
 # User callbacks are not intercepted by a broad warning handler.
 before<-FALSE
 got<-tryCatch(scope({activate();warning('unwired');before<-TRUE}),warning=identity)
 expect_s3_class(got,'warning');expect_false(before)
 for(mode in c(-1L,0L,1L)) {
  options(warn=mode);finished<-FALSE
  got<-tryCatch(scope({activate();defer(warning(w1));finished<-TRUE}),warning=identity)
  expect_true(finished);expect_identical(got,w1)
 }
 options(warn=2);finished<-FALSE
 got<-tryCatch(scope({activate();defer(warning('convert'));finished<-TRUE}),error=identity)
 expect_s3_class(got,'error');expect_false(finished);expect_null(attr(got,'.npRmpi.deferred.warnings'))
 options(warn=0)
 # An exiting first handler and a nested scope cannot contaminate another call.
 expect_identical(tryCatch(scope({activate();defer(warning(w1));defer(warning(w2))}),warning=identity),w1)
 expect_identical(scope(88L),88L)
 captured<-list()
 withCallingHandlers(scope({activate();defer(warning(w1));scope({activate();defer(warning(w2))})}),warning=function(w){captured[[length(captured)+1L]]<<-w;invokeRestart('muffleWarning')})
 expect_identical(captured,list(w2,w1))
})
