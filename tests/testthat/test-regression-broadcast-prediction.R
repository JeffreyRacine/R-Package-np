test_that("prediction and residuals retain the established manual broadcast owner", {
 ns<-asNamespace("npRmpi")
 for(pool in c(FALSE,TRUE))for(bcast in c(FALSE,TRUE))for(method in c("predict.npregression","residuals.npregression")) {
   e<-new.env(parent=ns);f<-get(method,ns);environment(f)<-e
   e$.npRmpi_has_active_slave_pool<-function(...)pool
   e$.npRmpi_autodispatch_called_from_bcast<-function(...)bcast
   localized<-0L
   e$.npRmpi_with_local_regression<-function(expr){localized<<-localized+1L;force(expr)}
   e$npreg<-function(...)structure(list(mean=1:3,resid=4:6),class="npregression")
   object<-list(bws=list(formula=quote(y~x)),residuals=FALSE)
   value<-f(object)
   expect_identical(value,if(method=="predict.npregression")1:3 else 4:6)
   expect_identical(localized,as.integer(!(pool||bcast)))
 }
})

test_that("manual broadcast prediction and residuals complete with direct-fit parity", {
 skip_on_cran()
 lines<-c(
  "library(npRmpi)",
  "options(np.messages=FALSE)",
  "npRmpi.init(nslaves=as.integer(Sys.getenv(\"N5_SLAVES\",\"1\")),quiet=TRUE)",
  "on.exit(npRmpi.quit(force=TRUE),add=TRUE)",
  "set.seed(821);dd<-data.frame(x=runif(80));dd$y<-sin(dd$x)+rnorm(80,sd=.2)",
  "nd<-data.frame(x=c(.2,.4,.6))",
  "mpi.bcast.Robj2slave(dd);mpi.bcast.Robj2slave(nd)",
  "mpi.bcast.cmd(bb<-npregbw(y~x,data=dd,bws=.2,bandwidth.compute=FALSE),caller.execute=TRUE)",
  "mpi.bcast.cmd(ff<-npreg(bws=bb,data=dd),caller.execute=TRUE)",
  "mpi.bcast.cmd(ref<-npreg(bws=bb,data=dd,newdata=nd,se=TRUE),caller.execute=TRUE)",
  "cat(\"DIRECT_CONTROL_COMPLETE\\n\");flush(stdout())",
  "mpi.bcast.cmd(pp<-predict(ff,newdata=nd,se.fit=TRUE),caller.execute=TRUE)",
  "stopifnot(identical(pp$fit,fitted(ref)),identical(pp$se.fit,se(ref)))",
  "cat(\"PREDICT_COMPLETE\\n\");flush(stdout())",
  "mpi.bcast.cmd(rr<-residuals(ff),caller.execute=TRUE)",
  "stopifnot(isTRUE(all.equal(rr,dd$y-fitted(ff),tolerance=1e-14)))",
  "cat(\"RESIDUALS_COMPLETE\\n\");flush(stdout())",
  "# Ordinary auto-dispatch remains usable after the manual calls.",
  "au<-predict(ff,newdata=nd,se.fit=TRUE)",
  "stopifnot(isTRUE(all.equal(au$fit,pp$fit,tolerance=1e-14)))",
  "npRmpi.quit(force=TRUE)",
  "cat(\"N5_PREDICTION_BCAST_COMPLETE\\n\")"
 )
 out<-npRmpi_run_rscript_subprocess(lines,timeout=60L)
 expect_identical(out$status,0L)
 expect_true(any(out$output=="N5_PREDICTION_BCAST_COMPLETE"))
})
