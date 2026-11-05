# Bounded external sentinel. Run one block per fresh process/pool.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("usage: block[1:9] nslaves output_directory")
block <- as.integer(args[1L]); nslaves <- as.integer(args[2L])
stopifnot(length(block)==1L,!is.na(block),block %in% 1:9,
          length(nslaves)==1L,!is.na(nslaves),nslaves>=1L)
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
script <- file.path(dirname(normalizePath(self)),"delegated_numerical",sprintf("gap%02d.R",block))
output <- normalizePath(args[3L],mustWork=TRUE)
setwd(output)
pkg <- "npRmpi"
suppressPackageStartupMessages(library(pkg,character.only=TRUE));library(testthat)
options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,rgl.useNULL=TRUE,
 device=function(...) grDevices::pdf(file=tempfile(fileext='.pdf'),...))
Sys.setenv(NP_NUMERICAL_RESULT_DIR=getwd())
closed<-FALSE
main <- function() {
 if(pkg=='npRmpi') {
  npRmpi.init(nslaves=nslaves,quiet=TRUE)
  on.exit(if(!closed)npRmpi.quit(force=TRUE),add=TRUE)
  options(npRmpi.reuse.slaves=FALSE)
 }
 cat('PACKAGE',find.package(pkg),'NATIVE',getLoadedDLLs()[[pkg]][['path']],'\n')
 collective <- function(expr) {
  if(pkg=='np') return(eval(expr,envir=.GlobalEnv))
  npRmpi:::.npRmpi_bcast_cmd_expr(expr,comm=1L,caller.execute=TRUE)
 }
 assign('collective',collective,envir=.GlobalEnv)
 if(pkg=='npRmpi') collective(quote(cat('NUMERICAL_WORKER',npRmpi::mpi.comm.rank(1),'of',npRmpi::mpi.comm.size(1),'\n')))
 testthat::set_max_fails(Inf)
 result<-testthat::test_file(script,reporter='summary',stop_on_failure=FALSE,stop_on_warning=FALSE)
 frame<-as.data.frame(result);saveRDS(frame,'results.rds')
 out<-frame[intersect(c('test','nb','failed','error','skipped','warning','passed','real'),names(frame))]
 out$messages<-lapply(frame$result,function(z)lapply(z,function(e)list(class=class(e)[1],message=if(is.null(e$message))'' else e$message)))
 jsonlite::write_json(out,'assertions.json',pretty=TRUE,auto_unbox=TRUE,dataframe='rows',na='null')
 if(pkg=='npRmpi'){npRmpi.quit(force=TRUE);closed<<-TRUE}
 if(!nrow(frame)||any(frame$failed>0|frame$error|frame$skipped))1L else 0L
}
status<-main();cat('BLOCK_COMPLETE status=',status,'\n');quit(save='no',status=status,runLast=FALSE)
