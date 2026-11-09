library(npRmpi);options(np.messages=FALSE);ns<-asNamespace('npRmpi')
mode<-Sys.getenv('R55_IDENTITY_MODE','reuse');mutant<-identical(Sys.getenv('R55_IDENTITY_MUTANT',''),'1')
if(mode=='reuse'){
 npRmpi.init(nslaves=1L,quiet=TRUE)
 npRmpi:::.npRmpi_bcast_cmd_expr(quote(if(mpi.comm.rank()>0L){e<-get('.npRmpi_worker_identity_state',asNamespace('npRmpi'));e$record$source$r<-paste(rep('0',32),collapse='')}),caller.execute=TRUE)
}else{
 if(mpi.comm.rank(0L)>0L){e<-get('.npRmpi_worker_identity_state',ns);e$record$source$r<-paste(rep('0',32),collapse='')}
}
master<-if(mode=='reuse')TRUE else mpi.comm.rank(0L)==0L
if(master&&mutant)utils::assignInNamespace('.npRmpi_worker_identity_admit',function(comm)invisible(TRUE),ns='npRmpi')
r<-tryCatch(npRmpi.init(mode=if(mode=='reuse')'spawn'else'attach',nslaves=1L,quiet=TRUE),error=identity)
if(master){
 detected<-inherits(r,'error')&&grepl('worker code identity mismatch',conditionMessage(r),fixed=TRUE)
 cat('IDENTITY_ADMISSION_RESULT',mode,'detected',detected,'mutant',mutant,'\n');flush.console()
 npRmpi.quit(force=TRUE,mode=if(mode=='reuse')'spawn'else'attach')
 if(!mutant)stopifnot(detected)else stopifnot(!detected)
 cat(if(mutant)'IDENTITY_NOOP_MUTANT_KILLED'else'IDENTITY_ADMISSION_PASS','\n');flush.console()
 mpi.quit()
}else mpi.quit()
