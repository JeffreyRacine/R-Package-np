# Invocation-owned files and explicit private R streams. The sampler child loads
# stats only; it never evaluates caller code, loads np/npRmpi or initializes MPI.
.np_udens_beta_store <- function() {
  path<-tempfile('np-beta-plot-')
  if(!dir.create(path,mode='0700'))stop('cannot create beta plot workspace',call.=FALSE)
  store<-new.env(parent=emptyenv());store$path<-path;store$closed<-FALSE
  store$budget<-512*1024^2;store$declared<-0;store$files<-list()
  store
}

.np_udens_beta_store_file <- function(store,name,bytes) {
  if(store$closed||!is.character(name)||length(name)!=1L||!grepl('^[a-z][a-z0-9_.-]*$',name)||
     name%in%names(store$files)||!is.numeric(bytes)||length(bytes)!=1L||!is.finite(bytes)||
     bytes<0||bytes!=floor(bytes)||store$declared+bytes>store$budget)
    stop('beta plot owned-file budget/name',call.=FALSE)
  path<-file.path(store$path,name)
  if(file.exists(path))stop('beta plot owned-file collision',call.=FALSE)
  store$declared<-store$declared+bytes;store$files[[name]]<-list(path=path,bytes=bytes)
  path
}

.np_udens_beta_store_close <- function(store) {
  if(!store$closed) {
    store$closed<-TRUE
    status<-unlink(store$path,recursive=TRUE)
    if(status!=0L||dir.exists(store$path))stop('cannot retire beta plot workspace',call.=FALSE)
  }
  invisible(NULL)
}

.np_udens_beta_random_request <- function(mode,n,B,seed,kind) {
  uniform<-c('Wichmann-Hill','Marsaglia-Multicarry','Super-Duper','Mersenne-Twister',
    'Knuth-TAOCP','Knuth-TAOCP-2002',"L'Ecuyer-CMRG")
  normal<-c('Buggy Kinderman-Ramage','Ahrens-Dieter','Box-Muller','Inversion','Kinderman-Ramage')
  if(!is.character(mode)||length(mode)!=1L||!mode%in%c('normal','multinomial')||
     !is.numeric(n)||length(n)!=1L||!is.finite(n)||n<1||n!=floor(n)||n>.Machine$integer.max||
     !is.numeric(B)||length(B)!=1L||!is.finite(B)||B<1||B!=floor(B)||B>.Machine$integer.max||
     !is.numeric(seed)||length(seed)!=1L||!is.finite(seed)||seed<0||seed!=floor(seed)||seed>.Machine$integer.max||
     !is.character(kind)||length(kind)!=3L||anyNA(kind)||!kind[1]%in%uniform||!kind[2]%in%normal||
     !kind[3]%in%c('Rounding','Rejection'))stop('invalid beta plot random stream request',call.=FALSE)
  bytes<-as.double(n)*as.double(B)*if(mode=='normal')8 else 4
  if(!is.finite(bytes)||bytes>512*1024^2)stop('beta plot random stream exceeds owned-disk budget',call.=FALSE)
  list(mode=mode,n=as.integer(n),B=as.integer(B),seed=as.integer(seed),kind=kind,
       bytes=bytes,version=R.version$version.string)
}

.np_udens_beta_sampler_child <- function() {
  args<-commandArgs(TRUE);request<-readRDS(args[1])
  stopifnot(identical(request$version,R.version$version.string))
  do.call(RNGkind,as.list(request$kind));set.seed(request$seed);entry<-.Random.seed
  connection<-file(args[2],'wb');on.exit(close(connection),add=TRUE)
  batch<-max(1L,min(256L,floor(8*1024^2/(8*as.double(request$n)))))
  for(first in seq.int(1L,request$B,by=batch)) {
    count<-min(batch,request$B-first+1L)
    if(request$mode=='normal')writeBin(stats::rnorm(as.double(request$n)*count),connection,size=8L,endian='little')
    else writeBin(as.vector(stats::rmultinom(count,request$n,rep(1/request$n,request$n))),connection,size=4L,endian='little')
  }
  flush(connection)
  saveRDS(list(request=request,entry=entry,end=.Random.seed,kind=RNGkind(),pid=Sys.getpid(),
    namespaces=loadedNamespaces()),args[3])
}

.np_udens_beta_random_stream <- function(store,mode,n,B,seed,kind=RNGkind()) {
  request<-.np_udens_beta_random_request(mode,n,B,seed,kind)
  output<-.np_udens_beta_store_file(store,'random.bin',request$bytes)
  # Charge a conservative bounded allowance for the request, script and receipt.
  script<-.np_udens_beta_store_file(store,'sampler.r',65536)
  input<-.np_udens_beta_store_file(store,'request.rds',65536)
  receipt<-.np_udens_beta_store_file(store,'sampler-receipt.rds',65536)
  log<-.np_udens_beta_store_file(store,'sampler.log',65536)
  writeLines(c('(',deparse(.np_udens_beta_sampler_child,width.cutoff=500L),')()'),script)
  saveRDS(request,input)
  executable<-file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
  if(!file.exists(executable))stop('Rscript is unavailable for beta plot draws',call.=FALSE)
  raw<-system2(executable,c('--vanilla','--default-packages=stats',shQuote(c(script,input,output,receipt))),
    stdout=log,stderr=log,timeout=600L)
  if(!identical(as.integer(raw),0L))stop(sprintf('beta plot random sampler failed (status %s)',raw),call.=FALSE)
  if(!file.exists(output)||!file.exists(receipt)||file.info(output)$size!=request$bytes)
    stop('incomplete beta plot random stream',call.=FALSE)
  for(f in store$files)if(file.exists(f$path)&&file.info(f$path)$size>f$bytes)
    stop('beta plot random stream exceeded its declared storage',call.=FALSE)
  info<-readRDS(receipt)
  if(!identical(info$request,request)||!identical(info$kind,kind)||!is.integer(info$pid)||
     length(info$pid)!=1L||info$pid==Sys.getpid()||any(c('np','npRmpi','Rmpi')%in%info$namespaces))
    stop('invalid beta plot random stream receipt',call.=FALSE)
  list(path=output,request=request,receipt=info)
}

.np_udens_beta_random_read <- function(stream,first,last) {
  request<-stream$request
  if(length(first)!=1L||length(last)!=1L||!is.finite(first)||!is.finite(last)||
     first!=floor(first)||last!=floor(last)||first<1L||last<first||last>request$B)
    stop('invalid beta plot random stream range',call.=FALSE)
  count<-as.double(request$n)*(last-first+1)
  if(count*8>32*1024^2)stop('beta plot random read exceeds batch budget',call.=FALSE)
  connection<-file(stream$path,'rb');on.exit(close(connection),add=TRUE)
  size<-if(request$mode=='normal')8L else 4L
  seek(connection,where=(first-1)*as.double(request$n)*size,origin='start')
  values<-readBin(connection,if(request$mode=='normal')'double' else 'integer',n=count,size=size,endian='little')
  if(length(values)!=count||any(!is.finite(values)))stop('incomplete beta plot random batch',call.=FALSE)
  result<-matrix(as.double(values),request$n)
  if(request$mode=='multinomial'&&(any(values<0L)||any(colSums(result)!=request$n)))
    stop('invalid beta plot multinomial batch',call.=FALSE)
  result
}
