.np_udens_beta_blocks <- function(n,block=64L) {
  lapply(seq.int(1L,n,by=block),function(first)seq.int(first,min(n,first+block-1L)))
}

.np_udens_beta_root_store <- function(store,m,B) {
  rows<-.np_udens_beta_blocks(m)
  files<-lapply(seq_along(rows),function(j) {
    file<-.np_udens_beta_store_file(store,paste0('roots-',j,'.bin'),8*as.double(length(rows[[j]]))*B)
    connection<-file(file,'wb');close(connection);file
  })
  list(rows=rows,files=files,B=B,m=m,next_draw=1L)
}

.np_udens_beta_root_write <- function(root,values,first) {
  if(!is.matrix(values)||nrow(values)!=root$m||first!=root$next_draw||
     first+ncol(values)-1L>root$B||any(!is.finite(values)))stop('invalid beta plot root batch',call.=FALSE)
  for(j in seq_along(root$rows)) {
    connection<-file(root$files[[j]],'ab')
    tryCatch(writeBin(as.vector(values[root$rows[[j]],,drop=FALSE]),connection,size=8L,endian='little'),
      finally=close(connection))
  }
  root$next_draw<-first+ncol(values);root
}

.np_udens_beta_root_read <- function(root,j) {
  count<-as.double(length(root$rows[[j]]))*root$B
  if(root$next_draw!=root$B+1L||count*8>32*1024^2||file.info(root$files[[j]])$size!=count*8)
    stop('incomplete beta plot root store',call.=FALSE)
  connection<-file(root$files[[j]],'rb');on.exit(close(connection),add=TRUE)
  x<-readBin(connection,'double',n=count,size=8L,endian='little')
  if(length(x)!=count||any(!is.finite(x)))stop('invalid beta plot stored roots',call.=FALSE)
  matrix(x,length(root$rows[[j]]))
}

.np_udens_beta_stored_reduce <- function(root,rhat,scale,region,center,alpha,band,calibration) {
  lower<-upper<-numeric(root$m)
  if(band=='simultaneous') {
    lo<-rep(Inf,root$B);hi<-rep(-Inf,root$B)
    for(j in seq_along(root$rows)) {
      roots<-.np_udens_beta_root_read(root,j)
      lo<-pmin(lo,apply(roots,2,min));hi<-pmax(hi,apply(roots,2,max))
    }
    # Two extreme rows retain exactly the same family extrema and absolute max.
    extremes<-rbind(lo,hi)
    k<-floor(root$B*alpha/4)+1L;ka<-ceiling(root$B*(1-alpha/2))
    if(calibration=='symmetric-sd') {
      critical<-sort(pmax(abs(lo),abs(hi)),partial=ka)[ka];dl<--critical;du<-critical
    } else {
      dl<--sort(hi,partial=root$B+1L-k)[root$B+1L-k]
      du<--sort(lo,partial=k)[k]
    }
    lower<-rhat+region$lower+scale*dl;upper<-rhat+region$upper+scale*du
    if(calibration=='symmetric-sd') {
      hull<-.np_udens_beta_hull(center,lower,upper);lower<-hull[,1];upper<-hull[,2]
    }
  } else {
    for(j in seq_along(root$rows)) {
      ix<-root$rows[[j]]
      out<-.np_udens_beta_sd_reduce(.np_udens_beta_root_read(root,j),rhat[ix],scale[ix],
        list(lower=region$lower[ix],upper=region$upper[ix]),center[ix],alpha,band,calibration)
      lower[ix]<-out$lower;upper[ix]<-out$upper
    }
  }
  if(any(!is.finite(c(lower,upper)))||any(lower>upper))stop('invalid beta plot reduced endpoints',call.=FALSE)
  list(lower=lower,upper=upper,B=root$B,k_signed=floor(root$B*alpha/4)+1L,
    k_absolute=ceiling(root$B*(1-alpha/2)))
}

.np_udens_beta_factor <- function(covariance) {
  if(!is.matrix(covariance)||nrow(covariance)!=ncol(covariance)||any(!is.finite(covariance))||
     any(diag(covariance)<=0))stop('invalid beta plot covariance',call.=FALSE)
  .Call(C_np_plot_udens_beta_factor,covariance)
}

.np_udens_beta_gaussian_critical <- function(covariance,B,alpha,stream,progress) {
  factor<-.np_udens_beta_factor(covariance);maxima<-numeric(B)
  for(ix in .np_udens_beta_blocks(B,64L)) {
    z<-.np_udens_beta_random_read(stream,min(ix),max(ix))
    maxima[ix]<-apply(abs(factor%*%z),2,max)
    progress('Gaussian calibration',max(ix),B)
  }
  if(any(!is.finite(maxima)))stop('invalid beta plot Gaussian maxima',call.=FALSE)
  unname(quantile(maxima,1-alpha/2,type=7))
}
