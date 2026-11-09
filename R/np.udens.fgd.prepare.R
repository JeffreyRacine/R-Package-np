# Gaussian operator adapter. Guide fitting, influence algebra and certificates
# retain the qualified common owners; no namespace rebinding is used.
.np_udens_fgd_write <- function(path,value) {
  connection<-file(path,'wb')
  tryCatch(writeBin(as.vector(value),connection,size=8L,endian='little'),finally=close(connection))
  invisible(NULL)
}

.np_udens_fgd_read <- function(path,nrow,ncol) {
  count<-as.double(nrow)*ncol
  if(!is.finite(count)||count*8>192*1024^2||file.info(path)$size!=count*8)
    stop('invalid Gaussian plot matrix store',call.=FALSE)
  connection<-file(path,'rb');on.exit(close(connection),add=TRUE)
  value<-readBin(connection,'double',count,size=8L,endian='little')
  if(length(value)!=count||any(!is.finite(value)))stop('incomplete Gaussian plot matrix store',call.=FALSE)
  matrix(value,nrow,ncol)
}

.np_udens_fgd_prepare <- function(x,h,grid,full,store=NULL,progress=function(...)NULL,B=9999L,workers=0L) {
  q<-ncol(x);n<-nrow(x);m<-nrow(grid)
  progress('guide',0L,1L)
  fit<-if(q==1L).np_udens_beta_q1_fit(x[,1]) else .np_udens_beta_q2_fit(x)
  theta<-fit$theta;progress('guide',1L,1L)
  # Reject oversized families before materializing the family quadrature plan.
  if(q==1L).np_udens_fgd_workspace(n,grid,theta,full,B,workers)
  variance.plan<-if(q==1L).np_udens_q1_operator_plan(theta,h,grid[,1L]) else NULL
  workspace<-.np_udens_fgd_workspace(n,grid,theta,full,B,workers,
    variance.nodes=if(q==1L)variance.plan$refined.nodes else 0L)
  out<-list(theta=theta,history=fit$history,center=numeric(m),fitted=numeric(m),mu=numeric(m),
    pilot_variance=numeric(m),empirical_variance=numeric(m),guide_variance=numeric(m),
    B=matrix(0,m,length(theta)),B_refined=matrix(0,m,length(theta)),
    observed_files=list(),checks=list(workspace=workspace),covariance_files=list())
  dense<-full&&q==1L
  blocks<-.np_udens_beta_blocks(m,if(dense)m else 64L)
  if(q==2L) {
    main<-.np_udens_beta_q2_moment_prepare(x,h,grid,theta,24L,32L,fit$state,.np_udens_fgd_axis_cache)
    refined<-.np_udens_beta_q2_moment_prepare(x,h,grid,theta,48L,64L,fit$state,.np_udens_fgd_axis_cache)
    if(full)empirical.files<-list()
  }
  for(j in seq_along(blocks)) {
    ix<-blocks[[j]];progress('moments/refinement',min(ix)-1L,m)
    if(q==1L) {
      a<-.np_udens_beta_q1_moments(x[,1],h,grid[ix,1],theta,24L,128L,full,.np_udens_fgd_score,variance.plan)
      b<-.np_udens_beta_q1_moments(x[,1],h,grid[ix,1],theta,48L,256L,full,.np_udens_fgd_score,variance.plan)
    } else {
      a<-.np_udens_beta_q2_moment_block(main,ix,FALSE)
      b<-.np_udens_beta_q2_moment_block(refined,ix,FALSE)
    }
    fields<-c('center','fitted','mu','B','observed','pilot_variance','empirical_variance','guide_variance','Veta')
    if(dense)fields<-c(fields,'covariance','empirical','guide_covariance')
    checks<-lapply(fields,function(name).np_udens_beta_check(a[[name]],b[[name]],
      paste('Gaussian original',name),if(name=='mu')1e-10 else 1e-9));names(checks)<-fields
    checks$observed_mean<-.np_udens_beta_check(rowMeans(a$observed),rowMeans(b$observed),'Gaussian observed mean',1e-10)
    out$checks[[j]]<-checks
    for(name in c('center','fitted','mu','pilot_variance','empirical_variance','guide_variance'))out[[name]][ix]<-a[[name]]
    out$B[ix,]<-a$B;out$B_refined[ix,]<-b$B;out$Veta<-a$Veta
    if(dense)for(name in c('covariance','empirical')) {
      path<-.np_udens_beta_store_file(store,paste0('fgd-',name,'.bin'),8*as.double(m)^2)
      .np_udens_fgd_write(path,a[[name]]);out$covariance_files[[name]]<-path
    }
    if(full&&q==2L) {
      paths<-character(2L)
      for(k in 1:2) {
        z<-if(k==1L)a else b;prep<-if(k==1L)main else refined
        residual<-z$observed-rowMeans(z$observed)-z$B%*%prep$Iinv%*%t(sweep(prep$Td,2,colMeans(prep$Td)))
        paths[k]<-.np_udens_beta_store_file(store,paste0('fgd-influence-',k,'-',j,'.bin'),8*as.double(length(ix))*n)
        .np_udens_fgd_write(paths[k],residual)
      }
      empirical.files[[j]]<-list(rows=ix,paths=paths)
      rm(z,prep,residual)
    }
    if(!is.null(store)&&!full) {
      path<-.np_udens_beta_store_file(store,paste0('observed-',j,'.bin'),8*as.double(length(ix))*n)
      .np_udens_fgd_write(path,a$observed)
      out$observed_files[[j]]<-list(path=path,rows=ix,n=n)
    }
    rm(a,b);progress('moments/refinement',max(ix),m)
  }
  if(full&&q==2L) {
    # Reuse one declared matrix slot. Guide covariance is certified, then the
    # same slot retains pilot covariance. No full matrix enters region fanout.
    path<-.np_udens_beta_store_file(store,'fgd-covariance.bin',8*as.double(m)^2)
    for(law in 0:1) {
      progress(paste(c('guide','pilot')[law+1L],'covariance'),0L,2L)
      a<-.np_udens_beta_covariance_native(main,out$B,law,progress)
      .np_udens_beta_check(diag(a),if(law==0L)out$guide_variance else out$pilot_variance,'Gaussian full/diagonal covariance')
      .np_udens_positive_gram_check(a)
      .np_udens_fgd_write(path,a);rm(a);gc(FALSE)
      b<-.np_udens_beta_covariance_native(refined,out$B_refined,law,progress)
      .np_udens_positive_gram_check(b)
      a<-.np_udens_fgd_read(path,m,m)
      out$checks[[paste0('covariance_law_',law)]]<-.np_udens_beta_dense_compare(a,b,'Gaussian full residual covariance')
      rm(a,b);gc(FALSE)
      progress(paste(c('guide','pilot')[law+1L],'covariance'),2L,2L)
    }
    out$covariance_files$covariance<-path
    rm(main,refined,fit);gc(FALSE)
    # Empirical features are read only after the quadrature caches are retired.
    path<-.np_udens_beta_store_file(store,'fgd-empirical.bin',8*as.double(m)^2)
    for(k in 1:2) {
      residual<-matrix(0,m,n)
      for(block in empirical.files)residual[block$rows,]<-.np_udens_fgd_read(block$paths[k],length(block$rows),n)
      covariance<-tcrossprod(residual)/(n*(n-1));rm(residual);gc(FALSE)
      .np_udens_positive_gram_check(covariance)
      if(k==1L) {
        .np_udens_beta_check(diag(covariance),out$empirical_variance,'Gaussian empirical full/diagonal covariance')
        .np_udens_fgd_write(path,covariance)
      } else {
        a<-.np_udens_fgd_read(path,m,m)
        out$checks$empirical_covariance<-.np_udens_beta_dense_compare(a,covariance,'Gaussian empirical covariance refinement')
        rm(a)
      }
      rm(covariance);gc(FALSE)
    }
    out$covariance_files$empirical<-path
  }
  out
}
