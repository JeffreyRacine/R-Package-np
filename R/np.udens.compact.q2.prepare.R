# Plot-only common-measure compact q2 preparation. No public admission here.
# Native quadrature leaves and the normal owner are reused without modification.
.np_udens_compact_q2_axes <- function(theta,grid,h,kernel,order) {
  lapply(1:2,function(j).np_udens_compact_rule(theta[(2*j-1):(2*j)],order,
    .np_udens_compact_operator_plan(sort(unique(grid[,j])),h[j],kernel)))
}
.np_udens_compact_q2_parts <- function(kernel) .np_udens_compact_normal_q2_axis(kernel)
.np_udens_compact_q2_tile <- function(x,h,grid,kernel,fit,axes,ids) {
  theta<-fit$theta;n<-nrow(x);ax<-lapply(1:2,function(j){z<-axes[[j]];for(k in c('lv','l1v','prob'))z[[k]]<-z[[k]][ids[[j]]];z$nodes<-length(ids[[j]]);z})
  Td<-.np_udens_beta_basis(x);lg<-drop(Td%*%theta)-fit$state$A;rho<-n^(-1/3)
  mats<-vector('list',2L);shift<-numeric(2L)
  for(j in 1:2){v<-exp(ax[[j]]$lv);a<-v/rho;b<-(1-v)/rho
    L<-outer(a,log(x[,j]))+outer(b,log1p(-x[,j]))-lbeta(a+1,b+1)
    L<-sweep(L,2,.5*lg);shift[j]<-max(L);mats[[j]]<-exp(L-shift[j])}
  ratio<-tcrossprod(mats[[1]],mats[[2]])/n
  if(any(!is.finite(ratio))||any(ratio<=0))stop('compact q2 ratio cache range')
  ratio<-log(ratio)+sum(shift);axis<-.np_udens_compact_q2_parts(kernel)
  list(x=x,h=h,grid=grid,theta=theta,state=fit$state,axes=ax,cache=axis(grid,ax,h),ratio=ratio,
    Td=Td,Iinv=solve(fit$state$I),axis.cache=axis,logz=fit$state$A-axes[[1]]$A-axes[[2]]$A)
}
.np_udens_compact_q2_raw_merge <- function(a,b,m) {
  if(is.null(a))return(b)
  peak<-max(a$logmass,b$logmass);lm<-peak+log(exp(a$logmass-peak)+exp(b$logmass-peak))
  w<-exp(b$logmass-lm);v<-exp(a$logmass-lm);d<-b$mean-a$mean
  list(mean=v*a$mean+w*b$mean,variance=v*a$variance+w*b$variance+v*w*d*d,
    cross=v*a$cross+w*b$cross+v*w*outer(d,d[m+seq_len(5L)]),cov=NULL,logmass=lm)
}
.np_udens_compact_q2_stream_moments <- function(x,h,grid,kernel,fit,order,tile=1280L,progress=function(...)NULL,physical=NULL) {
  m<-nrow(grid);n<-nrow(x);axes<-.np_udens_compact_q2_axes(fit$theta,grid,h,kernel,order)
  blocks<-lapply(axes,function(a).np_udens_beta_blocks(a$nodes,tile))
  raw<-list(NULL,NULL);saved<-list();done<-0L;total<-length(blocks[[1]])*length(blocks[[2]])
  for(j in seq_along(blocks[[2]]))for(i in seq_along(blocks[[1]])){
    ids<-list(blocks[[1]][[i]],blocks[[2]][[j]])
    p<-.np_udens_compact_q2_tile(x,h,grid,kernel,fit,axes,ids)
    r<-.np_udens_beta_joint_raw_native(p$cache,p$axes,p$theta[5],p$logz,p$ratio,FALSE)
    for(law in 1:2)raw[[law]]<-.np_udens_compact_q2_raw_merge(raw[[law]],r[[law]],m)
    saved[[length(saved)+1L]]<-list(ids=ids,mass=vapply(r,`[[`,numeric(1),'logmass'),mean=lapply(r,`[[`,'mean'))
    done<-done+1L;progress('joint moments',done,total);rm(p,r)
  }
  score<-m+seq_len(5L);ix<-seq_len(m);g<-raw[[1]];st<-fit$state;theta<-fit$theta
  if(abs(exp(g$logmass)-1)>1e-11)stop('compact q2 guide mass')
  .np_udens_beta_check(g$mean[score],st$mean,'stream score mean',1e-10)
  .np_udens_beta_check(g$cross[score,,drop=FALSE],st$I,'stream information')
  mu<-g$mean[ix];B<-g$cross[ix,,drop=FALSE];K<-B%*%solve(st$I)
  da<-list(list(lv=log(x[,1]),l1v=log1p(-x[,1])),list(lv=log(x[,2]),l1v=log1p(-x[,2])))
  obs<-.np_udens_beta_joint_block_native(.np_udens_compact_normal_q2_axis(kernel,physical)(grid,da,h),seq_len(n),seq_len(n))
  Td<-.np_udens_beta_basis(x);Td<-sweep(Td,2,colMeans(Td));residual<-obs-rowMeans(obs)-K%*%t(Td)
  ve<-rowSums(residual^2)/(n*(n-1))
  variance<-lapply(raw,function(r)(r$variance[ix]-2*rowSums(K*r$cross[ix,,drop=FALSE])+rowSums((K%*%r$cross[score,,drop=FALSE])*K))/n)
  if(any(!is.finite(unlist(variance)))||any(unlist(variance)<=0)||any(ve<=0))stop('compact q2 nonpositive variance')
  U<-Td%*%solve(st$I);Ueta<-sweep(U,2,c(1/(theta[1:4]+1),1),'*')
  f<-exp(drop(.np_udens_beta_basis(grid)%*%theta)-st$A)
  list(center=rowMeans(obs)-mu+f,fitted=f,mu=mu,B=B,observed=obs,
    pilot_variance=variance[[2L]],guide_variance=variance[[1L]],empirical_variance=ve,
    Veta=crossprod(Ueta)/(n*(n-1)),raw=raw,nodes=vapply(axes,`[[`,integer(1),'nodes'),
    x=x,h=h,grid=grid,kernel=kernel,fit=fit,axes=axes,saved=saved)
}

# Separate complete covariance stage; no second law or empirical dense array retained.
.np_udens_compact_q2_covariance <- function(prepared, progress=function(...)NULL) {
  z<-prepared;m<-nrow(z$grid);n<-nrow(z$x);score<-m+seq_len(5L);ix<-seq_len(m)
  K<-z$B%*%solve(z$fit$state$I);C<-NULL;lm<--Inf;mean<-NULL
  for(k in seq_along(z$saved)) {
    s<-z$saved[[k]];p<-.np_udens_compact_q2_tile(z$x,z$h,z$grid,z$kernel,z$fit,z$axes,s$ids)
    cv<-.np_udens_beta_covariance_native(p,z$B,1L,progress)
    mb<-s$mean[[2L]][ix]-drop(K%*%s$mean[[2L]][score]);lb<-s$mass[2L]
    rm(p);gc(FALSE)
    if(is.null(C)){C<-cv;lm<-lb;mean<-mb}else{
      peak<-max(lm,lb);next.mass<-peak+log(exp(lm-peak)+exp(lb-peak));w<-exp(lb-next.mass);v<-exp(lm-next.mass);delta<-mb-mean
      for(rows in .np_udens_beta_blocks(m,64L))C[rows,]<-v*C[rows,,drop=FALSE]+w*cv[rows,,drop=FALSE]+v*w*outer(delta[rows],delta)/n
      mean<-v*mean+w*mb;lm<-next.mass
    }
    rm(cv);gc(FALSE);progress('common-measure covariance merge',k,length(z$saved))
  }
  .np_udens_beta_check(diag(C),z$pilot_variance,'streamed full covariance diagonal')
  .np_udens_positive_gram_check(C);C
}
.np_udens_compact_q2_work_forecast <- function(n,m,nodes,unique.axes,held=0,full=TRUE) {
  if(any(nodes>8192)||prod(as.double(nodes))>16777216)stop('compact q2 common-measure work bound')
  t<-pmin(nodes,1280);axis<-24*sum(unique.axes*t);joint<-prod(as.double(t))
  stages<-8*c(builder=held/8+8*n*max(t)+axis+3*joint,
    covariance=held/8+2*m*m+2*m*256+axis+joint+10*m,
    merge=held/8+3*m*m+8*m*64,
    factor=held/8+2*m*m+8*m*64)
  if(!full)stages<-stages["builder"]
  bytes<-max(stages)+16*1024^2
  if(!is.finite(bytes)||bytes>192*1024^2)stop('compact q2 complete workspace budget')
  list(bytes=bytes,stages=stages,cap=192*1024^2,joint=prod(as.double(nodes)))
}

.np_udens_compact_q2_validate <- function(x,h,grid,kernel) {
  if(!is.matrix(x)||!is.numeric(x)||ncol(x)!=2L||nrow(x)<10L||
     !is.matrix(grid)||!is.numeric(grid)||ncol(grid)!=2L||!nrow(grid)||
     !is.numeric(h)||length(h)!=2L||any(!is.finite(c(x,h,grid)))||
     any(x<=0|x>=1)||any(grid<=0|grid>=1)||any(h<=0)||
     length(kernel)!=1L||!kernel%in%c('uniform','epanechnikov'))
    stop('invalid compact calibrated q2 request',call.=FALSE)
  if(8*(6*as.double(nrow(x))*nrow(grid)+100*nrow(grid))+32*1024^2>192*1024^2)
    stop('compact q2 observed-feature workspace budget',call.=FALSE)
  invisible(NULL)
}
.np_udens_compact_q2_prepare <- function(x,h,grid,kernel,full,store,progress=function(...)NULL,physical=NULL) {
  .np_udens_compact_q2_validate(x,h,grid,kernel)
  n<-nrow(x);m<-nrow(grid)
  progress('original joint guide',0L,1L);fit<-.np_udens_beta_q2_fit(x);progress('original joint guide',1L,1L)
  axes<-lapply(c(16L,32L),function(o).np_udens_compact_q2_axes(fit$theta,grid,h,kernel,o))
  held<-8*(2*as.double(n)*m+100*m)
  budget<-lapply(axes,function(ax).np_udens_compact_q2_work_forecast(n,m,
    vapply(ax,`[[`,integer(1),'nodes'),vapply(1:2,function(j)length(unique(grid[,j])),integer(1)),held,full))
  # The complete common quadrature is independent of covariance/storage tiles.
  a<-.np_udens_compact_q2_stream_moments(x,h,grid,kernel,fit,16L,1280L,progress,physical)
  b<-.np_udens_compact_q2_stream_moments(x,h,grid,kernel,fit,32L,1280L,progress,physical)
  keys<-c('center','fitted','mu','B','observed','pilot_variance','empirical_variance','guide_variance','Veta')
  checks<-setNames(lapply(keys,function(k).np_udens_beta_check(a[[k]],b[[k]],
    paste('compact q2 original',k),if(k=='mu')1e-10 else 1e-9)),keys)
  b$observed<-NULL;gc(FALSE)
  if(full) {
    C<-.np_udens_compact_q2_covariance(a,progress)
    path<-.np_udens_beta_store_file(store,'compact-covariance.bin',8*as.double(m)*m)
    con<-file(path,'wb');tryCatch(writeBin(as.vector(C),con,size=8L,endian='little'),finally=close(con));rm(C);gc(FALSE)
    C<-.np_udens_compact_q2_covariance(b,progress)
    con<-file(path,'rb');main<-tryCatch(matrix(readBin(con,'double',n=as.double(m)*m,size=8L,endian='little'),m),finally=close(con))
    checks$full_covariance<-.np_udens_beta_dense_compare(main,C,'compact q2 full covariance')
    rm(C,main);gc(FALSE)
    a$covariance_file<-path
  }
  rm(b);gc(FALSE)
  a$theta<-fit$theta;a$checks<-c(checks,list(workspace=budget));a
}
.np_udens_compact_q2_covariance_read <- function(prepared) {
  m<-nrow(prepared$grid);path<-prepared$covariance_file
  if(is.null(path)||!file.exists(path)||file.info(path)$size!=8*as.double(m)*m)
    stop('incomplete compact q2 covariance store',call.=FALSE)
  con<-file(path,'rb');on.exit(close(con),add=TRUE)
  C<-matrix(readBin(con,'double',n=as.double(m)*m,size=8L,endian='little'),m)
  if(any(!is.finite(C)))stop('invalid compact q2 covariance store',call.=FALSE)
  C
}
