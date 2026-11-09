# Hidden normalized-support providers. Public route admission is a separate seam.
.np_udens_beta_prepare <- function(x,h,grid,full,store=NULL,progress=function(...)NULL,B=9999L,workers=0L) {
  q<-ncol(x);n<-nrow(x);m<-nrow(grid)
  progress('guide',0L,1L)
  fit<-if(q==1L).np_udens_beta_q1_fit(x[,1]) else .np_udens_beta_q2_fit(x)
  theta<-fit$theta;progress('guide',1L,1L)
  # Reject oversized families before materializing the family quadrature plan.
  if(q==1L).np_udens_beta_workspace(n,grid,theta,full,B,workers)
  variance.plan<-if(q==1L).np_udens_q1_variance_plan(theta,h,grid[,1L]) else NULL
  workspace<-.np_udens_beta_workspace(n,grid,theta,full,B,workers,
    variance.nodes=if(q==1L)variance.plan$nodes else 0L)
  out<-list(theta=theta,history=fit$history,center=numeric(m),fitted=numeric(m),mu=numeric(m),
    pilot_variance=numeric(m),empirical_variance=numeric(m),guide_variance=numeric(m),
    B=matrix(0,m,length(theta)),B_refined=matrix(0,m,length(theta)),
    observed_files=list(),checks=list(workspace=workspace))
  dense.block<-full&&q==1L
  blocks<-.np_udens_beta_blocks(m,if(dense.block)m else 64L)
  if(dense.block) {
    nodes<-max(.np_udens_beta_marginal_rule(theta,256L)$nodes,variance.plan$nodes)
    if(8*(10*as.double(m)*nodes+5*as.double(m)^2+3*as.double(m)*n+10*nodes)>192*1024^2)
      stop('beta plot univariate dense workspace budget',call.=FALSE)
  }
  if(q==2L) {
    main<-.np_udens_beta_q2_moment_prepare(x,h,grid,theta,64L,32L,fit$state)
    refined<-.np_udens_beta_q2_moment_prepare(x,h,grid,theta,128L,64L,fit$state)
  }
  for(j in seq_along(blocks)) {
    ix<-blocks[[j]]
    progress('moments/refinement',min(ix)-1L,m)
    if(q==1L) {
      a<-.np_udens_beta_q1_moments(x[,1],h,grid[ix,1],theta,64L,128L,full,variance.plan=variance.plan)
      b<-.np_udens_beta_q1_moments(x[,1],h,grid[ix,1],theta,128L,256L,full)
    } else {
      a<-.np_udens_beta_q2_moment_block(main,ix,FALSE)
      b<-.np_udens_beta_q2_moment_block(refined,ix,FALSE)
    }
    names<-c('center','fitted','mu','B','observed','pilot_variance','empirical_variance','Veta')
    names<-c(names,'guide_variance')
    if(dense.block)names<-c(names,'covariance','empirical','guide_covariance')
    checks<-lapply(names,function(name).np_udens_beta_check(a[[name]],b[[name]],
      paste('original',name),if(name=='mu')1e-10 else 1e-9));names(checks)<-names
    checks$observed_mean<-.np_udens_beta_check(rowMeans(a$observed),rowMeans(b$observed),'original observed mean',1e-10)
    out$checks[[j]]<-checks
    for(name in c('center','fitted','mu','pilot_variance','empirical_variance','guide_variance'))out[[name]][ix]<-a[[name]]
    out$B[ix,]<-a$B;out$B_refined[ix,]<-b$B;out$Veta<-a$Veta
    if(dense.block)out$covariance<-a$covariance
    if(!is.null(store)&&!full) {
      path<-.np_udens_beta_store_file(store,paste0('observed-',j,'.bin'),8*as.double(length(ix))*n)
      connection<-file(path,'wb')
      tryCatch(writeBin(as.vector(a$observed),connection,size=8L,endian='little'),finally=close(connection))
      out$observed_files[[j]]<-list(path=path,rows=ix,n=n)
    }
    progress('moments/refinement',max(ix),m)
  }
  if(full&&q==2L) {
    for(law in 0:1) {
      progress(paste(c('guide','ratio')[law+1L],'covariance main'),0L,1L)
      a<-.np_udens_beta_covariance_native(main,out$B,law,progress)
      .np_udens_beta_check(diag(a),if(law==0L)out$guide_variance else out$pilot_variance,
        'full/diagonal residual covariance identity')
      .np_udens_positive_gram_check(a)
      path<-.np_udens_beta_store_file(store,paste0('covariance-',law,'.bin'),8*as.double(m)*m)
      connection<-file(path,'wb')
      tryCatch(writeBin(as.vector(a),connection,size=8L,endian='little'),finally=close(connection))
      rm(a);gc(FALSE)
      progress(paste(c('guide','ratio')[law+1L],'covariance refinement'),0L,1L)
      b<-.np_udens_beta_covariance_native(refined,out$B_refined,law,progress)
      .np_udens_positive_gram_check(b)
      connection<-file(path,'rb')
      a<-tryCatch(matrix(readBin(connection,'double',n=as.double(m)*m,size=8L,endian='little'),m),finally=close(connection))
      out$checks[[paste0('covariance_law_',law)]]<-.np_udens_beta_dense_compare(a,b,'full residual covariance')
      if(law==1L)out$covariance<-a
      rm(a,b);gc(FALSE)
    }
  }
  out
}

.np_udens_beta_observed_read <- function(block) {
  count<-as.double(length(block$rows))*block$n
  if(count*8>32*1024^2||file.info(block$path)$size!=count*8)stop('invalid beta plot observed store',call.=FALSE)
  connection<-file(block$path,'rb');on.exit(close(connection),add=TRUE)
  values<-readBin(connection,'double',n=count,size=8L,endian='little')
  if(length(values)!=count||any(!is.finite(values)))stop('incomplete beta plot observed store',call.=FALSE)
  matrix(values,length(block$rows))
}

.np_udens_beta_regions <- function(theta,Veta,grid,alpha,progress) {
  .np_udens_beta_region_dispatch(theta,Veta,grid,alpha,progress)
}

.np_udens_beta_sd_roots <- function(x,h,grid,prepared,stream,store,progress,
                                   score=.np_udens_beta_q1_score,
                                   axis.cache=.np_udens_beta_axis_cache,
                                   inner=c(64L,128L)) {
  n<-nrow(x);q<-ncol(x);m<-nrow(grid);B<-stream$request$B
  theta<-prepared$theta;p<-length(theta);thetas<-matrix(0,B,p);iterations<-numeric(B)
  draw.blocks<-.np_udens_beta_blocks(B,min(32L,max(1L,floor(8*1024^2/(n*8)))))
  refits<-.np_udens_beta_refit_dispatch(x,theta,stream,progress)
  thetas<-refits$theta;iterations<-refits$iterations
  roots<-.np_udens_beta_root_store(store,m,B)
  if(q==1L) {
    # The common rule uses the same componentwise minimum over ALL accepted fits.
    all.theta<-rbind(theta,thetas)
    rhat<-prepared$center-prepared$fitted
    scale<-sqrt(prepared$empirical_variance)
    for(j in seq_along(roots$rows)) {
      ix<-roots$rows[[j]]
      main<-.np_udens_beta_q1_common(x[,1],h,grid[ix,1],all.theta,inner[1L],128L,score)
      refined<-.np_udens_beta_q1_common(x[,1],h,grid[ix,1],all.theta,inner[2L],256L,score)
      .np_udens_beta_check(main$mu,refined$mu,'bootstrap common guide means',1e-10)
      .np_udens_beta_check(main$fitted,refined$fitted,'bootstrap common guide density',1e-10)
      .np_udens_beta_check(main$mu[1,],prepared$mu[ix],'bootstrap original mean',1e-10)
      rm(refined)
      values<-matrix(0,length(ix),B)
      for(ids in draw.blocks) {
        counts<-.np_udens_beta_random_read(stream,min(ids),max(ids))
        values[,ids]<-(main$observed%*%(counts/n)-t(main$mu[1L+ids,,drop=FALSE])-rhat[ix])/scale[ix]
      }
      if(any(!is.finite(values)))stop('nonfinite beta plot bootstrap roots',call.=FALSE)
      connection<-file(roots$files[[j]],'wb')
      tryCatch(writeBin(as.vector(values),connection,size=8L,endian='little'),finally=close(connection))
      progress('bootstrap projection',max(ix),m)
    }
    roots$next_draw<-B+1L
    return(list(roots=roots,theta=thetas,iterations=iterations))
  } else {
    main<-.np_udens_beta_series_projection(.np_udens_beta_series_rule(theta,32L),h,n,grid,inner[1L],axis.cache)
    refined<-.np_udens_beta_series_projection(.np_udens_beta_series_rule(theta,64L),h,n,grid,inner[2L],axis.cache)
    original<-.np_udens_beta_series_means_native(main,matrix(theta,1L))
    original.ref<-.np_udens_beta_series_means_native(refined,matrix(theta,1L))
    .np_udens_beta_check(original$mu,original.ref$mu,'bootstrap original mean refinement',1e-10)
    .np_udens_beta_check(drop(original$mu),prepared$mu,'bootstrap original mean',1e-10)
  }
  rhat<-prepared$center-prepared$fitted
  scale<-sqrt(prepared$empirical_variance)
  for(ix in draw.blocks) {
    counts<-.np_udens_beta_random_read(stream,min(ix),max(ix));values<-matrix(0,m,length(ix))
    for(block in prepared$observed_files)values[block$rows,]<-.np_udens_beta_observed_read(block)%*%(counts/n)
    projection<-.np_udens_beta_series_means_native(main,thetas[ix,,drop=FALSE])
    projection.ref<-.np_udens_beta_series_means_native(refined,thetas[ix,,drop=FALSE])
    .np_udens_beta_check(projection$mu,projection.ref$mu,'bootstrap guide means',1e-10)
    .np_udens_beta_check(projection$fitted,projection.ref$fitted,'bootstrap guide density',1e-10)
    delta<-values-t(projection$mu)-rhat
    roots<-.np_udens_beta_root_write(roots,delta/scale,min(ix))
    progress('bootstrap projection',max(ix),B)
  }
  list(roots=roots,theta=thetas,iterations=iterations)
}

.np_udens_beta_hidden <- function(x,h,grid,alpha,variant,band,B=9999L,seed=42L,
                                  progress=function(...)NULL) {
  if(!is.matrix(x)||!is.matrix(grid)||!ncol(x)%in%1:2||ncol(x)!=ncol(grid)||nrow(x)<10L||!nrow(grid)||
     any(!is.finite(c(x,grid,h)))||any(x<=0|x>=1)||any(grid<=0|grid>=1)||
     length(h)!=ncol(x)||any(h<=0)||length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1||
     length(variant)!=1L||!variant%in%c('asymptotic','symmetric-sd','asymmetric-sd')||
     length(band)!=1L||!band%in%c('pointwise','simultaneous','bonferroni','all'))stop('invalid beta plot hidden provider request',call.=FALSE)
  .npRmpi_bootstrap_fanout_enabled(n=nrow(x),B=nrow(grid),chunk.size=8L,what='beta density bands')
  # Canonical unique targets make a seeded family independent of caller ordering.
  original<-grid
  unique.grid<-unique(grid);ordering<-do.call(order,as.data.frame(unique.grid))
  grid<-unique.grid[ordering,,drop=FALSE]
  keys<-function(x)apply(x,1,function(z)paste(format(z,digits=17,scientific=TRUE),collapse=':'))
  restore<-match(keys(original),keys(grid))
  sd<-variant!='asymptotic';draws<-sd||band%in%c('simultaneous','all')
  if(draws).np_udens_beta_random_request(if(sd)'multinomial' else 'normal',if(sd)nrow(x) else nrow(grid),B,seed,RNGkind())
  # Preflight the whole owned-disk requirement before starting mathematical work.
  disk<-if(sd)4*as.double(nrow(x))*B+8*as.double(nrow(grid))*B+8*as.double(nrow(grid))*nrow(x)+262144 else
    if(draws)8*as.double(nrow(grid))*B+if(ncol(x)==2L)16*as.double(nrow(grid))^2+262144 else 262144 else 0
  if(disk>512*1024^2)stop('beta plot inference exceeds owned-disk budget',call.=FALSE)
  store<-if(draws).np_udens_beta_store() else NULL
  if(draws)on.exit(.np_udens_beta_store_close(store),add=TRUE)
  prepared<-.np_udens_beta_prepare(x,h,grid,!sd&&band%in%c('simultaneous','all'),store,progress,B,.npRmpi_bootstrap_worker_count())
  region<-.np_udens_beta_regions(prepared$theta,prepared$Veta,grid,alpha,progress)
  scale<-sqrt(prepared$pilot_variance);rhat<-prepared$center-prepared$fitted
  if(draws)stream<-.np_udens_beta_random_stream(store,if(sd)'multinomial' else 'normal',
    if(sd)nrow(x) else nrow(grid),B,seed)
  if(sd) {
    bootstrap<-.np_udens_beta_sd_roots(x,h,grid,prepared,stream,store,progress)
    result<-if(band%in%c('bonferroni','all'))
      .np_udens_bounded_sd(bootstrap$roots,rhat,list(scale),region,prepared$center,alpha,band,variant,progress) else
      .np_udens_beta_stored_reduce(bootstrap$roots,rhat,scale,region,prepared$center,alpha,band,variant)
  } else {
    critical<-if(band=='pointwise')qnorm(1-alpha/4) else if(band=='bonferroni')
      if(nrow(grid)==1L)qnorm(1-alpha/4) else qnorm(alpha/(4*nrow(grid)),lower.tail=FALSE) else
      .np_udens_beta_gaussian_critical(prepared$covariance,B,alpha,stream,progress)
    result<-list(lower=rhat+region$lower-critical*scale,upper=rhat+region$upper+critical*scale,
      critical=critical,B=if(draws)B else 0L)
  }
  if(band=='all'&&!sd) {
    companions<-.np_udens_bounded_analytic(rhat,region,list(scale),alpha,nrow(grid))
    result$bands<-list(pointwise=companions$pointwise,
      simultaneous=list(lower=result$lower,upper=result$upper,critical=result$critical),
      bonferroni=companions$bonferroni)
  }
  result$center<-prepared$center;result$pilot_se<-scale;result$empirical_se<-sqrt(prepared$empirical_variance)
  for(name in c('lower','upper','center','pilot_se','empirical_se'))result[[name]]<-result[[name]][restore]
  if(band%in%c('bonferroni','all'))result$family_size<-nrow(grid)
  if(band=='all')result$bands<-lapply(result$bands,function(x) {
    x$lower<-x$lower[restore];x$upper<-x$upper[restore];x
  })
  result$theta<-prepared$theta;result$variant<-variant;result$band<-band;result$alpha<-alpha
  result$construction<-switch(variant,asymptotic='residual-Gaussian-plus-guide-region',
    'symmetric-sd'='residual-absolute-containing-hull','asymmetric-sd'='signed-residual-inversion')
  result$checks<-prepared$checks
  if(sd&&band%in%c('bonferroni','all')) {
    if(result$sd_companion$resolution.limited)
      warning('beta plot Bonferroni calibration has low finite-draw tail resolution; returning finite empirical limits with the simultaneous hull',call.=FALSE)
  } else if(draws&&(if(sd&&variant=='asymmetric-sd')floor(B*alpha/4) else floor(B*alpha/2))<20)
    warning('beta plot calibration has low finite-draw tail resolution; returning the finite band',call.=FALSE)
  result
}
