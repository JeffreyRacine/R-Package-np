# Common quadrature for the original guide and every literal bootstrap refit.
# All target tiles retain the full support-knot union and one all-fit base.
.np_udens_compact_common <- function(plan,thetas,order,rows=seq_along(plan$target),
                                    held.bytes=0,fit.block=NULL,target.block=NULL) {
  if(!is.matrix(thetas)||ncol(thetas)!=2L||!nrow(thetas)||
     any(!is.finite(thetas))||any(thetas<= -1)||!length(rows)||
     any(!rows%in%seq_along(plan$target)))
    stop('invalid compact bootstrap common projection',call.=FALSE)
  base<-apply(thetas,2,min);st<-.np_udens_compact_rule(base,order,plan)
  nr<-length(rows);nb<-nrow(thetas);nodes<-st$nodes;cap<-256*1024^2
  # Includes the other refinement, output/root copies and the caller's live state.
  reserve<-held.bytes+8*(40*nodes+4*as.double(nb)*nr)+32*1024^2
  available<-cap-reserve
  nf<-min(nb,floor(min(32*1024^2,available/2)/(8*6*nodes)))
  if(!is.finite(nf)||nf<1L)stop('compact bootstrap common workspace exceeds 256 MiB',call.=FALSE)
  if(!is.null(fit.block)){
    if(length(fit.block)!=1L||!is.finite(fit.block)||fit.block<1L||fit.block!=floor(fit.block))
      stop('invalid compact bootstrap fit block',call.=FALSE)
    nf<-min(nf,fit.block)
  }
  nt<-min(nr,floor((available-8*6*nodes*nf)/(8*2*nodes)))
  if(!is.finite(nt)||nt<1L)stop('compact bootstrap target workspace exceeds 256 MiB',call.=FALSE)
  if(!is.null(target.block)){
    if(length(target.block)!=1L||!is.finite(target.block)||target.block<1L||target.block!=floor(target.block))
      stop('invalid compact bootstrap target block',call.=FALSE)
    nt<-min(nt,target.block)
  }
  mu<-matrix(0,nb,nr);mass.error<-score.error<-tail.max<-0
  TL<-max(40,(70+abs(lbeta(base[1]+1,base[2]+1)))/(base[1]+1))
  TR<-max(40,(70+abs(lbeta(base[1]+1,base[2]+1)))/(base[2]+1))
  for(jj in .np_udens_beta_blocks(nr,as.integer(nt))) {
    ix<-rows[jj];pl<-plan;pl$target<-plan$target[ix]
    pl$coefficients<-plan$coefficients[ix,,drop=FALSE]
    pl$mass<-plan$mass[ix];pl$condition<-plan$condition[ix]
    P<-.np_udens_compact_operator(pl,exp(st$lv))-1
    for(ii in .np_udens_beta_blocks(nb,as.integer(nf))) {
      th<-thetas[ii,,drop=FALSE];a<-th[,1]+1;b<-th[,2]+1;A<-lbeta(a,b)
      lw<-sweep(st$T%*%t(sweep(th,2,base)),1,log(st$prob),'+')
      lw<-sweep(lw,2,st$A-A,'+')
      zm<-apply(lw,2,.np_udens_beta_logsum);prob<-exp(sweep(lw,2,zm))
      exact<-cbind(digamma(a)-digamma(a+b),digamma(b)-digamma(a+b))
      mass<-max(abs(expm1(zm)));err<-max(abs(crossprod(prob,st$T)-exact))
      if(mass>1e-10||err>1e-12+1e-10*max(abs(exact)))
        stop('compact bootstrap common quadrature mass/score',call.=FALSE)
      tail<-exp(-A-a*TL)*((TL+1)^2/a+2*(TL+1)/a^2+2/a^3)+
        exp(-A-b*TR)*((TR+1)^2/b+2*(TR+1)/b^2+2/b^3)
      if(any(!is.finite(tail))||max(tail)>1e-12)
        stop('compact bootstrap per-fit tail certificate',call.=FALSE)
      mu[ii,jj]<-t(P%*%prob)
      mass.error<-max(mass.error,mass);score.error<-max(score.error,err);tail.max<-max(tail.max,tail)
    }
    rm(P)
  }
  list(mu=mu,diagnostics=list(order=order,nodes=nodes,fit.block=nf,target.block=nt,
    bytes=reserve+8*6*nodes*nf+8*2*nodes*nt,cap=cap,
    mass=mass.error,score=score.error,tail=tail.max,base=base))
}

.np_udens_compact_sd_roots <- function(x,h,grid,kernel,prepared,stream,store,progress) {
  n<-length(x);m<-length(grid);B<-stream$request$B
  refits<-.np_udens_beta_refit_dispatch(matrix(x,ncol=1L),prepared$theta,stream,progress)
  all.theta<-rbind(prepared$theta,refits$theta)
  plan<-.np_udens_compact_operator_plan(grid,h,kernel)
  roots<-.np_udens_beta_root_store(store,m,B)
  held<-8*(6*as.double(m)*n+12*m+8*as.double(B))+16*1024^2
  batch<-min(B,max(1L,floor(8*1024^2/(8*as.double(n)))))
  rhat<-prepared$center-prepared$fitted;scale<-sqrt(prepared$empirical_variance)
  # Compute the all-fit measures once; store-row boundaries are not quadrature boundaries.
  main<-.np_udens_compact_common(plan,all.theta,128L,held.bytes=held)
  refined<-.np_udens_compact_common(plan,all.theta,256L,held.bytes=held)
  checks<-list(refinement=.np_udens_beta_check(main$mu,refined$mu,'compact bootstrap common guide means',1e-10),
    original=.np_udens_beta_check(main$mu[1L,],prepared$mu,'compact bootstrap original mean',1e-10),
    main=main$diagnostics,refined=refined$diagnostics)
  rm(refined)
  for(j in seq_along(roots$rows)) {
    ix<-roots$rows[[j]]
    values<-matrix(0,length(ix),B)
    for(ids in .np_udens_beta_blocks(B,as.integer(batch))) {
      counts<-.np_udens_beta_random_read(stream,min(ids),max(ids))
      values[,ids]<-(prepared$observed[ix,,drop=FALSE]%*%(counts/n)-
        t(main$mu[1L+ids,ix,drop=FALSE])-rhat[ix])/scale[ix]
    }
    if(any(!is.finite(values)))stop('nonfinite compact bootstrap roots',call.=FALSE)
    connection<-file(roots$files[[j]],'wb')
    tryCatch(writeBin(as.vector(values),connection,size=8L,endian='little'),finally=close(connection))
    progress('compact bootstrap projection',max(ix),m)
  }
  roots$next_draw<-B+1L
  list(roots=roots,theta=refits$theta,iterations=refits$iterations,checks=checks)
}

.np_udens_compact_sd_hidden <- function(x,h,grid,kernel,alpha,variant,band,B,seed,progress) {
  if(!is.numeric(x)||length(x)<10L||!is.numeric(grid)||!length(grid)||
     any(!is.finite(c(x,grid,h)))||any(x<=0|x>=1)||any(grid<=0|grid>=1)||
     length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1||
     length(variant)!=1L||!variant%in%c('symmetric-sd','asymmetric-sd')||
     length(band)!=1L||!band%in%c('pointwise','bonferroni','simultaneous','all'))
    stop('invalid compact corrected SD request',call.=FALSE)
  original<-grid;grid<-sort(unique(grid));restore<-match(original,grid)
  n<-length(x);m<-length(grid)
  .np_udens_beta_random_request('multinomial',n,B,seed,RNGkind())
  disk<-4*as.double(n)*B+8*as.double(m)*B+262144
  # The maintained root reader admits at most 32 MiB per 64-target block.
  if(!is.finite(disk)||disk>512*1024^2||8*as.double(min(64L,m))*B>32*1024^2||
     8*(8*as.double(B)+6*as.double(m)*n+12*m)+48*1024^2>256*1024^2)
    stop('compact bootstrap request exceeds owned workspace/storage budget',call.=FALSE)
  store<-.np_udens_beta_store();on.exit(.np_udens_beta_store_close(store),add=TRUE)
  prepared<-.np_udens_compact_prepare(x,h,grid,kernel,FALSE,progress)
  region<-.np_udens_beta_regions(prepared$theta,prepared$Veta,matrix(grid,ncol=1L),alpha,progress)
  stream<-.np_udens_beta_random_stream(store,'multinomial',n,B,seed)
  bootstrap<-.np_udens_compact_sd_roots(x,h,grid,kernel,prepared,stream,store,progress)
  rhat<-prepared$center-prepared$fitted
  scales<-list(sqrt(prepared$pilot_variance),sqrt(prepared$empirical_variance))
  if(band%in%c('bonferroni','all'))
    result<-.np_udens_bounded_sd(bootstrap$roots,rhat,scales,region,prepared$center,alpha,band,variant,progress)
  else {
    parts<-lapply(scales,function(s).np_udens_beta_stored_reduce(bootstrap$roots,rhat,s,region,prepared$center,alpha,band,variant))
    result<-parts[[1L]]
    result$lower<-pmin(parts[[1L]]$lower,parts[[2L]]$lower)
    result$upper<-pmax(parts[[1L]]$upper,parts[[2L]]$upper)
  }
  result$center<-prepared$center;result$pilot_se<-scales[[1L]];result$empirical_se<-scales[[2L]]
  for(k in c('lower','upper','center','pilot_se','empirical_se'))result[[k]]<-result[[k]][restore]
  if(band=='all')result$bands<-lapply(result$bands,function(z){z$lower<-z$lower[restore];z$upper<-z$upper[restore];z})
  result$family_size<-m;result$theta<-prepared$theta
  result$checks<-c(prepared$checks,list(bootstrap=bootstrap$checks))
  result$construction<-if(variant=='symmetric-sd')'compact-residual-absolute-containing-hull-envelope'else'compact-signed-residual-inversion-envelope'
  low<-if(band%in%c('bonferroni','all'))result$sd_companion$resolution.limited else
    (if(variant=='asymmetric-sd')floor(B*alpha/4)else floor(B*alpha/2))<20L
  if(low)warning(if(band%in%c('bonferroni','all'))
    'compact plot Bonferroni calibration has low finite-draw tail resolution; returning finite empirical limits with the simultaneous hull'else
    'compact plot calibration has low finite-draw tail resolution; returning the finite band',call.=FALSE)
  result
}
