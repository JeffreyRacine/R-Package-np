# Plot-only literal counts/refits, complete common base and bounded root storage.
.np_udens_compact_q2_sd_roots <- function(z,stream,store,progress=function(...)NULL) {
 n<-nrow(z$x);m<-nrow(z$grid);B<-stream$request$B
 refits<-.np_udens_beta_refit_dispatch(z$x,z$fit$theta,stream,progress)
 all.theta<-rbind(z$fit$theta,refits$theta)
 base<-c((apply(all.theta[,1:4,drop=FALSE],2,min)+1)*.5-1,0)
 axes<-lapply(c(16L,32L),function(o).np_udens_compact_q2_axes(base,z$grid,z$h,z$kernel,o))
 nodes<-vapply(axes[[2L]],`[[`,integer(1),'nodes')
 unique.axes<-vapply(1:2,function(j)length(unique(z$grid[,j])),integer(1))
 # Two complete projections coexist; fit/root batches do not set their measures.
 forecast<-as.double(object.size(z))+as.double(object.size(all.theta))+8*(2*(24*sum(unique.axes*nodes)+257*sum(nodes)+24*257*sum(unique.axes))+12*32*m+40*B+8*32*n)+24*1024^2
 if(any(nodes>8192L)||!is.finite(forecast)||forecast>192*1024^2)stop('compact SD complete workspace budget')
 projections<-lapply(1:2,function(i){rule<-list(base=base,axes=axes[[i]],A=axes[[i]][[1]]$A+axes[[i]][[2]]$A);.np_udens_beta_series_projection(rule,z$h,n,z$grid,c(16L,32L)[i],.np_udens_compact_q2_parts(z$kernel))})
 original<-lapply(projections,function(p).np_udens_beta_series_means_native(p,matrix(z$fit$theta,1L)))
 for(k in c('mu','A','fitted')) .np_udens_beta_check(original[[1L]][[k]],original[[2L]][[k]],paste('SD original',k),1e-10)
 .np_udens_beta_check(drop(original[[1L]]$mu),z$mu,'SD original preparation mean',1e-10)
 roots<-.np_udens_beta_root_store(store,m,B);normal<-numeric(B);maximum<-c(mu=0,A=0,fitted=0)
 rhat<-z$center-z$fitted;scale<-sqrt(z$empirical_variance)
 for(ix in .np_udens_beta_blocks(B,32L)) {
  a<-.np_udens_beta_series_means_native(projections[[1L]],refits$theta[ix,,drop=FALSE])
  b<-.np_udens_beta_series_means_native(projections[[2L]],refits$theta[ix,,drop=FALSE])
  for(k in c('mu','A','fitted')){.np_udens_beta_check(a[[k]],b[[k]],paste('SD refit',k),1e-10);maximum[k]<-max(maximum[k],max(abs(a[[k]]-b[[k]])))}
  counts<-.np_udens_beta_random_read(stream,min(ix),max(ix))
  values<-(z$observed%*%(counts/n)-t(a$mu)-rhat)/scale
  roots<-.np_udens_beta_root_write(roots,values,min(ix));normal[ix]<-a$A
  progress('literal refitted SD roots',max(ix),B)
 }
 list(roots=roots,theta=refits$theta,iterations=refits$iterations,normalization=c(original[[1L]]$A,normal),checks=list(refinement=maximum,base=base,nodes=nodes,forecast=forecast))
}
