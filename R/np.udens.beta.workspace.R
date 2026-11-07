# Additional inference arrays, not whole-process RSS. Counts are doubles for
# arithmetic; both retained numerical rules and serialization copies are charged.
.np_udens_beta_workspace <- function(n,grid,theta,full,B=9999L,workers=0L) {
  m<-nrow(grid);q<-ncol(grid);tile<-min(m,64L);cap<-256*1024^2
  if(!is.numeric(n)||length(n)!=1L||!is.finite(n)||n<10||!is.numeric(B)||length(B)!=1L||!is.finite(B)||B<1||
     !is.numeric(workers)||length(workers)!=1L||!is.finite(workers)||workers<0)
    stop('invalid beta plot workspace geometry',call.=FALSE)
  if(q==1L) {
    nodes<-.np_udens_beta_marginal_rule(theta,256L)$nodes
    rows<-if(full)m else tile
    moments<-8*(10*rows*nodes+5*rows^2+6*rows*n+32*nodes+12*n+32*m)
    # Common-guide rules use bounded draw tiles, but retain B-by-target means.
    bootstrap<-8*(10*tile*8192+10*8192*32+8*B*tile+6*n*tile+20*B+12*n)
    factor<-8*(4*m^2+6*m*64+4*B+20*m)
    stages<-c(moments=moments,calibration=if(full)factor else bootstrap)
  } else {
    nodes<-vapply(c(32L,64L),function(order)vapply(1:2,function(j)
      .np_udens_beta_marginal_rule(theta[(2*j-1):(2*j)],order)$nodes,integer(1)),integer(2))
    unique.axes<-vapply(1:2,function(j)length(unique(grid[,j])),integer(1))
    if(any(nodes>1600L)||any(apply(nodes,2,prod)>2560000))stop('beta plot joint refinement node cap',call.=FALSE)
    # Main/refined cached axis factors, joint weights, guide/data/output arrays.
    retained<-8*(24*sum(unique.axes*rowSums(nodes))+sum(apply(nodes,2,prod))+16*n+40*m+16*sum(nodes))
    preparation<-retained+8*(5*n*max(nodes)+2*max(apply(nodes,2,prod)))
    moments<-retained+8*(3*24*n*sum(pmin(unique.axes,tile))+12*tile*n+12*tile*256+20*tile^2)
    covariance<-retained+8*(4*m^2+4*m*256+40*m)
    # Bound SD numerical-series factors independently of target Cartesian layout.
    series.nodes<-vapply(c(32L,64L),function(order)vapply(.np_udens_beta_series_rule(theta,order)$axes,`[[`,integer(1),'nodes'),integer(2))
    projection<-8*(24*sum(unique.axes*rowSums(series.nodes))+3*257*sum(series.nodes[,2])+3*24*257*sum(unique.axes)+16*B+12*m*32+12*n*32+40*m)
    stages<-c(preparation=preparation,moments=moments,calibration=if(full)covariance else projection)
  }
  # The fanout leaves receive small target tiles or bounded count waves. Charge
  # data/guide copies and native region/refit scratch over the complete pool.
  transport<-if(workers>0)32*1024^2+8*(workers+1)*(4*n*q+(if(q==2L)4*257*sum(series.nodes[,1]) else 4*n)+8192) else 0
  # MPI work occurs after preparation returns; its arrays do not overlap the
  # dense covariance preparations. Retained output and the B-by-guide fits do.
  retained.outputs<-8*(40*m+20*B+if(full)m^2 else 0)
  # Region mass-only certificates have at most 600,000 joint nodes. Seven
  # node-length double temporaries plus axes, optimizer histories/native stack
  # fit inside 40 MiB per rank. Collect after each sign so these certificates
  # do not accumulate across targets. Charge all ranks plus retained covariance.
  # Region fanout sends only theta/Veta/target tiles, not bootstrap series axes.
  region<-if(q==2L)40*1024^2*(workers+1L) else 0
  stages<-c(stages,transport=transport+retained.outputs,
            region=(if(workers>0)32*1024^2 else 0)+retained.outputs+region)
  bytes<-max(stages)+4*1024^2
  if(!is.finite(bytes)||bytes>cap)
    stop(sprintf('beta plot inference workspace forecast exceeds 256 MiB (%s stage)',names(which.max(stages))),call.=FALSE)
  list(bytes=bytes,stages=stages,cap=cap)
}
