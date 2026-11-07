# Signed density bootstrap-t. This entry is private and leaves v1 untouched.
.np_udens_b2_boot_request <- function(law, blocklen, n) {
  law <- match.arg(law,c("inid","fixed","geom"))
  if (law %in% c("fixed", "geom")) {
    blocklen <- .np_udens_b2_block_length(blocklen, n)
    if (law == "fixed" && blocklen == n)
      stop("fixed block length equal to sample size has no resampling variation", call. = FALSE)
    return(list(law = law, blocklen = blocklen))
  }
  if (!is.null(blocklen) && (!is.numeric(blocklen) || length(blocklen)!=1L ||
      !is.finite(blocklen) || blocklen!=1))
    stop("inid density bootstrap requires block length one",call.=FALSE)
  list(law=law,blocklen=1L)
}

.np_udens_b2_indices <- function(n, first, last, B, law, blocklen) {
  if(law!="inid"&&blocklen>1L){
    draw<-.np_block_replicate_drawer(n,n,blocklen,law,TRUE)
    return(lapply(seq_len(last-first+1L),function(j)draw()))
  }
  counts <- .np_bias_process_innovations(n,first,last,B,law,blocklen)+1
  lapply(seq_len(ncol(counts)),function(j)rep.int(seq_len(n),counts[,j]))
}

.np_udens_b2_variance <- function(F, law, blocklen) {
  if(law=="inid"||blocklen==1L)return(rowSums(F^2))
  .np_bias_score_variance(F,blocklen,law)
}

.np_udens_b2_boot_chunk <- function(state, grid, base, first, last, B, law, blocklen) {
  if(law!="inid"&&blocklen>1L)
    return(.np_udens_b2_block_chunk(state,grid,base,first,last,B,law,blocklen))
  n <- state$n
  counts <- .np_bias_process_innovations(n,first,last,B,law,blocklen)+1
  z <- .np_udens_bc_inid_pivot(base$factor,base$corrected,counts,
                              accurate=state$r+state$s>0L)
  signed <- z$absolute*sign(base$factor%*%counts)
  if (anyNA(signed)) stop("undefined signed density bootstrap pivot",call.=FALSE)
  list(signed=signed,degeneracy=z$degeneracy,cancellation=z$cancellation)
}

# New private MPI transport uses the incumbent chunk/task/RNG owners unchanged.
.np_udens_b2_boot_dispatch <- function(state, grid, base, B, law, blocklen) {
  .npRmpi_bootstrap_fanout_enabled(n=state$n,B=B,chunk.size=32L,
                                  what="signed density bootstrap")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B,32L,with.seeds=FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks,function(task)
    invisible(.np_udens_b2_indices(state$n,task$start,task$start+task$bsz-1L,
                                   B,law,blocklen)),
    normal.draws=identical(law,"geom") && blocklen>1L)
  worker <- function(task,state,grid,base,B,law,blocklen) {
    .npRmpi_bootstrap_use_task_rng(task)
    z <- .np_udens_b2_boot_chunk(state,grid,base,task$start,
                                 task$start+task$bsz-1L,B,law,blocklen)
    out <- t(z$signed)
    attr(out,"density.b2") <- z[c("degeneracy","cancellation")]
    out
  }
  reducer <- function(out,parts,tasks) {
    d <- lapply(parts,attr,which="density.b2",exact=TRUE)
    if(any(vapply(d,is.null,logical(1))))stop("density worker omitted diagnostics",call.=FALSE)
    attr(out,"density.b2") <- list(degeneracy=Reduce(`+`,lapply(d,`[[`,"degeneracy")),
      cancellation=Reduce(`+`,lapply(d,`[[`,"cancellation")))
    out
  }
  deps <- list(.np_udens_b2_boot_chunk=.np_udens_b2_boot_chunk,
    .np_udens_b2_indices=.np_udens_b2_indices,.np_udens_b2_variance=.np_udens_b2_variance,
    .np_udens_b2_block_chunk=.np_udens_b2_block_chunk,.np_udens_b2_block_zero=.np_udens_b2_block_zero)
  out <- .npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=nrow(grid),
    what="signed density bootstrap",profile.where="mpi.applyLB:density-b2",
    progress.label="Calibrating signed density bootstrap",progress.unit="draws",
    required.bindings=deps,metadata.reducer=reducer,state=state,grid=grid,base=base,
    B=B,law=law,blocklen=blocklen)
  d <- attr(out,"density.b2",exact=TRUE);attr(out,"density.b2") <- NULL
  list(signed=t(out),degeneracy=d$degeneracy,cancellation=d$cancellation)
}

.np_udens_b2_side_update <- function(side, draws, ii, alpha, m) {
  side$maximum <- pmax(side$maximum,apply(draws,2L,max))
  for (j in seq_along(ii)) {
    q <- unname(quantile(draws[j,],c(1-alpha/2,1-alpha/(2*m)),type=1L))
    side$point[ii[j],] <- q[seq_along(alpha)]
    side$bonf[ii[j],] <- q[length(alpha)+seq_along(alpha)]
  }
  side
}

.np_udens_b2_bootstrap <- function(state, grid, law="inid", blocklen=NULL,
                                   alpha=c(.05,.1), B=9999L, tile=25L,
                                   corrected=TRUE) {
  existed <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  previous <- if(existed)get(".Random.seed",.GlobalEnv,inherits=FALSE)else NULL
  kind <- RNGkind();success <- FALSE
  on.exit({if(!success){do.call(RNGkind,as.list(kind));if(existed)
    assign(".Random.seed",previous,.GlobalEnv) else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
      rm(list=".Random.seed",envir=.GlobalEnv)}},add=TRUE)
  request <- .np_udens_b2_boot_request(law,blocklen,state$n)
  law <- request$law;blocklen <- request$blocklen
  if(!is.numeric(alpha)||!length(alpha)||any(!is.finite(alpha))||
     any(alpha<=0|alpha>=.5)||anyDuplicated(alpha))stop("invalid density alpha",call.=FALSE)
  B <- .np_plot_gaussian_integer(B,"density draws",.Machine$integer.max)
  if(B<2L)stop("density calibration requires at least two draws",call.=FALSE)
  .np_udens_bc_grid(state,grid);m <- nrow(grid)
  budget <- .np_udens_b2_budget(state,m,B,tile,FALSE)
  if(!existed)runif(1L)
  start <- get(".Random.seed",.GlobalEnv,inherits=FALSE);end <- NULL
  .np_bias_normal_replay_check(law, blocklen, m > tile, seed = start)
  fresh <- function()list(maximum=rep(-Inf,B),point=matrix(0,m,length(alpha)),bonf=matrix(0,m,length(alpha)))
  plus <- minus <- fresh();absolute <- .np_bias_limits_prepare(rep("density",m),alpha,FALSE,B)
  center <- original <- se <- numeric(m);degenerate <- cancellation <- integer(m)
  for(block in .np_bias_target_blocks(m,1L,anchors=tile)) {
    ii <- block$grid;g <- grid[ii,,drop=FALSE]
    value <- .np_udens_bc_contribution(state,g,corrected)
    scale <- sqrt(.np_udens_b2_variance(value$factor,law,blocklen))
    if(any(!is.finite(scale))||any(scale<=0))stop("density bootstrap requires positive observed scales",call.=FALSE)
    stream <- .np_bias_stream(start,.np_udens_b2_boot_dispatch(state,g,value,B,law,blocklen))
    if(!is.null(end)&&!identical(end,stream$end))stop("density bootstrap tile RNG mismatch",call.=FALSE)
    end <- stream$end;z <- stream$value
    plus <- .np_udens_b2_side_update(plus,z$signed,ii,alpha,m)
    minus <- .np_udens_b2_side_update(minus,-z$signed,ii,alpha,m)
    absolute <- .np_udens_bc_limits_update(absolute,abs(z$signed),ii,z$degeneracy)
    center[ii] <- value$center;original[ii] <- value$original.center;se[ii] <- scale
    degenerate[ii] <- z$degeneracy;cancellation[ii] <- z$cancellation
  }
  # Both constructions use precisely the same draws, centers and scales.
  make.intervals <- function(lower,upper) {
    out <- lapply(names(lower),function(b){z<-array(0,c(m,2L,length(alpha)),
      dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
      z[,1L,]<-center-se*lower[[b]];z[,2L,]<-center+se*upper[[b]];z})
    names(out)<-names(lower)
    for(b in c("simultaneous","bonferroni")){
      prev<-if(b=="simultaneous")"pointwise"else"simultaneous"
      out[[b]][,1L,]<-pmin(out[[b]][,1L,],out[[prev]][,1L,])
      out[[b]][,2L,]<-pmax(out[[b]][,2L,],out[[prev]][,2L,])
    }
    if(any(vapply(out,function(z)anyNA(z)||any(z[,1L,]>z[,2L,]),logical(1))))
      stop("invalid density bootstrap interval ordering",call.=FALSE)
    out
  }
  side.finish <- function(s)list(pointwise=s$point,
    simultaneous=matrix(rep(unname(quantile(s$maximum,1-alpha/2,type=1L)),each=m),m),bonferroni=s$bonf)
  signed <- make.intervals(side.finish(plus),side.finish(minus))
  fin <- .np_bias_limits_finish(absolute)
  critical <- lapply(fin$limits,function(z)matrix(z[,2L,],m,length(alpha)))
  abs.intervals <- make.intervals(critical,critical)
  assign(".Random.seed",end,.GlobalEnv);success<-TRUE
  list(center=center,original.center=original,se=se,intervals=signed,
    absolute.intervals=abs.intervals,degeneracy=degenerate,cancellation=cancellation,
    settings=list(construction="density-bc-v2-signed",law=law,blocklen=blocklen,
      alpha=alpha,B=B,tile=tile,budget=budget,corrected=corrected),
    rng=list(kind=kind,entry=previous,start=start,end=end))
}
