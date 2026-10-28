# Exact GNN geometry and center algebra on bounded marginal-coordinate blocks.
.np_udist_gnn_blocks <- function(state, grid, size=25L) {
 numeric<-which(state$numeric);categories<-vapply(state$category,`[[`,integer(1),'index')
 groups<-if(length(categories))split(seq_len(nrow(grid)),interaction(grid[categories],drop=TRUE))else list(seq_len(nrow(grid)))
 out<-list()
 for(ii in groups) {
  axes<-lapply(numeric,function(j)unique(grid[[j]][ii]))
  chunks<-lapply(axes,function(a)split(seq_along(a),ceiling(seq_along(a)/size)))
  for(a in chunks[[1]])for(b in if(length(axes)==2)chunks[[2]]else list(1L)) {
   u<-list(axes[[1]][a]);if(length(axes)==2)u[[2]]<-axes[[2]][b]
   keep<-grid[[numeric[1]]][ii]%in%u[[1]]
   if(length(axes)==2)keep<-keep&grid[[numeric[2]]][ii]%in%u[[2]]
   jj<-ii[keep];if(!length(jj))next
   out[[length(out)+1]]<-list(rows=jj,axes=u,map=lapply(seq_along(numeric),function(k)match(grid[[numeric[k]]][jj],u[[k]])),profile=grid[jj[1],,drop=FALSE])
  }
 }
 out
}
.np_udist_gnn_block <- function(state, block, corrected = TRUE) {
 numeric<-which(state$numeric)
 cache<-lapply(seq_along(numeric),function(k) {
  delta<-outer(block$axes[[k]],state$x[[numeric[k]]],'-')
  order<-t(apply(abs(delta),1L,order))
  if(nrow(delta)==1)order<-matrix(order,1L)
  list(delta=delta,order=order,k=state$rank$lookup[k])
 })
 cat<-rep(1,state$n)
 for(c in state$category)cat<-cat*.np_udist_bc_ordered(
   as.numeric(as.character(block$profile[[c$index]]))-c$training,c$lambda,corrected)
 list(cache=cache,category=cat,map=block$map,n=state$n,corrected=corrected)
}
.np_udist_gnn_center <- function(prepared,counts) {
 stopifnot(length(counts)==prepared$n,sum(counts)==prepared$n,all(counts>=0))
 h<-list();margins<-lapply(prepared$cache,function(c) {
  radii<-vapply(seq_len(nrow(c$delta)),function(i){o<-c$order[i,];abs(c$delta[i,o[which(cumsum(counts[o])>=c$k)[1]]])},numeric(1))
  if(any(!is.finite(radii))||any(radii<=0))stop('exact GNN center unavailable: nonpositive physical radius')
  h[[length(h)+1]]<<-radii
  z<-c$delta/radii;D<-z*dnorm(z);D[is.infinite(z)]<-0;list(H=pnorm(z),D=D)
 })
 w<-counts*prepared$category
 A<-margins[[1]]
 if(length(margins)==1)center<-drop((if(prepared$corrected)A$H+.5*A$D else A$H)%*%w)/prepared$n
 else {
  B<-margins[[2]];AH<-sweep(A$H,2L,w,'*');AD<-sweep(A$D,2L,w,'*')
  center<-(AH%*%t(B$H)+if(prepared$corrected).5*(AD%*%t(B$H)+AH%*%t(B$D))else 0)/prepared$n
 }
 list(center=if(length(margins)==1)center[prepared$map[[1]]]else center[cbind(prepared$map[[1]],prepared$map[[2]])],
      radii=do.call(cbind,lapply(seq_along(h),function(k)h[[k]][prepared$map[[k]]])))
}

# Preserve the observed center's original recurrence and shifted summation.
# Radius computation here needs no spacing pilot or bootstrap studentizer.
.np_udist_gnn_main <- function(state, grid, corrected = TRUE) {
  n <- state$n; numeric <- which(state$numeric)
  P <- matrix(1, nrow(grid), n); Q <- P * 0
  for (a in seq_along(numeric)) {
    j <- numeric[a]; k <- state$rank$lookup[a]; distinct <- unique(grid[[j]])
    radius <- vapply(distinct, function(t)
      sort.int(abs(state$x[[j]] - t), partial = k)[k], numeric(1))
    h <- radius[match(grid[[j]], distinct)]
    if (any(!is.finite(h)) || any(h <= 0))
      stop("exact GNN center unavailable: nonpositive physical radius", call. = FALSE)
    u <- outer(grid[[j]], state$x[[j]], "-") / h
    H <- pnorm(u); D <- u * dnorm(u); D[is.infinite(u)] <- 0
    Q <- Q * H + P * D; P <- P * H
  }
  original <- P; G <- if (corrected) P + .5 * Q else P
  for (cat in state$category) {
    distance <- outer(as.numeric(as.character(grid[[cat$index]])), cat$training, "-")
    original <- original * .np_udist_bc_ordered(distance, cat$lambda, FALSE)
    G <- G * .np_udist_bc_ordered(distance, cat$lambda, corrected)
  }
  list(center = G[, 1L] + rowMeans(G - G[, 1L]),
    original.center = original[, 1L] + rowMeans(original - original[, 1L]))
}

.np_udist_gnn_budget <- function(state, m, B, processes = 1L) {
  n <- as.double(state$n); p <- ncol(state$x)
  base <- 8 * (16*n*p + 40*n*25 + 8*n*32 + 12*B + 12*m)
  capacity <- floor((state$max.bytes / processes - base) / (6 * 8 * B))
  side <- min(25L, floor(sqrt(max(0, capacity))))
  if (!is.finite(side) || side < 1L)
    stop("exact GNN resident memory forecast exceeds budget", call. = FALSE)
  per.process <- base + 6 * 8 * side^2 * B
  list(bytes = processes * per.process, per.process = per.process,
    processes = processes, marginal.tile = as.integer(side),
    target.tile = as.integer(side^2), draw.chunk = 32L, covariance.cap = 32L)
}
