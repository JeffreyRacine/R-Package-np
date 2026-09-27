ann_compact_literal <- function(x, y, kx, ky, degree, kernel) {
  ker <- if(kernel == "uniform") function(u) 0.5*(abs(u)<1) else
    function(u) 3/(4*sqrt(5))*pmax(0,1-u^2/5)
  support <- if(kernel == "uniform") 1 else sqrt(5)
  mean(vapply(seq_along(y), function(i) {
    donor <- setdiff(seq_along(y),i)
    hx <- vapply(donor,function(j) sort(abs(x[-c(i,j)]-x[j]))[kx],0)
    hy <- vapply(donor,function(j) sort(abs(y[-c(i,j)]-y[j]))[ky],0)
    w <- 3/(4*sqrt(5))*pmax(0,1-((x[donor]-x[i])/hx)^2/5)/hx
    Z <- vapply(0:degree,function(p) (x[donor]-x[i])^p,numeric(length(donor)))
    influence <- drop(w*(Z %*% solve(crossprod(Z,w*Z),c(1,rep(0,degree)))))
    density <- function(q) vapply(q,function(t) sum(influence*ker((t-y[donor])/hy)/hy),0)
    cuts <- sort(unique(c(y[donor]-support*hy,y[donor]+support*hy)))
    # Support cuts make density squared a polynomial of degree at most four.
    # Three-point Gauss is exact here and avoids adaptive-integrator roundoff
    # warnings on nearly coincident uniform-kernel endpoints.
    integral <- sum(vapply(seq_len(length(cuts)-1L),function(j) {
      half <- (cuts[j+1L]-cuts[j])/2
      mid <- cuts[j]+half
      half*sum(c(5/9,8/9,5/9)*density(mid+half*c(-sqrt(3/5),0,sqrt(3/5)))^2)
    },0))
    2*density(y[i])-integral
  },0))
}

test_that("ANN compact overlap blocks agree with direct deleted-density integration", {
  old <- options(np.messages=FALSE,np.extendednn=FALSE,np.macMseries.accelerate=FALSE)
  on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  set.seed(4717)
  x <- runif(17,-1,1); y <- rnorm(17)
  for(kernel in c("epanechnikov","uniform")) for(degree in 0:2) {
    reference <- ann_compact_literal(x,y,8L,4L,degree,kernel)
    args <- list(xdat=data.frame(x=x),ydat=data.frame(y=y),bws=c(4,8),
                 bwtype="adaptive_nn",bwmethod="cv.ls",
                 regtype=if(degree==0L) "lc" else "lp",
                 cxkertype="epanechnikov",cykertype=kernel,bandwidth.compute=FALSE)
    if(degree>0L) {
      args$degree <- degree; args$degree.select <- "manual"; args$bernstein.basis <- TRUE
    }
    b <- do.call(npcdensbw,args)
    for(tree in c(FALSE,TRUE)) {
      options(np.tree=tree)
      call.args <- list(xdat=args$xdat,ydat=args$ydat,bws=b,invalid.penalty="dbmax")
      call.args$force.local <- FALSE
      command <- substitute({
        options(np.messages=FALSE,np.tree=TREE,np.extendednn=FALSE,np.macMseries.accelerate=FALSE)
        do.call(get(".npcdensbw_eval_only",asNamespace("npRmpi")),ARGS)
      },list(TREE=tree,ARGS=call.args))
      observed <- getFromNamespace(".npRmpi_bcast_cmd_expr","npRmpi")(
        command,comm=1L,caller.execute=TRUE)$objective
      expect_lte(abs(observed-reference),1e-10+1e-10*abs(reference))
    }
  }
})
