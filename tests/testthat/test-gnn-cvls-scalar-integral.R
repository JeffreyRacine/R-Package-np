gnn_scalar_literal <- function(y, k, kernel) {
  ker <- if(kernel=="gaussian") dnorm else
    function(u) 3/(4*sqrt(5))*pmax(0,1-u^2/5)
  radius <- function(q,v) sort(abs(q-v))[k]
  density <- function(q) vapply(q,function(t) {
    h <- radius(t,y); mean(ker((t-y)/h)/h)
  },0)
  cuts <- c(-Inf,sort(unique(c(y,as.vector(outer(y,y,"+")/2)))),Inf)
  if(kernel=="epanechnikov")
    cuts <- sort(unique(c(cuts,
      as.vector(outer(y,y,function(a,b) a+(a-b)/(sqrt(5)-1))),
      as.vector(outer(y,y,function(a,b) a+(a-b)/(-sqrt(5)-1))))))
  I1 <- sum(vapply(seq_len(length(cuts)-1L),function(j)
    integrate(function(q) density(q)^2,cuts[j],cuts[j+1L],
              abs.tol=1e-12/length(cuts),rel.tol=1e-10,subdivisions=1000)$value,0))
  I2 <- mean(vapply(seq_along(y),function(i) {
    h <- radius(y[i],y[-i]); mean(ker((y[i]-y[-i])/h)/h)
  },0))
  2*I2-I1
}

test_that("scalar GNN CVLS integrates the literal query-dependent density", {
  old <- options(np.messages=FALSE,np.extendednn=FALSE,np.macMseries.accelerate=FALSE)
  on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  y <- c(-2,-1.4,-.8,-.5,-.1,.2,.55,.8,1.2,1.7,2.1,3)
  for(kernel in c("gaussian","epanechnikov")) {
    reference <- gnn_scalar_literal(y,3L,kernel)
    for(tree in c(FALSE,TRUE)) for(scale in c(1,4)) {
      options(np.tree=tree)
      dat <- data.frame(y=1024+scale*y)
      b <- npudensbw(dat=dat,bws=3,bwtype="generalized_nn",bwmethod="cv.ls",
                     ckertype=kernel,bandwidth.compute=FALSE)
      call.args <- list(dat=dat,bws=b,eval.only=TRUE,invalid.penalty="dbmax",nmulti=1L)
      command <- substitute({
        options(np.messages=FALSE,np.tree=TREE,np.extendednn=FALSE,np.macMseries.accelerate=FALSE)
        do.call(get("npudensbw.bandwidth",asNamespace("npRmpi")),ARGS)
      },list(TREE=tree,ARGS=call.args))
      observed <- getFromNamespace(".npRmpi_bcast_cmd_expr","npRmpi")(
        command,comm=1L,caller.execute=TRUE)$fval
      expect_lte(abs(observed-reference/scale),1e-10+1e-10*abs(reference/scale))
    }
  }
})
