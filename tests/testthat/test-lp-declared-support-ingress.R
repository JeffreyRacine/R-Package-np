# C179 support qualification; deleted-NN geometry is separately tracked.
lp_declared_support_contract <- function(package) {
  ns <- asNamespace(package)
  bwfun <- get("npregbw", ns); fitfun <- get("npreg", ns)
  hatfun <- get("npreghat", ns)
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  x <- data.frame(x=c(.04,.12,.19,.27,.34,.41,.48,.53,.62,.71,.78,.91),
    u=factor(rep(c(1,2,4),4),levels=1:16),
    o=ordered(rep(c(.5,1.5,3.5),4),levels=seq(.5,15.5)))
  y <- sin(4*x$x)+as.integer(x$u)/8+cos(seq_len(nrow(x)))/9
  responses <- cbind(y, y^2+.2)
  ex <- x[c(11,2,8,4),,drop=FALSE]
  for (type in c("fixed","generalized_nn","adaptive_nn"))
    for (degree in c(0L,1L,2L)) for (tree in c(FALSE,TRUE)) {
      options(np.tree=tree)
      b <- bwfun(xdat=x,ydat=y,bws=c(if(type=="fixed").3 else 7,.25,.3),
        bwtype=type,regtype="lp",degree=degree,okertype="racineliyan",
        bandwidth.compute=FALSE)
      H <- hatfun(bws=b,txdat=x,exdat=ex)
      applied <- hatfun(bws=b,txdat=x,exdat=ex,y=responses,output="apply")
      direct <- vapply(seq_len(ncol(responses)),function(j)
        fitted(fitfun(bws=b,txdat=x,tydat=responses[,j],exdat=ex)),
        numeric(nrow(ex)))
      expect_equal(unname(H%*%responses),unname(direct),tolerance=2e-10)
      expect_equal(unname(applied),unname(direct),tolerance=2e-10)
      # Continuous derivative with identical category support.
      G <- hatfun(bws=b,txdat=x,exdat=ex,s=1L)
      derivative <- fitfun(bws=b,txdat=x,tydat=y,exdat=ex,
                           gradients=TRUE,gradient.order=1L)$grad[,1]
      expect_equal(as.double(G%*%y),as.double(derivative),tolerance=2e-10)
    }
  options(np.tree=FALSE)
  b <- bwfun(xdat=x,ydat=y,bws=c(.3,.25,.3),regtype="lp",degree=2L,
    okertype="racineliyan",bandwidth.compute=FALSE)
  tile <- get(".np_npsig_streamed_iid_tile",ns)
  statistic <- get(".np_npsig_statistic",ns)
  for (coordinate in 1:3) {
    expected <- vapply(seq_len(ncol(responses)),function(j)
      statistic(fitfun(bws=b,txdat=x,tydat=responses[,j],gradients=TRUE,se=TRUE),
                coordinate,FALSE),numeric(1))
    result <- tile(b,x,coordinate,response.matrix=responses,
      null.mean=y,residual.pool=y,pivotal=FALSE)
    expect_equal(result,expected,tolerance=2e-10)
  }
}
test_that("LP operators and tiles preserve declared training categories", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  lp_declared_support_contract("npRmpi")
})
