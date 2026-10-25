test_that("pooled auxiliary capacity is checked before allocating the basis", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  local_mocked_bindings(.np_bias_profile_basis=function(...)
    stop("reached bounded basis allocation",call.=FALSE),.package="npRmpi")
  for(P in c(255L,256L,257L,383L,384L,385L)) {
    G <- P-5L;set.seed(50000L+P)
    x <- data.frame(x=runif(2L*G),z=factor(rep(seq_len(G),each=2L)))
    y <- rnorm(nrow(x))
    b <- npregbw(xdat=x,ydat=y,bws=c(.5,.2),regtype="lc",bandwidth.compute=FALSE)
    expected <- if(P<=384L)"reached bounded basis allocation" else "384-column workspace"
    expect_error(npRmpi:::.np_bias_regression_prepare(b,x,y,np_bias_control(degree=5L)),expected)
  }
})
