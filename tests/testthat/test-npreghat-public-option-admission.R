test_that("public npreghat matrix and apply retain option admission and pool reuse", {
  skip_on_cran()
  if(!spawn_mpi_slaves(1L))skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  withr::local_options(np.messages=FALSE,np.tree=FALSE,npRmpi.autodispatch=TRUE,
                      npRmpi.autodispatch.option.sync="onchange")
  x <- data.frame(x=seq(.03,.97,length.out=25)); y <- sin(3*x$x)
  b <- npregbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE)
  K <- exp(-outer(x$x,x$x,"-")^2/(2*.3^2)); oracle <- K/rowSums(K)
  ns <- asNamespace("npRmpi")
  witness <- new.env(parent=emptyenv()); witness$n <- 0L
  trace(".npRmpi_autodispatch_materialize_call",tracer=function()witness$n<-witness$n+1L,
        where=ns,print=FALSE)
  on.exit(untrace(".npRmpi_autodispatch_materialize_call",where=ns),add=TRUE)
  for (output in c("matrix","apply")) {
    options(npRmpi.autodispatch.option.sync="never")
    before <- witness$n
    expect_error(npreghat(b,txdat=x,y=y,output=output),"verified option agreement",fixed=TRUE)
    expect_identical(witness$n,before)
    for (mode in c("always","onchange")) {
      options(npRmpi.autodispatch.option.sync=mode)
      actual <- npreghat(b,txdat=x,y=y,output=output)
      expected <- if(output=="matrix")oracle else drop(oracle%*%y)
      expect_equal(as.numeric(actual),as.numeric(expected),tolerance=2e-12)
    }
  }
  expect_true(.mpi_pool_active())
})
