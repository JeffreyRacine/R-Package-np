test_that("adaptive derivative hats retain absolute donor weights", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(-1.2,-.9,-.5,-.2,.15,.5,.9,1.4,2))
  e <- data.frame(x=c(-.45,.33,1.1))
  y <- sin(x$x)+x$x^2/7
  b <- npregbw(xdat=x,ydat=y,bws=4,bwtype="adaptive_nn",
                bandwidth.compute=FALSE,regtype="lc")
  h <- vapply(x$x,function(t)sort(abs(x$x-t))[5L],0)
  delta <- outer(x$x,e$x,"-")
  w <- dnorm(delta/h)/h
  dw <- delta*w/h^2
  expected <- t(sweep(dw,2L,colSums(w),"/")-
    sweep(w,2L,colSums(dw)/colSums(w)^2,"*"))
  got <- .npreghat_exact_lc_derivative_matrix_from_npksum_chunked(b,x,e,s=1L)
  expect_equal(unname(got),unname(expected),tolerance=3e-12)
  fit <- npreg(b,txdat=x,tydat=y,exdat=e,gradients=TRUE,se=FALSE)
  expect_equal(as.double(got %*% y),as.double(gradients(fit)),tolerance=3e-12)
})

test_that("frozen adaptive LP bootstrap retains the donor normalization", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(-1.2,-.9,-.5,-.2,.15,.5,.9,1.4,2))
  e <- data.frame(x=c(-.45,.33,1.1))
  y <- sin(x$x)+x$x^2/7
  b <- npregbw(xdat=x,ydat=y,bws=4,bwtype="adaptive_nn",
                bandwidth.compute=FALSE,regtype="ll")
  h <- vapply(x$x,function(t)sort(abs(x$x-t))[5L],0)
  weights <- dnorm(outer(x$x,e$x,"-")/h)/h
  counts <- cbind(rep(1,9),rep(c(0,2,1),3),rep(c(2,1,0),3))
  design <- cbind(1,x$x)
  for(gradient in c(FALSE,TRUE)) {
    rhs <- if(gradient)cbind(0,rep(1,nrow(e)))else cbind(1,e$x)
    expected <- vapply(seq_len(nrow(e)),function(j) {
      vapply(seq_len(ncol(counts)),function(i) {
        ww <- weights[,j]*counts[,i]
        beta <- solve(crossprod(design,ww*design),crossprod(design,ww*y))
        drop(rhs[j,] %*% beta)
      },0)
    },numeric(ncol(counts)))
    got <- .np_inid_boot_from_regression_localpoly_frozen(x,e,b,y,3L,
      counts=counts,gradients=gradient)
    expect_equal(got$t,expected,tolerance=3e-12)
    expect_equal(as.double(got$t0),as.double(expected[1,]),tolerance=3e-12)
  }
})
