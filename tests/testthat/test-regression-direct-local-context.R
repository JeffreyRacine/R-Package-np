b1_context_snapshot <- function() {
  list(options = lapply(c("npRmpi.local.regression.mode",
    "npRmpi.autodispatch.disable","npRmpi.autodispatch.context"),getOption),
    size = mpi.comm.size(1L), rank = mpi.comm.rank(1L))
}
b1_ann_oracle <- function(x,y,ex,k,lambda) {
  n <- nrow(x)
  h <- vapply(seq_len(n),function(i)
    sort(abs(x$x[-i]-x$x[i]))[k],numeric(1))
  t(vapply(seq_len(nrow(ex)),function(i) {
    dx <- x$x-ex$x[i]
    w <- dnorm(dx/h)/h * ifelse(x$u==ex$u[i],1-lambda,lambda/2)
    wp <- dx/h^2*w
    m <- sum(w*y)/sum(w)
    c(mean=m,gradient=(sum(wp*y)-m*sum(wp))/sum(w))
  },c(mean=0,gradient=0)))
}

test_that("direct local ANN corrections share and restore R/native context", {
  if(!spawn_mpi_slaves()) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  direct <- getFromNamespace(".np_regression_direct","npRmpi")
  local <- getFromNamespace(".npRmpi_with_local_regression","npRmpi")
  set.seed(8331)
  x <- data.frame(x=runif(45),u=factor(rep(letters[1:3],15)))
  y <- sin(3*x$x)+.3*(x$u=="b")+rnorm(45,sd=.1)
  ex <- x[c(7,15,22,30),,drop=FALSE]
  bw <- npregbw(xdat=x,ydat=y,bws=c(12,.3),bandwidth.compute=FALSE,
    regtype="lc",bwtype="adaptive_nn")
  before <- b1_context_snapshot(); rng <- .Random.seed
  expected <- b1_ann_oracle(x,y,ex,12,.3)
  for(flag in c(TRUE,FALSE)) {
    actual <- direct(bw,x,y,ex,gradients=TRUE,local.mode=flag)
    expect_equal(actual$mean,expected[,"mean"],tolerance=5e-10)
    expect_equal(actual$grad[,1],expected[,"gradient"],tolerance=5e-10)
    expect_identical(.Random.seed,rng)
    expect_identical(b1_context_snapshot(),before)
  }
  local({
    outer <- b1_context_snapshot()
    expect_identical(outer$size,1L)
    expect_true(isTRUE(outer$options[[1L]]))
    direct(bw,x,y,ex,gradients=TRUE,local.mode=TRUE)
    expect_identical(b1_context_snapshot(),outer)
    bad <- x; bad$x[1:20] <- bad$x[1]
    expect_error(direct(bw,bad,y,ex,gradients=TRUE,local.mode=TRUE),"radius")
    expect_identical(b1_context_snapshot(),outer)
  })
  expect_identical(b1_context_snapshot(),before)
  bad <- x; bad$x[1:20] <- bad$x[1]
  expect_error(direct(bw,bad,y,ex,gradients=TRUE,local.mode=TRUE),"radius")
  expect_identical(b1_context_snapshot(),before)
  expect_equal(direct(bw,x,y,ex,gradients=TRUE)$grad[,1],
    expected[,"gradient"],tolerance=5e-10)
  expect_identical(b1_context_snapshot(),before)
  # Public bootstrap route, formerly errors despite a live pool.
  # Other plot scopes historically normalize absent flags to FALSE; this gate
  # checks their configured-state contract, separately from the exact absent-
  # option restoration of the direct/canonical scope tested above and below.
  old.plot <- options(npRmpi.local.regression.mode=FALSE)
  on.exit(options(old.plot),add=TRUE)
  before <- b1_context_snapshot()
  set.seed(834)
  out <- plot(bw,xdat=x,ydat=y,output="data",gradients=TRUE,
    errors="bootstrap",bootstrap="wild",B=5,neval=5)
  expect_length(out,2L)
  expect_true(all(is.finite(gradients(out[[1L]]))))
  expect_identical(b1_context_snapshot(),before)
})

test_that("canonical local scopes preserve absent configured and nested state", {
  if(!spawn_mpi_slaves()) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(),add=TRUE)
  local <- getFromNamespace(".npRmpi_with_local_regression","npRmpi")
  keys <- c("npRmpi.local.regression.mode","npRmpi.autodispatch.disable",
            "npRmpi.autodispatch.context")
  saved <- setNames(lapply(keys,getOption),keys)
  on.exit(options(saved),add=TRUE)
  for(value in list(NULL,FALSE,TRUE)) {
    options(setNames(rep(list(value),length(keys)),keys))
    before <- b1_context_snapshot()
    expect_error(local({
      expect_true(all(vapply(keys,function(k)isTRUE(getOption(k)),logical(1))))
      expect_identical(mpi.comm.size(1),1L)
      local(stop("local scope sentinel"))
    }),"local scope sentinel")
    expect_identical(b1_context_snapshot(),before)
  }
})
