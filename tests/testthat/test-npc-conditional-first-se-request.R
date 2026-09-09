l1_conditional_fixture <- function(cdf = FALSE) {
  i <- seq_len(31L)
  x <- data.frame(x = .03+.94*(i/32)^1.2,
                  z = .04+.92*((i*7L) %% 37L)/37)
  y <- data.frame(y = .03+.94*((i*11L) %% 37L)/37)
  e <- c(3L,8L,15L,22L,29L)
  args <- list(xdat=x,ydat=y,bws=c(.2,.3,.35),bandwidth.compute=FALSE,
    bwscaling=FALSE,regtype="lp",degree=c(1L,0L),basis="additive")
  b <- do.call(if(cdf)npcdistbw else npcdensbw,args)
  list(x=x,y=y,ex=x[e,,drop=FALSE],ey=y[e,,drop=FALSE],b=b)
}

test_that("partial conditional LP restores only requested first derivative errors", {
  if(!spawn_mpi_slaves(1L)) skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  for(cdf in c(FALSE,TRUE)) {
    d <- l1_conditional_fixture(cdf)
    fun <- if(cdf)npcdist else npcdens
    a <- list(bws=d$b,txdat=d$x,tydat=d$y,exdat=d$ex,eydat=d$ey,gradients=TRUE)
    none <- suppressWarnings(do.call(fun,c(a,list(.np_lp_first_se_demand=FALSE))))
    fit <- suppressWarnings(do.call(fun,a))
    expected <- vapply(seq_len(nrow(d$ex)),function(j) {
      w <- dnorm((d$x$x-d$ex$x[j])/.3)*dnorm((d$x$z-d$ex$z[j])/.35)
      z <- if(cdf) pnorm((d$ey$y[j]-d$y$y)/.2) else
        dnorm((d$ey$y[j]-d$y$y)/.2)/.2
      p <- cbind(1,d$x$x-d$ex$x[j])
      inv <- solve(crossprod(p,w*p))
      sigma2 <- sum(w*z^2)/sum(w)-(sum(w*z)/sum(w))^2
      v <- sigma2*inv %*% crossprod(p,w^2*p) %*% inv
      sqrt(max(0,v[2L,2L]))
    },0.0)
    expect_identical(fitted(fit),fitted(none))
    expect_identical(se(fit),se(none))
    expect_identical(fit$congrad,none$congrad)
    expect_identical(fit$gradient.order,none$gradient.order)
    expect_true(all(is.na(none$congerr)))
    expect_true(all(is.na(fit$congerr[,2L])))
    expect_equal(unname(fit$congerr[,1L]),expected,tolerance=2e-10)
    expect_identical(gradients(fit,se=TRUE),fit$congerr)
    suppressed <- suppressWarnings(do.call(fun,c(a,list(
      .np_lp_first_se_demand=c(FALSE,TRUE)))))
    expect_identical(suppressed$congerr,none$congerr)
  }
})

test_that("private conditional first-SE demand keeps exact dots validation", {
  if(!spawn_mpi_slaves(1L)) skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  d <- l1_conditional_fixture()
  expect_error(npcdens(bws=d$b,txdat=d$x,tydat=d$y,
    .np_lp_first_se_deman=stop("must not evaluate")),
    "unused|unrecognized|unknown")
  expect_error(npcdens(bws=d$b,txdat=d$x,tydat=d$y,
    .np_lp_first_se_demand=TRUE),"invalid internal conditional first-SE demand")
})

test_that("partial conditional LP masks two directions in original column order", {
  if(!spawn_mpi_slaves(1L)) skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  d <- l1_conditional_fixture(TRUE)
  x <- data.frame(x=d$x$x,cat=factor(rep(c("a","b","c"),length.out=nrow(d$x))),
    z=d$x$z,v=.1+.8*((seq_len(nrow(d$x))*13L) %% 41L)/41)
  e <- c(3L,8L,15L,22L,29L)
  b <- npcdistbw(xdat=x,ydat=d$y,bws=c(.2,.3,.2,.35,.31),
    bandwidth.compute=FALSE,bwscaling=FALSE,regtype="lp",
    degree=c(1L,0L,1L),basis="additive")
  args <- list(bws=b,txdat=x,tydat=d$y,exdat=x[e,,drop=FALSE],
    eydat=d$ey,gradients=TRUE)
  full <- suppressWarnings(do.call(npcdist,args))
  subset <- suppressWarnings(do.call(npcdist,c(args,list(
    .np_lp_first_se_demand=c(FALSE,FALSE,TRUE)))))
  expect_true(all(is.finite(full$congerr[,c(1L,4L)])))
  expect_true(all(is.na(full$congerr[,c(2L,3L)])))
  expect_true(all(is.na(subset$congerr[,1L:3L])))
  expect_identical(full$congerr[,4L],subset$congerr[,4L])
  expect_identical(fitted(full),fitted(subset))
  expect_identical(se(full),se(subset))
  expect_identical(full$congrad,subset$congrad)
  args$eydat$y <- 10
  zero <- suppressWarnings(do.call(npcdist,args))
  expect_identical(unname(zero$congerr[,c(1L,4L)]),matrix(0,nrow(d$ey),2L))
  expect_true(all(is.na(zero$congerr[,c(2L,3L)])))
})
