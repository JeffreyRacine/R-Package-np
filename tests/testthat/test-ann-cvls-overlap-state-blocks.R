ann_overlap_reference <- function(x, y, kx, ky, degree, epan=FALSE) {
  n <- length(y)
  rows <- numeric(n)
  integrated <- numeric(n)
  for(i in seq_len(n)) {
    donor <- setdiff(seq_len(n), i)
    hx <- vapply(donor, function(j) sort(abs(x[-c(i,j)]-x[j]), partial=kx)[kx], 0)
    hy <- vapply(donor, function(j) sort(abs(y[-c(i,j)]-y[j]), partial=ky)[ky], 0)
    u <- (x[i]-x[donor])/hx
    w <- if(epan) (3/(4*sqrt(5)))*pmax(0,1-u^2/5)/hx else dnorm(u)/hx
    Z <- vapply(0:degree, function(p) (x[donor]-x[i])^p, numeric(n-1L))
    coefficient <- solve(crossprod(Z,w*Z),c(1,rep(0,degree)))
    influence <- drop(w*(Z %*% coefficient))
    h <- sqrt(outer(hy^2,hy^2,"+"))
    overlap <- dnorm(outer(y[donor],y[donor],"-")/h)/h
    I1 <- sum(outer(influence,influence)*overlap)
    integrated[i] <- I1
    I2 <- sum(influence*dnorm((y[i]-y[donor])/hy)/hy)
    rows[i] <- 2*I2-I1
  }
  # The incumbent analytical overlap has a rounded Gaussian constant. Account
  # for that known arithmetic separately from the independent literal formula;
  # do not loosen the comparison or change the shared production primitive.
  list(exact=mean(rows),
       incumbent=mean(rows)+mean(integrated)*(1-0.3989422803/dnorm(0)))
}

test_that("ANN overlap states preserve signed influence and tile boundaries", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE,
                 np.macMseries.accelerate=FALSE)
  on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  cases <- rbind(expand.grid(n=c(17L,129L),degree=0:2,epan=FALSE),
                 data.frame(n=129L,degree=c(0L,2L),epan=TRUE))
  for(case in seq_len(nrow(cases))) {
    n <- cases$n[case]; degree <- cases$degree[case]; epan <- cases$epan[case]
    set.seed(4700+n)
    x <- runif(n,-1,1)
    y <- rnorm(n)
    kx <- as.integer(max(8,floor(n/3)))
    ky <- as.integer(max(4,floor(n/10)))
    expected <- ann_overlap_reference(x,y,kx,ky,degree,epan)
    args <- list(xdat=data.frame(x=x),ydat=data.frame(y=y),bws=c(ky,kx),
                 bwtype="adaptive_nn",bwmethod="cv.ls",regtype=if(degree==0L)"lc"else"lp",
                 cxkertype=if(epan)"epanechnikov"else"gaussian",
                 bandwidth.compute=FALSE)
    if(degree>0L) {
      args$degree <- degree
      args$degree.select <- "manual"
      args$bernstein.basis <- TRUE
    }
    b <- do.call(npcdensbw,args)
    score <- function(xx,yy,bb) {
      call.args <- list(xdat=xx,ydat=yy,bws=bb,invalid.penalty="dbmax",force.local=FALSE)
      command <- substitute({
        options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE,
                np.macMseries.accelerate=FALSE)
        do.call(get(".npcdensbw_eval_only",asNamespace("npRmpi")),ARGS)
      },list(ARGS=call.args))
      getFromNamespace(".npRmpi_bcast_cmd_expr","npRmpi")(
        command,comm=1L,caller.execute=TRUE)$objective
    }
    observed <- score(args$xdat,args$ydat,b)
    expect_lte(abs(observed-expected$incumbent),
               1e-10+1e-10*abs(expected$incumbent))
    p <- sample.int(n)
    expect_lte(abs(score(args$xdat[p,,drop=FALSE],args$ydat[p,,drop=FALSE],b)-
                   observed),1e-10+1e-10*abs(observed))
  }
})

test_that("ANN tree blocks preserve literal signed influence and original identities", {
  old <- options(np.messages=FALSE,np.tree=TRUE,np.extendednn=FALSE,
                 np.macMseries.accelerate=FALSE)
  on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  cases <- rbind(expand.grid(n=c(17L,129L),degree=0:2,epan=FALSE),
                 data.frame(n=129L,degree=c(0L,2L),epan=TRUE))
  for(case in seq_len(nrow(cases))) {
    n <- cases$n[case]; degree <- cases$degree[case]; epan <- cases$epan[case]
    set.seed(4700+n)
    x <- runif(n,-1,1)
    y <- rnorm(n)
    kx <- as.integer(max(8,floor(n/3)))
    ky <- as.integer(max(4,floor(n/10)))
    expected <- ann_overlap_reference(x,y,kx,ky,degree,epan)
    args <- list(xdat=data.frame(x=x),ydat=data.frame(y=y),bws=c(ky,kx),
                 bwtype="adaptive_nn",bwmethod="cv.ls",regtype=if(degree==0L)"lc"else"lp",
                 cxkertype=if(epan)"epanechnikov"else"gaussian",
                 bandwidth.compute=FALSE)
    if(degree>0L) {
      args$degree <- degree
      args$degree.select <- "manual"
      args$bernstein.basis <- TRUE
    }
    b <- do.call(npcdensbw,args)
    score <- function(xx,yy,bb) {
      call.args <- list(xdat=xx,ydat=yy,bws=bb,invalid.penalty="dbmax",force.local=FALSE)
      command <- substitute({
        options(np.messages=FALSE,np.tree=TRUE,np.extendednn=FALSE,
                np.macMseries.accelerate=FALSE)
        do.call(get(".npcdensbw_eval_only",asNamespace("npRmpi")),ARGS)
      },list(ARGS=call.args))
      getFromNamespace(".npRmpi_bcast_cmd_expr","npRmpi")(
        command,comm=1L,caller.execute=TRUE)$objective
    }
    observed <- score(args$xdat,args$ydat,b)
    expect_lte(abs(observed-expected$incumbent),
               1e-10+1e-10*abs(expected$incumbent))
    p <- sample.int(n)
    expect_lte(abs(score(args$xdat[p,,drop=FALSE],args$ydat[p,,drop=FALSE],b)-
                   observed),1e-10+1e-10*abs(observed))
  }
})
