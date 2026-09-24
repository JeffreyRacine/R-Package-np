test_that("smooth-coefficient raw CV certifies the deleted endpoint", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 24L
  x <- data.frame(x=seq_len(n)/n)
  z <- data.frame(z=c(rep(0,22L),1:2))
  y <- 1+x$x+sin(seq_len(n))/10
  ctx <- .npscoefbw_nomad_context_prepare(x,y,z)
  for(type in c("generalized_nn","adaptive_nn")) {
    b <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=22,bwtype=type,
                  regtype="lc",bandwidth.compute=FALSE)
    raw <- .npscoefbw_nomad_eval_direct(ctx,b)
    literal <- vapply(seq_len(n),function(i) {
      bb <- npscoefbw(xdat=x[-i,,drop=FALSE],ydat=y[-i],zdat=z[-i,,drop=FALSE],
        bws=22,bwtype=type,regtype="lc",bandwidth.compute=FALSE)
      as.numeric(fitted(npscoef(bws=bb,exdat=x[i,,drop=FALSE],ezdat=z[i,,drop=FALSE])))
    },0.0)
    expect_true(raw$raw.valid)
    expect_equal(raw$objective,mean((y-literal)^2),tolerance=2e-12)
    for(k in c(3,12,20,21,23)) {
      bb <- .npscoef_apply_bw_to_scbw(b,k)
      expect_false(.npscoefbw_nomad_eval_direct(ctx,bb)$raw.valid)
    }
  }
})

test_that("automatic smooth-coefficient recovery reaches the deleted endpoint", {
  skip_on_cran()
  skip_if_not_installed("crs",minimum_version="0.15.46")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 24L
  x <- data.frame(x=seq_len(n)/n)
  z <- data.frame(z=c(rep(0,22L),1:2))
  y <- 1+x$x+sin(seq_len(n))/10
  ctx <- .npscoefbw_nomad_context_prepare(x,y,z)
  for(type in c("generalized_nn","adaptive_nn"))
    for(engine in c("R","nomad","nomad+powell")) {
      a <- list(xdat=x,ydat=y,zdat=z,bwtype=type,nmulti=1L,
                optim.maxit=20L,optim.maxattempts=1L,powell.remin=FALSE)
      if(engine=="R") a$regtype <- "lc" else {
        a$nomad <- TRUE;a$search.engine <- engine;a$degree.max <- 1L
        a$nomad.opts <- list(MAX_BB_EVAL=1L)
      }
      set.seed(42)
      b <- do.call(npscoefbw,a)
      expect_identical(as.numeric(b$bw),22)
      raw <- .npscoefbw_nomad_eval_direct(ctx,b)
      expect_true(raw$raw.valid)
      expect_equal(as.numeric(b$fval),raw$objective,tolerance=2e-12)
      set.seed(42)
      a$bws <- 3
      expect_error(do.call(npscoefbw,a),if(engine=="R")
        "no feasible bandwidths found" else "did not return a raw-valid solution")
    }
})

test_that("smooth-coefficient recovery leaves explicit extended counts alone", {
  skip_on_cran()
  skip_if_not_installed("crs",minimum_version="0.15.46")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=TRUE)
  on.exit(options(old),add=TRUE)
  n <- 24L
  x <- data.frame(x=seq_len(n)/n)
  z <- data.frame(z=c(rep(0,22L),1:2))
  y <- 1+x$x+sin(seq_len(n))/10
  ctx <- .npscoefbw_nomad_context_prepare(x,y,z)
  for(type in c("generalized_nn","adaptive_nn")) {
    set.seed(42)
    b <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=26,bwtype=type,
      nomad=TRUE,search.engine="nomad",degree.max=1L,nmulti=1L,
      nomad.opts=list(MAX_BB_EVAL=1L))
    expect_identical(as.numeric(b$bw),26)
    expect_length(b$degree.search$restart.starts,1L)
    raw <- .npscoefbw_nomad_eval_direct(ctx,b)
    expect_true(raw$raw.valid)
    expect_equal(as.numeric(b$fval),raw$objective,tolerance=2e-12)
  }
})
