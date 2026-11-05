test_that("deterministic D8 dispatch restores tags on local success and unwind", {
  ns<-asNamespace("npRmpi");f<-get(".np_udens_d8_pilot_dispatch",ns)
  e<-new.env(parent=ns);environment(f)<-e
  e$.npRmpi_bootstrap_fanout_enabled<-function(...)TRUE
  e$.npRmpi_bootstrap_run_fanout<-function(...) {
    runif(1)
    if(fail)stop("ordinary local failure")
    matrix(0,1,2)
  }
  old<-RNGkind();on.exit(do.call(RNGkind,as.list(old)),add=TRUE)
  for(present in c(TRUE,FALSE)) for(fail in c(FALSE,TRUE)) {
    RNGkind("Mersenne-Twister","Box-Muller","Rejection");set.seed(5)
    invisible(rnorm(1));expected<-list(normal=rnorm(4),uniform=runif(4))
    set.seed(5);invisible(rnorm(1));seed<-.Random.seed
    if(!present)rm(".Random.seed",envir=.GlobalEnv)
    if(fail)expect_error(f(diag(2),matrix(1,2,1)),"ordinary local failure") else
      expect_identical(f(diag(2),matrix(1,2,1)),matrix(0,2,1))
    expect_identical(exists(".Random.seed",.GlobalEnv,inherits=FALSE),present)
    if(present) {
      expect_identical(.Random.seed,seed)
      expect_identical(list(normal=rnorm(4),uniform=runif(4)),expected)
    }
  }
})

test_that("D8 transport preserves caller uniform and cached-normal continuation", {
  if(TRUE) {
    skip_on_cran();if(!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES","1"))))stop("MPI unavailable")
    on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  }
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  kinds<-RNGkind();on.exit(do.call(RNGkind,as.list(kinds)),add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(611)
  d<-data.frame(x=rnorm(120));f<-npudens(bws=npudensbw(dat=d,bws=.4,bandwidth.compute=FALSE),tdat=d)
  full<-identical(Sys.getenv("N5_FULL"),"true")
  payloads<-list()
  for(kind in if(full)c("Box-Muller","Inversion") else "Box-Muller")
    for(pre in if(full)0:1 else 1L) for(B in if(full)c(257L,288L,289L) else 257L) {
      RNGkind(normal.kind=kind);set.seed(5);if(pre)invisible(rnorm(pre))
      expected<-list(normal=rnorm(4),uniform=runif(4))
      set.seed(5);if(pre)invisible(rnorm(pre));seed<-.Random.seed
      value<-expect_warning(plot(f,errors="asymptotic",center="bias-corrected",band="all",neval=11L,
        gaussian_draws=B,output="data"),"below the recommended")
      expect_identical(.Random.seed,seed)
      expect_identical(list(normal=rnorm(4),uniform=runif(4)),expected)
      payloads[[paste(kind,pre,B)]]<-value
    }
  if(full) {
    RNGkind(normal.kind="Box-Muller");set.seed(5);invisible(rnorm(1));seed<-.Random.seed
    payloads$default <- plot(f,errors="asymptotic",center="bias-corrected",band="all",neval=11L,output="data")
    expect_identical(.Random.seed,seed)
  }
  RNGkind(normal.kind="Box-Muller");set.seed(5);rm(".Random.seed",envir=.GlobalEnv)
  expect_warning(plot(f,errors="asymptotic",center="bias-corrected",band="all",neval=11L,
    gaussian_draws=257L,output="data"),"below the recommended")
  expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
  if(nzchar(Sys.getenv("N5_SNAPSHOT")))saveRDS(payloads,Sys.getenv("N5_SNAPSHOT"))
})
