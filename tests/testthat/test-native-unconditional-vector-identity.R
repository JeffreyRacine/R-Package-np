test_that("observed native vector roles survive unconditional fits and plots", {
  if (environmentName(environment(npudist)) == "npRmpi") {
    skip_on_cran()
    if (!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES", "1")))) stop("MPI unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages=FALSE); on.exit(options(old), add=TRUE)
  set.seed(4121); original <- rnorm(60)
  for (family in c("npudist", "npudens")) {
    bwfun <- get(paste0(family,"bw")); fitfun <- get(family)
    bw <- bwfun(dat=original, bws=.4, bandwidth.compute=FALSE)
    named <- setNames(data.frame(original),bw$xnames)
    bwref <- bwfun(dat=named,bws=.4,bandwidth.compute=FALSE)
    reference <- fitfun(bws=bwref,tdat=named)
    observed <- fitfun(bws=bw,tdat=original)
    expect_identical(observed$bws[[".np.native.training"]]$dat,named)
    expect_equal(unname(fitted(observed)),unname(fitted(reference)),tolerance=0)
    expected <- plot(reference,errors="asymptotic",center="bias-corrected",band="pointwise",neval=6L,output="data")
    frozen <- serialize(observed,NULL)
    for (object in list(observed,unserialize(frozen),bw)) {
      actual <- plot(object,errors="asymptotic",center="bias-corrected",band="pointwise",neval=6L,output="data")
      expect_identical(actual$panels,expected$panels)
    }
    direct <- plot(bw,xdat=original,errors="asymptotic",center="bias-corrected",band="pointwise",neval=6L,output="data")
    expect_identical(direct$panels,expected$panels)
    expect_identical(serialize(observed,NULL),frozen)
    expect_identical(observed$bws[[".np.native.training"]]$dat,named)
    altered <- named; names(altered) <- "genuinely.different"
    wrong <- fitfun(bws=bw,tdat=altered)
    expect_identical(wrong$bws[[".np.native.training"]]$dat,altered)
    if (family=="npudist") expect_error(plot(wrong,errors="asymptotic",center="bias-corrected",neval=6L,output="data"),"invalid CDF correction training frame")
    # A saved object with an ambiguous old frame is not retroactively renamed.
    legacy <- observed; names(legacy$bws[[".np.native.training"]]$dat) <- "tdat"
    if (family=="npudist") expect_error(plot(legacy,errors="asymptotic",center="bias-corrected",neval=6L,output="data"),"invalid CDF correction training frame")
    if (identical(Sys.getenv("N5_FULL"),"true")) {
      args <- list(errors="bootstrap",bootstrap="inid",B=999L,random_seed=14L,center="bias-corrected",band="all",neval=6L,output="data")
      expect_identical(do.call(plot,c(list(observed),args))$panels,do.call(plot,c(list(reference),args))$panels)
    }
  }
})

test_that("native frame identity is never inferred from values or old names", {
  if (environmentName(environment(npudist)) == "npRmpi") {
    skip_on_cran()
    if (!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES", "1")))) stop("MPI unavailable")
    on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  }
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(4122); x <- rnorm(60)
  for (family in c("npudist","npudens")) {
    bwfun <- get(paste0(family,"bw")); fitfun <- get(family)
    frame <- data.frame(`a b`=x,check.names=FALSE)
    bw <- bwfun(dat=frame,bws=.4,bandwidth.compute=FALSE)
    for (given in list(frame,setNames(frame,"wrong"),as.matrix(frame))) {
      fit <- fitfun(bws=bw,tdat=given)
      expect_identical(fit$bws[[".np.native.training"]]$dat,as.data.frame(given))
    }
  }
  bw <- npudistbw(~x,bws=.4,bandwidth.compute=FALSE)
  fit <- npudist(bw)
  expected <- plot(fit,errors="asymptotic",center="bias-corrected",neval=6L,band="pointwise",output="data")
  x[] <- 999
  expect_identical(plot(fit,errors="asymptotic",center="bias-corrected",neval=6L,band="pointwise",output="data")$panels,expected$panels)
})
