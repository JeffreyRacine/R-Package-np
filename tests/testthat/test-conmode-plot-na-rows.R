test_that("conmode uncertainty uses complete inputs and restores omitted plot rows", {
  skip_on_cran() # Multiple bootstrap laws and interval/gradient branches.
  old <- options(np.messages=FALSE); on.exit(options(old), add=TRUE)
  set.seed(11); d <- data.frame(x=runif(80))
  d$g <- factor(ifelse(d$x+rnorm(80,sd=.3)>.5,"hi","lo"))
  d$x[17] <- NA; d$g[32] <- NA
  keep <- complete.cases(d); clean <- d[keep, , drop=FALSE]
  count <- 0L
  omit <- function(z) {count <<- count+1L; na.omit(z)}
  fit <- npconmode(g~x,data=d,bws=c(.2,.3),probabilities=TRUE,se=TRUE,
                  gradients=TRUE,na.action=omit)
  expect_identical(count,1L)
  ref <- npconmode(g~x,data=clean,bws=c(.2,.3),probabilities=TRUE,se=TRUE,gradients=TRUE)
  saved <- fit
  # Primitive probability oracle on the retained sample.
  K <- outer(clean$x,clean$x,function(a,b)dnorm((a-b)/.3))
  truth <- as.vector(K %*% ifelse(clean$g==levels(clean$g)[1L],.8,.2))/rowSums(K)
  expect_equal(unname(fit$probabilities[keep,1L]),truth,tolerance=2e-13)
  compare <- function(fit, ref, rows, errors, band="pointwise", gradients=FALSE,
                      bootstrap="inid", view="sample") {
    args <- list(errors=errors,band=band,output="data",gradients=gradients,view=view)
    if(errors=="bootstrap") {
      # Equality/ownership fixture, not a coverage-precision experiment.
      args$B<-7L;args$bootstrap<-bootstrap
      if(bootstrap!="inid")args$boot.control<-np_boot_control(blocklen=3L)
    }
    if(errors=="asymptotic" && band %in% c("all","simultaneous"))args$gaussian.draws<-99L
    if(view=="fixed")args$neval<-9L
    seed <- .Random.seed
    a <- suppressWarnings(do.call(plot,c(list(fit),args)))
    b <- suppressWarnings(do.call(plot,c(list(ref),args)))
    expect_identical(.Random.seed,seed)
    expect_identical(length(a$panels),length(b$panels))
    bands<-if(band=="all")c("pointwise","bonferroni","simultaneous")else band
    for(bd in bands) for(j in seq_along(a$panels)) {
      aa<-as.data.frame(a,panel=j,band=bd);bb<-as.data.frame(b,panel=j,band=bd)
      expect_identical(nrow(aa),length(rows))
      cols<-c("estimate","bias_corrected","std_error","lower","upper")
      for(nm in cols)expect_equal(aa[[nm]][rows],bb[[nm]],tolerance=2e-12,info=paste(errors,bd,view,nm))
      expect_true(all(is.na(aa$lower[!rows])))
      expect_true(all(is.na(aa$upper[!rows])))
    }
    invisible(a)
  }
  for(method in c("inid","fixed","geom"))compare(fit,ref,keep,"bootstrap",bootstrap=method)
  compare(fit,ref,keep,"bootstrap",band="all")
  compare(fit,ref,keep,"asymptotic",band="simultaneous")
  compare(fit,ref,keep,"asymptotic",band="all",gradients=TRUE)
  compare(fit,ref,rep(TRUE,9),"bootstrap",view="fixed")
  compare(fit,ref,rep(TRUE,9),"asymptotic",band="simultaneous",view="fixed")
  native<-npconmode(txdat=d["x"],tydat=d["g"],bws=c(.2,.3),probabilities=TRUE,se=TRUE)
  compare(native,ref,keep,"bootstrap")
  ex<-data.frame(x=c(.1,NA,.6,.8));ek<-complete.cases(ex)
  fe<-npconmode(txdat=d["x"],tydat=d["g"],exdat=ex,bws=c(.2,.3),probabilities=TRUE,se=TRUE)
  re<-npconmode(txdat=clean["x"],tydat=clean["g"],exdat=ex[ek,,drop=FALSE],bws=c(.2,.3),probabilities=TRUE,se=TRUE)
  compare(fe,re,ek,"bootstrap")
  compare(fe,re,ek,"asymptotic",band="simultaneous")
  # The chunked bootstrap owner uses the same compact/padded boundary.
  old.block <- getOption("np.plot.conmode.apply.block.bytes")
  options(np.plot.conmode.apply.block.bytes=8*sum(keep)*3)
  on.exit(options(np.plot.conmode.apply.block.bytes=old.block),add=TRUE)
  chunk.calls <- 0L
  chunk.owner <- getFromNamespace(".np_plot_conmode_bootstrap_values_chunked", "npRmpi")
  local_mocked_bindings(.np_plot_conmode_bootstrap_values_chunked = function(...) {
    chunk.calls <<- chunk.calls + 1L
    chunk.owner(...)
  }, .package = "npRmpi")
  compare(fit,ref,keep,"bootstrap")
  expect_identical(chunk.calls,2L)
  expect_identical(fit,saved)
  expect_identical(count,1L)
})
