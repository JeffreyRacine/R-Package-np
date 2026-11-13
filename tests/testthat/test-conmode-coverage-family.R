test_that("class-probability families exclude only omitted coordinates", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  x<-seq(.02,.98,length.out=40);g<-factor(rep(c("a","b"),20))
  for(p in 1:2) {
    tr<-data.frame(x=x);e<-data.frame(x=c(.1,NA,.3,.5,NA,.7,.9))
    if(p==2){tr$z<-rev(x)^2;e$z<-c(.8,.7,NA,.4,.3,.2,.1)}
    keep<-complete.cases(e);bw<-c(.2,rep(.3,p))
    fit<-npconmode(txdat=tr,tydat=g,exdat=e,bws=bw,probabilities=TRUE,se=TRUE,gradients=TRUE)
    ref<-npconmode(txdat=tr,tydat=g,exdat=e[keep,,drop=FALSE],bws=bw,probabilities=TRUE,se=TRUE,gradients=TRUE)
    saved<-fit
    for(band in c("pointwise","bonferroni","all")) {
      set.seed(72);seed<-.Random.seed
      args<-list(view="sample",errors="asymptotic",band=band,output="data")
      if(band=="all")args$gaussian_draws<-99L
      if(band=="all") {
        expect_warning(a<-do.call(plot,c(list(fit),args)),"below the recommended")
        expect_warning(b<-do.call(plot,c(list(ref),args)),"below the recommended")
      } else {
        a<-do.call(plot,c(list(fit),args));b<-do.call(plot,c(list(ref),args))
      }
      expect_identical(.Random.seed,seed)
      for(j in seq_along(a$panels))for(bd in if(band=="all")c("pointwise","bonferroni","simultaneous")else band){
        aa<-as.data.frame(a,panel=j,band=bd);bb<-as.data.frame(b,panel=j,band=bd)
        for(nm in c("estimate","std_error","lower","upper"))expect_equal(aa[[nm]][keep],bb[[nm]],tolerance=2e-13)
        expect_true(all(is.na(aa$lower[!keep])))
        if(bd=="bonferroni")expect_equal((aa$upper-aa$estimate)[keep]/aa$std_error[keep],rep(qnorm(1-.05/(2*sum(keep))),sum(keep)),tolerance=2e-13)
      }
      for(fam in a$coverage_families)expect_true(all(fam$members$row_id %in% which(keep)))
      for(panel in a$panels)expect_true(all(is.na(panel$intervals$family_id[!keep[panel$intervals$row_id]])))
    }
    # Gradient arithmetic already has the complete-row family; metadata must agree.
    a<-plot(fit,view="sample",gradients=TRUE,errors="asymptotic",band="bonferroni",output="data")
    b<-plot(ref,view="sample",gradients=TRUE,errors="asymptotic",band="bonferroni",output="data")
    for(j in seq_along(a$panels))expect_equal(as.data.frame(a,panel=j)$upper[keep],as.data.frame(b,panel=j)$upper,tolerance=2e-13)
    for(fam in a$coverage_families)expect_true(all(fam$members$row_id %in% which(keep)))
    expect_identical(fit,saved)
  }
})

test_that("family exclusions do not discard defined but unavailable uncertainty", {
  ctor<-getFromNamespace(".np_plot_data_family_events","np")
  x<-list(panels=list(panel_1=list(intervals=data.frame(
    row_id=1:3,target_id="t",band="bonferroni",alpha=.05,
    lower=c(0,NA,NA),upper=c(1,NA,NA),status=c("available","undefined","undefined"),family_id="f"))),
    coverage_families=list(f=list(scope="test scope")))
  unchanged<-ctor(x);explicit<-ctor(x,list(panel_1=2L))
  expect_identical(unchanged$coverage_families[[1]]$members$row_id,1:3)
  expect_identical(explicit$coverage_families[[1]]$members$row_id,c(1L,3L))
  expect_identical(explicit$panels$panel_1$intervals$status,x$panels$panel_1$intervals$status)
  expect_true(is.na(explicit$panels$panel_1$intervals$family_id[2]))
})
