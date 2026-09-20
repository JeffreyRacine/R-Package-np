test_that("mixed beta NN operators carry the categorical compression policy", {
  old <- options(np.messages=FALSE,
                 np.categorical.compress=getOption("np.categorical.compress"))
  on.exit(options(old),add=TRUE)
  set.seed(92035)
  n <- 25L
  x <- data.frame(x=runif(n,.02,.98),g=factor(rep(letters[1:3],length.out=n)))
  ex <- x[c(2,7,15),,drop=FALSE]
  ex$x <- c(.04,.41,.94)
  y <- data.frame(y=runif(n,.02,.98))
  ey <- data.frame(y=c(.1,.5,.9))
  for (compress in c(FALSE,TRUE)) {
    options(np.categorical.compress=compress)
    for (type in c("fixed","generalized_nn","adaptive_nn")) {
      h <- if (type=="fixed") .2 else 8
      for (family in c("npudens","npudist")) {
        xx <- x
        ee <- ex
        if (family=="npudist") {
          xx$g <- ordered(xx$g)
          ee$g <- ordered(ee$g,levels=levels(xx$g))
        }
        bw <- do.call(get(paste0(family,"bw")),list(dat=xx,bws=c(h,.15),
          bandwidth.compute=FALSE,bwtype=type,ckertype="beta",
          ckerbound="fixed",ckerlb=0,ckerub=1))
        hat <- get(paste0(family,"hat"))
        H <- hat(bw,tdat=xx,edat=ee,output="matrix")
        rhs <- cbind(rep(1,n),seq_len(n)/n)
        expect_equal(hat(bw,tdat=xx,edat=ee,y=rhs,output="apply"),
                     H %*% rhs,tolerance=1e-12,ignore_attr=TRUE)
        fit <- do.call(get(family),list(bws=bw,tdat=xx,edat=ee))
        expect_equal(rowSums(H),as.double(fitted(fit)),tolerance=1e-11)
      }
      for (family in c("npcdens","npcdist")) {
        bw <- do.call(get(paste0(family,"bw")),list(xdat=x,ydat=y,
          bws=c(h,h,.15),bandwidth.compute=FALSE,bwtype=type,regtype="lc",
          cxkertype="beta",cykertype="beta",cxkerbound="fixed",cxkerlb=0,
          cxkerub=1,cykerbound="fixed",cykerlb=0,cykerub=1))
        args <- list(bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey)
        H <- do.call(get(paste0(family,"hat")),c(args,list(output="matrix")))
        fit <- do.call(get(family),args)
        expect_equal(rowSums(H),as.double(fitted(fit)),tolerance=1e-11)
      }
    }
  }
})
