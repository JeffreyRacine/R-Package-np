test_that("proper certificates reflect kernel signs and response geometry", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920351)
  x <- data.frame(x=runif(24,.05,.95))
  y <- data.frame(y=runif(24,.05,.95))
  ex <- data.frame(x=rep(.5,31L))
  for (family in c("npcdens","npcdist")) {
    for (kernel in c("gaussian","epanechnikov","uniform","beta")) {
      for (type in c("fixed","generalized_nn","adaptive_nn")) {
        bounds <- if (kernel=="beta") "fixed" else "none"
        args <- list(txdat=x,tydat=y,exdat=ex,
          eydat=data.frame(y=seq(0,1,length.out=31L)),
          bws=if(type=="fixed") c(.2,.25) else c(8,9),bwtype=type,
          cykertype=kernel,cykerbound=bounds,cykerlb=if(bounds=="fixed") 0 else NULL,
          cykerub=if(bounds=="fixed") 1 else NULL,se=FALSE)
        ctrl <- list(fail.on.unsupported=TRUE)
        if (family=="npcdens") ctrl$mass.warn.tol <- 0
        fit <- do.call(get(family),c(args,list(proper=TRUE,proper.control=ctrl)))
        certified <- type!="generalized_nn" &&
          (family=="npcdist" || kernel!="beta")
        expect_identical(fit$proper.info$reason=="already_proper",certified)
        expect_identical(fit$proper.applied,!certified)
        raw <- do.call(get(family),args)
        if (certified)
          expect_identical(fitted(fit),fitted(raw))
        if (family=="npcdist") {
          expect_true(all(fitted(fit)>=0 & fitted(fit)<=1))
          expect_true(all(diff(fitted(fit))>=-1e-12))
        } else if (!certified) {
          weights <- c(.5,rep(1,29),.5)/30
          expect_true(all(fitted(fit)>=0))
          expect_lt(abs(sum(weights*fitted(fit))-1),1e-9)
        }
      }
    }
    # A signed predictor kernel is independently disqualifying.
    for (side in c("x","y")) for (order in c(4L,6L,8L)) {
      a <- list(txdat=x,tydat=y,exdat=ex,eydat=data.frame(y=seq(-2,3,length.out=31)),
                bws=c(.2,.25),se=FALSE,proper=TRUE,proper.control=ctrl)
      a[[paste0("c",side,"kerorder")]] <- order
      fit <- do.call(get(family),a)
      expect_true(fit$proper.applied)
      expect_false(fit$proper.info$reason=="already_proper")
    }
    # Finite Y bounds disqualify target-centred densities, not CDF mixtures.
    for (bound in list(c(0,1),c(0,Inf),c(-Inf,1))) {
      fit <- do.call(get(family),list(txdat=x,tydat=y,exdat=ex,
        eydat=data.frame(y=seq(0,1,length.out=31)),bws=c(.2,.25),
        cykerbound="fixed",cykerlb=bound[1],cykerub=bound[2],proper=TRUE,
        proper.control=ctrl))
      expect_identical(fit$proper.info$reason=="already_proper",family=="npcdist")
    }
  }
})

test_that("declared categorical response normalization controls certification", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(.1,.9,length.out=30))
  for (ordered in c(FALSE,TRUE)) {
    y <- data.frame(y=factor(rep(c("a","b","c"),10),ordered=ordered))
    ey <- data.frame(y=factor(c("a","b","c"),levels=levels(y$y),ordered=ordered))
    args <- list(txdat=x,tydat=y,exdat=data.frame(x=rep(.5,3)),eydat=ey,
                  bws=c(.1,.3),proper=TRUE)
    args[[if(ordered) "oykertype" else "uykertype"]] <-
      if(ordered) "racineliyan" else "aitchisonaitken"
    fit <- do.call(npcdens,args)
    expect_identical(fit$proper.info$reason,"already_proper")
    expect_lt(abs(sum(fitted(fit))-1),1e-12)
    args[[if(ordered) "oykertype" else "uykertype"]] <- "liracine"
    unnormalized <- do.call(npcdens,args)
    expect_false(unnormalized$proper.info$supported)
    expect_false(unnormalized$proper.info$reason=="already_proper")
    args$proper.control <- list(fail.on.unsupported=TRUE)
    expect_error(do.call(npcdens,args),"univariate continuous")
  }
})
