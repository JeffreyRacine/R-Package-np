test_that("beta CVLS search admission is lazy and restricted to its policy cone", {
  admit <- getFromNamespace(".np_beta_cvls_search_admission", "np")
  b <- list(type="fixed", method="cv.ls", ckertype="beta", icon=c(FALSE,TRUE),
            ckerlb=c(-Inf,0), ckerub=c(Inf,1))
  y <- data.frame(group=factor(c("a","a","b","b")), y=c(0,0,.5,1))
  expect_error(admit(b,y,0,TRUE), "positive 'scale.factor.search.lower'", fixed=TRUE)
  expect_null(admit(b,stop("data must not be evaluated"),.1,TRUE))
  expect_null(admit(b,stop("data must not be evaluated"),0,FALSE))
  for (type in c("generalized_nn","adaptive_nn")) {
    control <- b; control$type <- type
    expect_null(admit(control,stop("data must not be evaluated"),0,TRUE))
  }
  for (method in c("cv.ml","normal-reference")) {
    control <- b; control$method <- method
    expect_null(admit(control,stop("data must not be evaluated"),0,TRUE))
  }
  control <- b; control$ckertype <- "gaussian"
  expect_null(admit(control,stop("data must not be evaluated"),0,TRUE))
  expect_null(admit(b,y[c(1,3,4),],0,TRUE))
  upper <- y; upper$y <- c(0,.5,1,1)
  expect_error(admit(b,upper,0,TRUE), "repeated response endpoints", fixed=TRUE)
  interior <- y; interior$y <- c(.1,.1,.5,.9)
  expect_null(admit(b,interior,0,TRUE))
  translated <- b; translated$ckerlb <- c(-Inf,-2); translated$ckerub <- c(Inf,5)
  translated.y <- y; translated.y$y <- -2+7*y$y
  expect_error(admit(translated,translated.y,0,TRUE),"repeated response endpoints",fixed=TRUE)
  conditional <- list(type="fixed",method="cv.ls",cykertype="beta",
    iycon=b$icon,cykerlb=b$ckerlb,cykerub=b$ckerub)
  expect_error(admit(conditional,y,0,TRUE,TRUE,data.frame(x=1:4)),
               "repeated response endpoints",fixed=TRUE)
  expect_null(admit(conditional,y,0,TRUE,TRUE,data.frame(x=c(1,NA,3,4))))
  missing <- y; missing$group[2] <- NA
  expect_null(admit(b,missing,0,TRUE))
})

test_that("automatic beta density search owners enforce the same endpoint policy", {
  skip_on_cran()
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  y <- data.frame(y=c(0,0,seq(.1,.9,length.out=12),1,1))
  x <- data.frame(x=seq(-1,1,length.out=nrow(y)))
  u <- list(dat=y, bws=.2, bwmethod="cv.ls", ckertype="beta",
            ckerbound="fixed", ckerlb=0, ckerub=1, scale.factor.search.lower=0)
  cargs <- list(xdat=x, ydat=y, bws=c(.2,.4), bwmethod="cv.ls",
                cykertype="beta", cykerbound="fixed", cykerlb=0, cykerub=1,
                scale.factor.search.lower=0)
  message <- "positive 'scale.factor.search.lower'"
  ub <- do.call(npudensbw,c(u,list(bandwidth.compute=FALSE)))
  cb <- do.call(npcdensbw,c(cargs,list(bandwidth.compute=FALSE)))
  expect_s3_class(ub,"bandwidth")
  expect_s3_class(cb,"conbandwidth")
  for (solver in c("powell","mads","mads+powell")) {
    expect_error(do.call(npudensbw,c(u,list(bwsolver=solver))),message,fixed=TRUE)
    expect_error(do.call(npcdensbw,c(cargs,list(bwsolver=solver))),message,fixed=TRUE)
    expect_error(npudensbw(bws=ub,dat=y,bwsolver=solver),message,fixed=TRUE)
    expect_error(npcdensbw(bws=cb,xdat=x,ydat=y,bwsolver=solver),message,fixed=TRUE)
  }
  for (regtype in c("ll","lp")) {
    extra <- if(regtype=="lp") list(degree=2) else list()
    expect_error(do.call(npcdensbw,c(cargs,list(regtype=regtype),extra)),message,fixed=TRUE)
  }
  expect_error(do.call(npcdensbw,c(cargs,list(nomad=TRUE))),message,fixed=TRUE)
  expect_error(npudensbw(~y,data=y,bws=.2,bwmethod="cv.ls",ckertype="beta",
    ckerbound="fixed",ckerlb=0,ckerub=1,scale.factor.search.lower=0),message,fixed=TRUE)
  d <- cbind(x,y)
  expect_error(npcdensbw(y~x,data=d,bws=c(.2,.4),bwmethod="cv.ls",
    cykertype="beta",cykerbound="fixed",cykerlb=0,cykerub=1,
    scale.factor.search.lower=0),message,fixed=TRUE)
  ue <- getFromNamespace("npudensbw.bandwidth","np")(dat=y,bws=ub,
    bandwidth.compute=TRUE,eval.only=TRUE,nmulti=1)
  ce <- getFromNamespace(".npcdensbw_eval_only","np")(xdat=x,ydat=y,bws=cb)
  expect_true(is.finite(ue$fval))
  expect_true(is.finite(ce$objective))
  expect_length(fitted(npudens(bws=ub,tdat=y)),nrow(y))
  expect_length(fitted(npcdens(bws=cb,txdat=x,tydat=y)),nrow(y))
})
