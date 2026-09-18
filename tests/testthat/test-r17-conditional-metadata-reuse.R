test_that("conditional X metadata adapter retains canonical finalized fields", {
  adapter <- getFromNamespace(".npcdhat_make_xbw", "np")
  specfun <- getFromNamespace("npConditionalRegEngineSpec", "np")
  incumbent <- function(b, x) {
    s <- specfun(b, where = "conditional hat")
    do.call(npregbw, list(xdat=x, ydat=rep.int(0.0,nrow(x)), bws=b$xbw,
      regtype=s$reg.engine, basis=s$basis.engine, degree=s$degree.engine,
      bernstein.basis=s$bernstein.engine, bwtype=b$type, bandwidth.compute=FALSE,
      ckertype=b$cxkertype, ckerorder=b$cxkerorder, ckerbound=b$cxkerbound,
      ckerlb=b$cxkerlb, ckerub=b$cxkerub, ukertype=b$uxkertype, okertype=b$oxkertype))
  }
  capture <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning=function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    list(value=value, warnings=warnings)
  }
  clean <- function(z) {
    z$call <- NULL
    z$ynames <- NULL
    z$varnames$y <- NULL
    z$total.time <- NULL
    z$timing.profile <- NULL
    # Distinct live MPI publications necessarily have distinct lease identities.
    for(nm in c("npRmpi.autodispatch.fingerprint","npRmpi.autodispatch.lease",
                "npRmpi.autodispatch.remote")) attr(z,nm) <- NULL
    z[[".np.native.training"]] <- NULL
    # Retention is attached by the public wrapper, not the numerical finalizer.
    attr(z, ".np.native.training") <- NULL
    z
  }
  set.seed(17129)
  x <- data.frame(u=factor(rep(c("a","b","c"),20), levels=c("a","b","c","d")),
                  x=runif(60,.05,.95), o=ordered(rep(1:3,20)),
                  z=runif(60,.05,.95))
  y <- data.frame(y=runif(60,.05,.95))
  for (bt in c("fixed","generalized_nn","adaptive_nn")) {
    for (ker in c("gaussian","epanechnikov","uniform","beta")) {
      if (ker=="beta" && bt!="fixed") next
      for (rt in c("lc","ll","lp0","lp")) for (bern in c(FALSE,TRUE)) {
        if (rt %in% c("lc","ll") && bern) next
        args <- list(xdat=x,ydat=y,
          bws=c(if(bt=="fixed") .4 else 35,.2,if(bt=="fixed") .4 else 35,
                 .2,if(bt=="fixed") .4 else 35),
          regtype=if(rt=="lp0") "lp" else rt, bwtype=bt, cxkertype=ker,
          cxkerbound=if(ker=="beta") "fixed" else "none",
          cxkerlb=if(ker=="beta") 0 else NULL,
          cxkerub=if(ker=="beta") 1 else NULL, bandwidth.compute=FALSE)
        if(rt %in% c("lp","lp0")) {
          args$degree <- if(rt=="lp0") c(0L,0L) else c(2L,1L)
          args$bernstein.basis <- bern
        }
        b <- suppressWarnings(do.call(npcdensbw,args))
        a <- capture(incumbent(b,x)); z <- capture(adapter(b,x))
        expect_identical(clean(z$value),clean(a$value),
          info=paste(bt,ker,rt,bern))
        expect_identical(z$warnings,a$warnings,info=paste(bt,ker,rt,bern))
      }
    }
  }
})

test_that("conditional adapter bypasses default setup but executes the finalizer", {
  ns <- asNamespace("np")
  adapter <- get(".npcdhat_make_xbw",ns)
  set.seed(17130)
  x <- data.frame(x=runif(30)); y <- data.frame(y=rnorm(30))
  b <- npcdensbw(xdat=x,ydat=y,bws=c(.4,.4),bandwidth.compute=FALSE)
  count <- new.env(parent=emptyenv()); count$default <- 0L; count$final <- 0L
  trace("npregbw.default", where=ns, tracer=substitute({
    assign("default",get("default",envir=COUNTER)+1L,envir=COUNTER)
  },list(COUNTER=count)), print=FALSE)
  on.exit(untrace("npregbw.default",where=ns),add=TRUE)
  trace("npregbw.rbandwidth", where=ns, tracer=substitute({
    assign("final",get("final",envir=COUNTER)+1L,envir=COUNTER)
  },list(COUNTER=count)), print=FALSE)
  on.exit(untrace("npregbw.rbandwidth",where=ns),add=TRUE)
  z <- adapter(b,x)
  expect_identical(count$default,0L)
  expect_gte(count$final,1L)
  expect_identical(z$bw,b$xbw)
})
